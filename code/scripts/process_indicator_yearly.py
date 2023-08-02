"""
Process W@H data from zip archive into parquet files of yearly data.


Requires:
---------
./config.json
tempdir/years
datadir/.../w@h/<scenario>.zip

Outputs:
--------
tempdir/<scenario>/yearly/wah_*.parquet
"""
import os
os.environ['USE_PYGEOS'] = '0'
import utils
import pandas as pd
import geopandas as gpd
import xarray as xr
import zipfile
import warnings
from tqdm import tqdm

# set scenario using command likne
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scenario', type=str, help='what scenario to process', default='Hist', choices=['Hist', 'BS', 'NF', 'FF', 'Test'])
parser.add_argument('-v', '--var', type=str, help='what indicator variable to process', default='rainfall', choices=['rainfall', 'ep'])
args = parser.parse_args()
SCENARIO= args.scenario
INDICATOR = args.var
print(f'\nProcessing indicator data into years for {SCENARIO.lower()} and {INDICATOR}')


def condition(filename, year):
    if SCENARIO == "Hist":
        return condition_historical(filename, year)
    if SCENARIO in ["BS", "NF", "FF", "Test"]:
        return condition_scenario(filename, year)
    else:
        raise Exception("Scenario not in ['Hist', 'BS', 'NF', 'FF', 'Test']")

def condition_historical(filename, year):
    if (filename.split('.')[-1] == 'nc') and (year in filename):
        return True
    else:
        return False

def condition_scenario(filename, year):
    if (filename.split('.')[-1] == 'csv') and (year in filename):
        return True
    else:
        return False


def process_wah_file(filename, wah_archive):
    if SCENARIO == "Hist":
        return process_wah_file_historical(filename, wah_archive)
    if SCENARIO in ["BS", "NF", "FF", "Test"]:
        return process_wah_file_scenarios(filename, wah_archive)
    else:
        raise Exception("Scenario not in ['Hist', 'BS', 'NF', 'FF', 'Test']")


def process_wah_file_historical(filename, wah_archive):
    """Process historical data stored in NetCDF format."""
    file = wah_archive.open(filename)
    df = xr.open_dataset(file)
    df = df.to_dataframe().reset_index()
    df['lat'] = df['y'].round(2)
    df['lon'] = df['x'].round(2)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    df['Month'] = df['time'].dt.month
    df['ensemble'] = [f'{SCENARIO}1'] * len(df)
    df = df.rename(columns={f'{INDICATOR}_amount': INDICATOR})
    df_grouped = df.groupby(['lon','lat','Year','Month', 'ensemble'])[INDICATOR].sum().to_frame().reset_index()
    return df_grouped, "EPSG:27700"  # might need to change this, this is for old data

def process_wah_file_scenarios(filename, wah_archive):
    """Process scenario data stored in CSV format."""
    file = wah_archive.open(filename)
    df = pd.read_csv(file)
    df =df.reset_index()
    df['lat'] = df['lat'].round(2)
    df['lon'] = df['lon'].round(2)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    df['Month'] = df['time'].dt.month
    df['ep'] = 86400 * df['prbc'] - df['pepm']
    df[INDICATOR] = -df[INDICATOR]  # need to be precip - evapotranspiration, not reverse
    df_grouped = df.groupby(['lon','lat','Year','Month', 'ensemble'])['ep'].sum().to_frame().reset_index()
    return df_grouped, "EPSG:4326"


def main(config):
    datadir = config['paths']["datadir"]
    outdir = config['paths']["tempdir"]
    uk_crs = config["config"]["uk_crs"]
    YEARS = utils.load_list(os.path.join(outdir, SCENARIO.lower(), "years"))
    wahpath = os.path.join(datadir, 'w@h', f"{SCENARIO.lower()}.zip")

    with zipfile.ZipFile(wahpath, 'r') as archive:
        for year in (pbar:= tqdm(YEARS)):
            pbar.set_description(f"Processing {year}")
            savedir = os.path.join(outdir, SCENARIO.lower(), 'yearly')
            if not os.path.exists(savedir):
                os.makedirs(savedir)
            savepath = os.path.join(savedir, f'wah_{year}.parquet')

            if not os.path.exists(savepath) | utils.first_file_is_newer(wahpath, savepath):
                try:
                    df_years = []
                    files_year = [file for file in archive.namelist() if condition(file, year)]
                    if len(files_year) > 0:
                        for file in files_year:
                            df_grouped, crs = process_wah_file(file, archive)
                            df_years.append(df_grouped)
                        df_year = pd.concat(df_years)
                        gdf_year = gpd.GeoDataFrame(df_year, geometry=gpd.points_from_xy(df_year.lon, df_year.lat)).set_crs(crs).to_crs(4326)
                        if len(gdf_year) > 0:
                            gdf_year['lon'] = gdf_year['geometry'].x
                            gdf_year['lat'] = gdf_year['geometry'].y
                            gdf_year.to_parquet(savepath)
                        else:
                            warnings.warn(f"No rows created for year {year}.")
                    else:
                        pass
                except Exception as e:
                    print(e)


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
