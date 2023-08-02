"""
Load the basins and create the union of all basins intersecting each WRZ. Gather the W@H data for each union of basins and save on a per-wrz basis.

Requires:
---------
./config.json
tempdir/wrz.gpkg
tempdir/years
datdir/river_basins/ukcp18-uk-land-river-hires.gpkg
tempdir/<scenario>/yearly/wah_*.parquet

Output:
-------
tempdir/<scenario>/by_basin/wah_*.parquet

"""
import os
os.environ['USE_PYGEOS'] = '0'
import utils
import warnings
import pandas as pd
import geopandas as gpd
import dask.dataframe as dd
from tqdm import tqdm

# set scenario using command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scenario', type=str, help='what scenario to process', default='BS', choices=['Hist', 'BS', 'NF', 'FF', 'Test'])
parser.add_argument('-v', '--var', type=str, help='what indicator variable to process', default='ep', choices=['rainfall', 'ep'])
args = parser.parse_args()
SCENARIO= args.scenario
INDICATOR = args.var

print(f'\nProcessing indicator data into basins for {SCENARIO.lower()} and {INDICATOR}')

def main(config):
    datadir = config['paths']["datadir"]
    outdir = config['paths']["tempdir"]
    YEARS = utils.load_list(os.path.join(outdir, SCENARIO.lower(), "years"))
    ENSEMBLES = utils.load_list(os.path.join(outdir, SCENARIO.lower(), "ensembles"))

    wrz = gpd.read_file(os.path.join(outdir, 'wrz.gpkg'))
    basins = gpd.read_file(os.path.join(datadir, 'river_basins', 'ukcp18-uk-land-river-hires.gpkg')).to_crs(4326)
    assert basins.crs == wrz.crs
    wrz_basins = wrz[['RZ_ID', 'geometry']].sjoin(basins[['id', 'geometry']], how='left', predicate='intersects') # 36 -> 15
    basin_geoms = {idx: geom for idx, geom in zip(basins.id, basins.geometry)}
    wrz_basins['geometry'] = wrz_basins['id'].apply(lambda x: basin_geoms[x])
    wrz_basins = wrz_basins.set_geometry('geometry')
    basins_dissolved = wrz_basins.dissolve(by='RZ_ID').reset_index()

    rz_basin_geoms = {rz_id: basin_geom for rz_id, basin_geom in zip(basins_dissolved['RZ_ID'], basins_dissolved['geometry'])}
    assert basins_dissolved['RZ_ID'].nunique() == wrz['RZ_ID'].nunique(), 'Number of RZ IDs changed.'
    wrz = wrz[['RZ_ID', 'geometry']].drop_duplicates()
    wrz['bounds'] = wrz.bounds.apply(lambda row: (row['minx'], row['miny'], row['maxx'], row['maxy']), axis=1)

    for ensemble in (pbar:=tqdm(ENSEMBLES)):
        pbar.set_description(f"Processing ensemble {ensemble} into basins")
        for wrz_row in (sub_pbar:=tqdm(wrz.itertuples(), total=len(wrz), leave=False)):
            sub_pbar.set_description(f"WRZ {wrz_row.RZ_ID}")
            savedir = os.path.join(outdir, SCENARIO.lower(), "by_basin", f"wrz_{wrz_row.RZ_ID}")
            savepath = os.path.join(savedir, f"{ensemble}.parquet")

            if not os.path.exists(savepath) | utils.first_file_is_newer(os.path.join(outdir, SCENARIO.lower(), 'yearly', f'wah_{YEARS[0]}.parquet'), savepath):
                bounds = rz_basin_geoms[wrz_row.RZ_ID].bounds
                wah_dfs = []
                for year in YEARS:
                    path_yearly = os.path.join(outdir, SCENARIO.lower(), 'yearly', f'wah_{year}.parquet')
                    if os.path.exists(path_yearly):
                        ddf = dd.read_parquet(path_yearly,
                                            filters=[[('lon', '>=', bounds[0]),
                                                ('lat', '>=', bounds[1]),
                                                ('lon', '<=', bounds[2]),
                                                ('lat', '<=', bounds[3]),
                                                ('ensemble', '==', ensemble)]],
                                            columns=['lon', 'lat', 'Year', 'Month', INDICATOR])
                        df = ddf.compute()
                        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326).clip(rz_basin_geoms[wrz_row.RZ_ID])
                        if len(gdf) > 0:
                            wah_dfs.append(gdf)
                
                if len(wah_dfs) > 0:
                    wah_df = pd.concat(wah_dfs)
                    # calculate all the metrics
                    q50 = wah_df.groupby("Month")[INDICATOR].quantile(.5).to_frame()
                    q75 = wah_df.groupby("Month")[INDICATOR].quantile(.75).to_frame()
                    q90 = wah_df.groupby("Month")[INDICATOR].quantile(.9).to_frame()
                    wah_df = pd.merge(wah_df, q50, on='Month', how='left', suffixes=('', '_q50'))
                    wah_df = pd.merge(wah_df, q75, on='Month', how='left', suffixes=('', '_q75'))
                    wah_df = pd.merge(wah_df, q90, on='Month', how='left', suffixes=('', '_q90'))
                    if len(wah_df) > 0:
                        if not os.path.exists(savedir):
                            os.makedirs(savedir)
                        sub_pbar.set_description(f"Saving WRZ {wrz_row.RZ_ID} results to {savepath}")
                        wah_df[['Year', 'Month', INDICATOR, f'{INDICATOR}_q50',f'{INDICATOR}_q75', f'{INDICATOR}_q90', 'geometry']].to_parquet(savepath) # check subset ok
                else:
                    warnings.warn(f"No files found for wrz {wrz.RZ_ID}")


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
