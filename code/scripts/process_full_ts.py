"""
Create full time series for each water resource zone across all ensembles.

Requires:
---------
tempdir/wrz.gpkg
tempdir/wrz_buffer.gpkg
tempdir/<SCENARIO.lower()>/monthly_los_melted.csv
tempdir/<SCENARIO.lower()>/hist/by_basin/*.parquet

Output:
--------
resultsdir/<SCENARIO.lower()>/buffers/f'wah_{wrz_row.RZ_ID}.parquet'
resultsdir/<SCENARIO.lower()>/full_timeseries/f'wah_{wrz_row.RZ_ID}.parquet'
"""
import os
os.environ['USE_PYGEOS'] = '0'
import utils
import time
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

# set SCENARIO using command likne
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scenario', type=str, help='what scenario to process', default='Hist', choices=['Hist', 'BS', 'NF', 'FF', 'Test'])
parser.add_argument('-v', '--var', type=str, help='what indicator variable to process', default="rainfall", choices=["rainfall", 'ep'])
args = parser.parse_args()
SCENARIO = args.scenario
INDICATOR = args.var
print(f'\nProcessing data into full time series {SCENARIO.lower()} and {INDICATOR}')

agg_kwargs= {
    f"{INDICATOR}_total": (INDICATOR, sum),
    f"{INDICATOR}_mean": (INDICATOR, np.mean),
    f"{INDICATOR}_q50_total": (f"{INDICATOR}_q50", sum),
    f"{INDICATOR}_q75_total": (f"{INDICATOR}_q75", sum),
    f"{INDICATOR}_q90_total": (f"{INDICATOR}_q90", sum),
}

def get_full_ts(wah, monthly_los):
    wah["RZ_ID"] = wah['RZ_ID'].astype(int)
    monthly_los["RZ_ID"] = monthly_los['RZ_ID'].astype(int)
    full_ts = wah.merge(monthly_los, on=['RZ_ID', 'Year', 'Month'], how='inner')
    full_ts = full_ts.sort_values(by=['RZ_ID', 'Year', 'Month', 'LoS', 'buffer'])
    full_ts = full_ts[['RZ_ID', 'Year', 'Month', 'LoS', 'buffer', f'{INDICATOR}_total', f'{INDICATOR}_mean', 'q50_anomaly_total', 'q75_anomaly_total', 'q90_anomaly_total', 'q50_deficit_total']]
    return full_ts


def main(config):
    datadir = config['paths']["datadir"]
    tempdir = config['paths']["tempdir"]
    outdir = config['paths']["resultsdir"]

    wrz = gpd.read_file(os.path.join(tempdir, "wrz.gpkg"))
    wrz_buffer = gpd.read_file(os.path.join(tempdir, "wrz_buffer.gpkg"))
    monthly_los_melted = pd.read_csv(os.path.join(tempdir, SCENARIO.lower(), "monthly_los_melted.csv"))
    ENSEMBLES = utils.load_list(os.path.join(tempdir, SCENARIO.lower(), "ensembles"))
    YEARS = utils.load_list(os.path.join(tempdir, SCENARIO.lower(), "years"))

    times = []
    for wrz_row in (pbar:=tqdm(wrz.itertuples(), total=len(wrz))):
        pbar.set_description(f"Processing WRZ {wrz_row.RZ_ID}")
        for ensemble in (sub_pbar:=tqdm(ENSEMBLES, leave=False)):
            sub_pbar.set_description(f"Processing ensemble'{ensemble}'")
            pathdir_ts = os.path.join(outdir, SCENARIO.lower(), 'full_timeseries', f"wrz_{wrz_row.RZ_ID}")
            pathname_ts = os.path.join(pathdir_ts, f'{ensemble}.parquet')

            if not os.path.exists(pathname_ts):
                start = time.time()
                buffers_wrz = wrz_buffer[wrz_buffer['RZ_ID']==wrz_row.RZ_ID]

                basinpath = os.path.join(tempdir, SCENARIO.lower(), "by_basin", f"wrz_{wrz_row.RZ_ID}", f"{ensemble}.parquet")
                if os.path.exists(basinpath):
                    if not os.path.exists(pathdir_ts):
                        os.makedirs(pathdir_ts)
                    # only do overlay for one month to save time
                    wah_coords = gpd.read_parquet(basinpath, filters=[("Year", "=", int(YEARS[0])), ("Month", "=", 1)])
                    assert wah_coords.crs == buffers_wrz.crs, "Dataframes have different coordinate reference systems"
                    wah_map = gpd.overlay(wah_coords, buffers_wrz, how='intersection')
                    wah_map = wah_map[['geometry', 'RZ_ID', 'buffer']].set_index('geometry')
                    wah_gdf = gpd.read_parquet(os.path.join(tempdir, SCENARIO.lower(), "by_basin", f"wrz_{wrz_row.RZ_ID}", f"{ensemble}.parquet"))
                    assert len(wah_gdf) > 0, "Basins data has zero rows."
                    wah_buffer = wah_gdf.join(wah_map, on='geometry', how='outer', lsuffix='_', rsuffix='_r')
                    wah_buffer = wah_buffer[['RZ_ID', 'Year', 'Month', 'buffer', INDICATOR, f'{INDICATOR}_q50', f'{INDICATOR}_q75', f'{INDICATOR}_q90']].groupby(['RZ_ID', 'buffer', 'Year', 'Month']).agg(**agg_kwargs).reset_index()

                    # calculate metrics
                    wah_buffer['q50_anomaly_total'] = wah_buffer[f'{INDICATOR}_total'] - wah_buffer[f'{INDICATOR}_q50_total']
                    wah_buffer['q75_anomaly_total'] = wah_buffer[f'{INDICATOR}_total'] - wah_buffer[f'{INDICATOR}_q75_total']
                    wah_buffer['q90_anomaly_total'] = wah_buffer[f'{INDICATOR}_total'] - wah_buffer[f'{INDICATOR}_q90_total']
                    wah_buffer['q50_deficit_total'] = (wah_buffer[f'{INDICATOR}_q50_total'] - wah_buffer[f'{INDICATOR}_total']).apply(lambda x: np.max([x, 0.]))

                    # save
                    wah_buffer = wah_buffer[['RZ_ID', 'Year', 'Month', 'buffer', f'{INDICATOR}_total', f'{INDICATOR}_mean', 'q50_anomaly_total', 'q75_anomaly_total', 'q90_anomaly_total', 'q50_deficit_total']]
                    full_ts = get_full_ts(wah_buffer, monthly_los_melted[monthly_los_melted['ensemble'] == ensemble].copy())
                    full_ts.to_parquet(pathname_ts)

                    # time update
                    duration = (time.time() - start) / 60
                    times.append(duration)
                else:
                    warnings.warn(f"No basin file for wrz {wrz_row.RZ_ID}")

    if len(times) > 0:
        print(f"{np.mean(times):.2f} minutes per wrz")
    else:
        print(f"There was nothing to process.")


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
