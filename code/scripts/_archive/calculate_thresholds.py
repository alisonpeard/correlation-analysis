"""
Calculate thresholds (per grid cell) for anomaly/deficit indicators by pooling the specified WAH scenarios.

"""
import os
import zipfile

import numpy as np
import pandas as pd
from tqdm import tqdm

import utils
from process_indicator_yearly import condition_scenario


YEARS = {'BS': range(1975, 2004+1), 'NF': range(2020, 2049+1), 'FF': range(2070, 2099+1)}


def quantile(q):
    def _quantile(x):
        return np.quantile(x, q)
    _quantile.__name__ = "q{:d}".format(int(q * 100))  # assumes only whole number percentiles are used
    return _quantile


def read_wah_file(filename, wah_archive):
    """Read WAH scenario data stored in CSV format."""
    file = wah_archive.open(filename)
    df = pd.read_csv(file, index_col=0)
    df['lat'] = df['lat'].round(2)
    df['lon'] = df['lon'].round(2)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    df['Month'] = df['time'].dt.month
    df['prbc'] *= 86400  # mm/s to mm/d
    df['ep'] = df['prbc'] - df['pepm']  # correction to data in input files

    # Try to keep size down to read all scenarios/years
    df = df.drop(columns=['time', 'pepm'])

    df = pd.melt(
        df, id_vars=['lat', 'lon', 'Year', 'Month', 'ensemble'], value_vars=['prbc', 'ep'], var_name='Variable',
        value_name='Value',
    )

    return df


def main(config, scenarios=["BS", "NF", "FF"]):
    datadir = config['paths']["datadir"]
    outdir = config['paths']["tempdir"]

    # List will contain one dataframe per year for each WAH scenario
    dfs = []

    for scenario in (pbar := tqdm(scenarios)):
        pbar.set_description(f"Reading {scenario}")

        wahpath = os.path.join(datadir, 'w@h', f"{scenario.lower()}.zip")

        with zipfile.ZipFile(wahpath, 'r') as archive:
            for year in YEARS[scenario]:
                files_year = [file for file in archive.namelist() if condition_scenario(file, str(year))]
                if len(files_year) == 1:
                    df = read_wah_file(files_year[0], archive)
                    dfs.append(df)

    df = pd.concat(dfs)

    df_thresholds = df.groupby(['lat', 'lon', 'Month', 'Variable'])['Value'].agg([
        'mean', quantile(0.5), quantile(0.25), quantile(0.1),
    ])
    df_thresholds = df_thresholds.reset_index()

    # Change naming convention of quantiles to be "hydrological"!
    df_thresholds = df_thresholds.rename(columns={'q25': 'q75', 'q10': 'q90'})

    # Save thresholds
    savedir = os.path.join(outdir, 'thresholds')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    savepath = os.path.join(savedir, 'wah_thresholds.csv')
    df_thresholds.to_csv(savepath, index=False)


if __name__ == "__main__":
    config = utils.load_config()
    # main(config, scenarios=['NF'])  # TESTING
    main(config, scenarios=['BS', 'NF', 'FF'])
