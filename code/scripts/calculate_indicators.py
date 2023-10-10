import os

import numpy as np
import pandas as pd
from tqdm import tqdm

import utils


def calculate_thresholds(tempdir, scenario, df):
    df_thresholds = df.groupby(['RZ_ID', 'buffer', 'Variable', 'Month'])['Value'].agg([
        'mean', utils.quantile(0.5), utils.quantile(0.25), utils.quantile(0.1),
    ])
    df_thresholds = df_thresholds.reset_index()

    # Change naming convention of quantiles to be "hydrological"!
    df_thresholds = df_thresholds.rename(columns={'q25': 'q75', 'q10': 'q90'})

    thresholds_path = os.path.join(tempdir, scenario.lower(), 'thresholds.csv')
    df_thresholds.to_csv(thresholds_path, index=False)

    return df_thresholds


def calculate_indicators(tempdir, scenario, df, df_thresholds):
    # monthly anomaly and deficit time series

    df = df.merge(df_thresholds, on=['RZ_ID', 'buffer', 'Variable', 'Month'])
    df = df.sort_values(['RZ_ID', 'buffer', 'Variable', 'ensemble', 'Year', 'Month'])

    for statistic in ['mean', 'q50', 'q75', 'q90']:
        df[f'anomaly_{statistic}'] = df['Value'] - df[statistic]
        df[f'deficit_{statistic}'] = np.where(
            df[f'anomaly_{statistic}'] < 0.0, df[f'anomaly_{statistic}'], 0.0
        )
        if statistic in ['q75', 'q90']:
            df = df.drop(columns=f'anomaly_{statistic}')

    df = df.drop(columns=['mean', 'q50', 'q75', 'q90'])

    output_path = os.path.join(tempdir, scenario.lower(), 'indicator_series.parquet')
    df.to_parquet(output_path, index=False)


def main(config, scenarios=["BS", "NF", "FF"]):
    tempdir = config['paths']["tempdir"]

    for scenario in (pbar := tqdm(scenarios)):
        pbar.set_description(f"Calculating indicators {scenario}")

        input_path = os.path.join(tempdir, scenario.lower(), 'aggregated_series.parquet')
        df_agg = pd.read_parquet(input_path)
        df_thresholds = calculate_thresholds(tempdir, scenario, df_agg)
        calculate_indicators(tempdir, scenario, df_agg, df_thresholds)


if __name__ == '__main__':
    config = utils.load_config()
    main(config)
