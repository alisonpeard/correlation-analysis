import os
import os
os.environ['USE_PYGEOS'] = '0'
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import utils
plot_kwargs = {'bbox_inches': "tight", 'dpi': 200}

config = utils.load_config()
resultsdir = config['paths']['resultsdir']


# do assertions on a single buffer
if __name__ == "__main__":
    df = pd.read_parquet(os.path.join(resultsdir, "test", "full_timeseries", "wrz_1", "Test1.parquet"))
    df = df[df['buffer'] == 25.].reset_index()
    assert df.query("Year == 2020 and Month == 4")['ep_total'].values == 0.
    assert df.query("Year == 2020 and Month == 4")['ep_mean'].values == 0.
    assert df.query("Year == 2020 and Month == 4")['q50_anomaly_total'].values == 0. -1.
    assert df.query("Year == 2020 and Month == 4")['q75_anomaly_total'].values == 0. - np.quantile([1, 1, 0], q=.75)
    assert df.query("Year == 2020 and Month == 4")['q90_anomaly_total'].values == 0. - np.quantile([1, 1, 0], q=.9)
    assert df.query("Year == 2020 and Month == 4")['q50_deficit_total'].values == 0. + np.quantile([1, 1, 0], q=.9)
    assert df.query("Year == 2020 and Month == 8")['ep_total'].values == 2.
    assert df.query("Year == 2020 and Month == 8")['ep_mean'].values == 2.
    assert df.query("Year == 2020 and Month == 8")['q50_anomaly_total'].values == 2. - np.quantile([1, 1, 2], q=.5)
    assert df.query("Year == 2020 and Month == 8")['q75_anomaly_total'].values == 2. - np.quantile([1, 1, 2], q=.75)
    assert df.query("Year == 2020 and Month == 8")['q90_anomaly_total'].values == 2. - np.quantile([1, 1, 2], q=.9)
    assert df.query("Year == 2020 and Month == 8")['q50_deficit_total'].values == 0.

    # test event set
    df = pd.read_csv("/Users/alison/Documents/RAPID/correlation_analysis/data_results/test/events/wrz_1.csv")
    df = df[df['buffer'] == 25.].copy()
    df = df[df['backcast'] == 3].copy()
    assert len(df) == 1
    assert df['severity'].values == 3.
    assert int(df['duration'].values[0].split(' ')[0]) == 61  # check this shouldn't be 92 later
    assert df['ep_total'].values == 5.
    assert df['ep_mean'].values == 5.  # check this, currently summing monthly averages
    assert df['q50_anomaly_total'].values == 0. - np.quantile([1, 1, 0], q=.5)
    assert df['q75_anomaly_total'].values == 0. - np.quantile([1, 1, 0], q=.75)
    assert df['q90_anomaly_total'].values == 0. - np.quantile([1, 1, 0], q=.9)
    assert df['q50_deficit_total'].values == 0. + np.quantile([1, 1, 0], q=.9)


    print("All tests passed.")