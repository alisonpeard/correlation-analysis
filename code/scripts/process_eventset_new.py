"""
Use the LoS data to identify LoS events then use the timeseries to create a dataset of events for each WRZ across all ensembles.

Requires:
---------
* monthly_los_melted
* BACKCASTS
* full_timeseries

Output:
-------
"""
import os
from dateutil.relativedelta import relativedelta

import numpy as np
import pandas as pd
from tqdm import tqdm

import utils


# TODO: Check - should buffer be used as a factor while getting events?
# - if sort by date then different buffers are placed next to each other
# - so then at diff.cumsum() each buffer gets the same event number - seems right
# - so probably ok then...


def get_events(df):
    """Process melted LoS dataframe to extract events"""
    df = df.drop_duplicates()
    events_ind = df[df['LoS'] > 0].index
    df = df.loc[events_ind]
    df = df.sort_values(by=['RZ_ID', 'Year', 'Month'], ascending=True)
    
    df['Day'] = [16] * len(df)  # check w@h input data and what day of month is being used, currently it is 16th
    df['date'] = pd.to_datetime(df[['Year', 'Month', 'Day']])
    df['diff'] = df['date'].dt.to_period('M').astype(np.int64).diff().replace(np.nan, 0).astype(np.int64)
    df['event'] = (df['diff'] != 1).cumsum()

    agg_kwargs = {'start': ('date', min), 'end': ('date', max), 'severity': ('LoS', sum)}
    df = df[['RZ_ID', 'event', 'date', 'LoS']].groupby(['RZ_ID', 'event']).agg(**agg_kwargs).reset_index()
    df['duration'] = df['end'] - df['start']  # + relativedelta(months=1) (ask Jim/Anna)
    return df.reset_index()


def get_predictors_in_window(df, start, end, backcast):
    # gather all entries within time window
    backcast = start - relativedelta(months=backcast)
    end = end + relativedelta(months=1)
    date_list = pd.date_range(backcast, end, freq='m')
    year_list = [date.year for date in date_list]
    month_list = [date.month for date in date_list]
    date_df = pd.DataFrame({'Year': year_list, 'Month': month_list})
    new_df = df.merge(date_df, on=['Year', 'Month'], how='inner')
    return new_df


def sum_(x):
    if np.any(~np.isfinite(x)):
        return np.nan
    else:
        return np.sum(x)


def main(config, scenarios=['BS', 'NF', 'FF'], variables=['ep', 'prbc']):
    tempdir = config['paths']["tempdir"]
    outdir = config['paths']["resultsdir"]
    backcasts = config['config']['backcasts']

    for scenario in scenarios:

        if not os.path.exists(os.path.join(outdir, scenario.lower(), "events")):
            os.makedirs(os.path.join(outdir, scenario.lower(), "events"))

        ensembles = utils.load_list(os.path.join(tempdir, scenario.lower(), "ensembles"))

        df1 = pd.read_parquet(os.path.join(tempdir, scenario.lower(), 'indicator_series.parquet'))
        df2 = pd.read_parquet(os.path.join(tempdir, scenario.lower(), 'standardised_series.parquet'))
        wah_df = df1.merge(df2)

        wah_dfs = {}
        for variable in variables:
            wah_dfs[variable] = wah_df.loc[wah_df['Variable'] == variable]

        monthly_los_melted = pd.read_csv(os.path.join(tempdir, scenario.lower(), "monthly_los_melted.csv"))
        # for RZ_ID in [117]:  # TESTING
        for RZ_ID in (pbar:=tqdm(monthly_los_melted['RZ_ID'].unique(), total=monthly_los_melted['RZ_ID'].nunique())):
            pbar.set_description(f"Creating event set for scenario {scenario} and WRZ {RZ_ID}")

            if not os.path.exists(os.path.join(outdir, scenario.lower(), "events")):
                os.makedirs(os.path.join(outdir, scenario.lower(), "events"))

            all_events = []
            for variable in variables:

                # Subset on variable/WRZ here and on ensemble member in loop below - faster to use successive subsets
                wah_df_rz = wah_dfs[variable].loc[wah_dfs[variable]['RZ_ID'] == RZ_ID]

                # all_events = []
                nevents = 0
                # for ensemble in (sub_pbar:=tqdm(ensembles, leave=False)):
                for ensemble in ensembles:
                    # sub_pbar.set_description(f"Processing ensemble {ensemble}")

                    monthly_los = monthly_los_melted[
                        (monthly_los_melted['RZ_ID'] == RZ_ID) & (monthly_los_melted['ensemble'] == ensemble)
                    ].copy()
                    event_df = get_events(monthly_los)
                    event_df['event'] += nevents  # so not giving different events same numbers
                    nevents += len(event_df)

                    wah_buffer = wah_df_rz.loc[wah_df_rz['ensemble'] == ensemble]

                    ensemble_events = []
                    for event in event_df.itertuples():
                        for bk in backcasts:
                            wah_sub = get_predictors_in_window(wah_buffer, event.start, event.end, bk)
                            n = len(wah_sub)
                            wah_sub['backcast'] = [int(bk)] * n
                            wah_sub['start'] = [event.start] * n
                            wah_sub['end'] = [event.end] * n
                            wah_sub['event'] = [event.event] * n
                            wah_sub['severity'] = [event.severity] * n
                            wah_sub['duration'] = [event.duration] * n
                            wah_sub['ensemble'] = [ensemble] * n
                            ensemble_events.append(wah_sub)

                    if len(ensemble_events) > 0:
                        ensemble_events = pd.concat(ensemble_events)[[
                            'Variable', 'RZ_ID', 'ensemble', 'event', 'start', 'end', 'severity', 'duration',
                            'backcast', 'buffer',
                            'Value', 'anomaly_mean',  'anomaly_q50',
                            'deficit_mean', 'deficit_q50', 'deficit_q75', 'deficit_q90',
                            'si6', 'si12', 'si24',
                        ]]
                        ensemble_events = ensemble_events.groupby([
                            'Variable', 'RZ_ID', 'ensemble', 'event', 'start', 'end', 'severity', 'duration',
                            'backcast', 'buffer',
                        ]).agg(sum_).reset_index()  # ]).sum().reset_index()
                        all_events.append(ensemble_events)

            if len(all_events) > 0:
                all_events = pd.concat(all_events)
                all_events.to_csv(
                    os.path.join(outdir, scenario.lower(), "events", f"wrz_{RZ_ID}.csv"), index=False, na_rep='NA',
                )


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
