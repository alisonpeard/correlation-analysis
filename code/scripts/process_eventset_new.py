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
import numpy as np
import utils
from dateutil.relativedelta import relativedelta
import pandas as pd
from tqdm import tqdm


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

    agg_kwargs = {'start': ('date', min),
                  'end': ('date', max),
                  'severity': ('LoS', sum)}
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

# set scenario using command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scenario', type=str, help='what scenario to process', default='Hist', choices=['Hist', 'BS', 'NF', 'FF', 'Test'])
parser.add_argument('-v', '--var', type=str, help='what indicator variable to process', default="rainfall", choices=["rainfall", 'ep'])
args = parser.parse_args()
SCENARIO = args.scenario
INDICATOR = args.var
print(f'\nProcessing events for {SCENARIO.lower()} and {INDICATOR}')


def main(config):
    datadir = config['paths']["datadir"]
    tempdir = config['paths']["tempdir"]
    outdir = config['paths']["resultsdir"]
    backcasts = config['config']['backcasts']
    ENSEMBLES = utils.load_list(os.path.join(tempdir, SCENARIO.lower(), "ensembles"))

    monthly_los_melted = pd.read_csv(os.path.join(tempdir, SCENARIO.lower(), "monthly_los_melted.csv"))
    # for RZ_ID in (pbar := tqdm([117], total=1)):  # TESTING
    for RZ_ID in (pbar:=tqdm(monthly_los_melted['RZ_ID'].unique(), total=monthly_los_melted['RZ_ID'].nunique())):
        pbar.set_description(f"Creating event set for WRZ {RZ_ID}")
        if not os.path.exists(os.path.join(outdir, SCENARIO.lower(), "events")):
            os.makedirs(os.path.join(outdir, SCENARIO.lower(), "events"))

        all_events = []
        nevents = 0
        for ensemble in (sub_pbar:=tqdm(ENSEMBLES, leave=False)):
            sub_pbar.set_description(f"Processing ensemble {ensemble}")

            wah_buffer = pd.read_parquet(os.path.join(outdir, SCENARIO.lower(), 'full_timeseries', f"wrz_{RZ_ID}", f'{ensemble}.parquet'))
            monthly_los = monthly_los_melted[(monthly_los_melted['RZ_ID'] == RZ_ID) & (monthly_los_melted['ensemble'] == ensemble)].copy()
            event_df = get_events(monthly_los)
            event_df['event'] += nevents  # so not giving different events same numbers
            nevents += len(event_df)

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
                    'RZ_ID', 'ensemble', 'event', 'start', 'end', 'severity', 'duration', 'backcast', 'buffer',
                    f'{INDICATOR}_total', f'{INDICATOR}_mean', 'mean_anomaly_total',  'q50_anomaly_total',
                    'mean_deficit_total', 'q50_deficit_total', 'q75_deficit_total', 'q90_deficit_total',
                ]]
                ensemble_events = ensemble_events.groupby([
                    'RZ_ID', 'ensemble', 'event', 'start', 'end', 'severity', 'duration', 'backcast', 'buffer'
                ]).sum().reset_index()
                all_events.append(ensemble_events)
        
        if len(all_events) > 0:
            all_events = pd.concat(all_events)
            # assert check same number of events and max event number
            all_events.to_csv(os.path.join(outdir, SCENARIO.lower(), "events", f"wrz_{RZ_ID}.csv"), index=False)


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
