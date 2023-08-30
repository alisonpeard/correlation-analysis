#!/Users/alison/mambaforge/envs/snakemake/bin/python
"""
Requires:
---------
./config.json
tempdir/wrz.gpkg
datadir/weighted_monthly_LoS_{era}.csv

Outputs:
--------
tempdir/years
tempdir/<scenario>/monthly_los_melted.csv


Note: Using Phase3_v3 for all except Phase3v1 for hist.
TODO: Modify for future scenarios

"""
import os
os.environ['USE_PYGEOS'] = '0'
import utils
import pandas as pd

# set scenario using command likne
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scenario', type=str, help='what scenario to process', default='Hist', choices=['Hist', 'BS', 'NF', 'FF', 'Test'])
args = parser.parse_args()
SCENARIO= args.scenario
print(f'\nProcessing WREW data for {SCENARIO.lower()}')


def main(config):
    datadir = config['paths']['datadir']
    outdir = config['paths']['tempdir']

    out_subdir = os.path.join(outdir, SCENARIO.lower())
    if not os.path.exists(out_subdir):
        os.makedirs(out_subdir)

    inpath = os.path.join(datadir, 'los', f'weighted_monthly_LoS_{SCENARIO}.csv')
    outpath = os.path.join(out_subdir, 'monthly_los_melted.csv')

    wrz_code = pd.read_excel(os.path.join(datadir, 'WRZ', 'wrz_code.xlsx'))
    wrz_code['RZ ID'] = pd.to_numeric(wrz_code['RZ ID'])
    wrz_dict = {wrz_code: f'{rz_id:.0f}' for wrz_code, rz_id in zip(wrz_code['WREW Code'], wrz_code['RZ ID'])}

    # load LoS data
    monthly_los = pd.read_csv(inpath, low_memory=False)
    monthly_los.columns = monthly_los.iloc[0]
    cols = ['Ensemble','Year', 'Month'] + list(monthly_los.columns)[3: len(list(monthly_los.columns))]
    monthly_los.columns = cols
    monthly_los = monthly_los.iloc[3:].reset_index(drop=True)
    monthly_los['Year'] = pd.to_numeric(monthly_los['Year'], downcast='integer')
    monthly_los['Month'] = pd.to_numeric(monthly_los['Month'], downcast='integer')

    cols = [*monthly_los.columns]
    monthly_los = monthly_los.groupby(cols).sum().reset_index()

    # melt dataframe for merging with pet data
    rz_id_cols = monthly_los.columns[3:]
    monthly_los_melted = monthly_los.melt(id_vars=['Year', 'Month', 'Ensemble'], value_vars=rz_id_cols, var_name='RZ_ID', value_name='LoS', ignore_index=True)
    monthly_los_melted = monthly_los_melted[['RZ_ID', 'Year', 'Month', 'Ensemble', 'LoS']]
    monthly_los_melted['RZ_ID'] = monthly_los_melted['RZ_ID'].map(wrz_dict)
    monthly_los_melted['LoS'] = monthly_los_melted['LoS'].astype(float)
    monthly_los_melted['Ensemble'] = monthly_los_melted['Ensemble'].apply(lambda x: f'{SCENARIO}{int(x) + 1}')
    monthly_los_melted = monthly_los_melted[monthly_los_melted['RZ_ID'] != 'nan']
    monthly_los_melted = monthly_los_melted.drop_duplicates()

    # grab how many years and ensembles we have WREW data for
    YEARS = [str(year) for year in monthly_los_melted['Year'].unique()]
    ENSEMBLES = [*monthly_los_melted['Ensemble'].unique()]
    utils.save_list(YEARS, os.path.join(out_subdir, "years"))
    utils.save_list(ENSEMBLES, os.path.join(out_subdir, 'ensembles'))
    
    print(f"Number of LoS days: {monthly_los_melted['LoS'].sum():,.0f}")
    monthly_los_melted = monthly_los_melted.rename(columns={'Ensemble': 'ensemble'})
    monthly_los_melted.to_csv(outpath)
    print(f"Saved LoS data to {outpath}")


if __name__ == "__main__":
    config = utils.load_config()
    main(config)
