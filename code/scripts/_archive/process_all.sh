#!/bin/bash

# SLURM COMMANDS HERE

#INDEX=1  # using $1 instead
INDICATOR="ep"
SCENARIOS=("Hist" "BS" "NF" "FF" "Test")
SCENARIO=${SCENARIOS[$1]}

source /Users/alison/mambaforge/bin/activate general

python process_wrz_data.py
python process_los_data.py -s $SCENARIO
python process_indicator_yearly.py -s $SCENARIO -v $INDICATOR
python process_indicator_for_basins.py -s $SCENARIO -v $INDICATOR
python process_full_ts.py -s $SCENARIO -v $INDICATOR
python process_eventset_new.py -s $SCENARIO -v $INDICATOR
