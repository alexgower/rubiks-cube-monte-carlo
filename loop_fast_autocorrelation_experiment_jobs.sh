#!/bin/bash


temperatures=(0.9 0.91 0.92 0.95)
number_of_trials=50
autocorrelation_window_length=250000


for temperature in ${temperatures[@]}; do
    for trial_number in $(seq 1 $number_of_trials); do
        echo "Temperature: $temperature, Window Length: $autocorrelation_window_length, Trial Number: $trial_number"
        sbatch single_fast_autocorrelation_experiment_sbatch_job_SL2.sh "$temperature" "$autocorrelation_window_length" "$trial_number"
    done
done
