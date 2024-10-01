#!/bin/bash
temperatures=(1.4)
number_of_trials=50
autocorrelation_window_length=100000

for temperature in ${temperatures[@]}; do
    for trial_number in $(seq 15 $number_of_trials); do # Change back later
        echo "Temperature: $temperature, Window Length: $autocorrelation_window_length, Trial Number: $trial_number"
        sbatch single_fast_autocorrelation_experiment_sbatch_job_SL2_slice_L_9.sh "$temperature" "$autocorrelation_window_length" "$trial_number"
    done
done
