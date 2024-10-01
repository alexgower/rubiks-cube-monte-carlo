#!/bin/bash

number_of_trials=100

for trial_number in $(seq 1 $number_of_trials); do
    echo "Trial Number: $trial_number"
    sbatch single_energy_connectivity_disorder_average.sh 11 $trial_number "HPC"
done
