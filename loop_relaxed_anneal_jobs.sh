#!/bin/bash


models_to_test=("inherent_disorder")
L_values=(11)
swap_move_probabilities=(1.0)
number_of_trials=50
# trial_numbers=(17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44)
# trial_numbers=(1 2 3 4 5 6 7 8 9 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50)
trial_numbers=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50)

for model in "${models_to_test[@]}"; do
    for L in "${L_values[@]}"; do
        for swap_move_probability in "${swap_move_probabilities[@]}"; do
            # for trial_number in $(seq 1 "$number_of_trials"); do
            for trial_number in "${trial_numbers[@]}"; do
                echo "Submitting job for model: $model, L: $L, swap_move_probability: $swap_move_probability, trial_number: $trial_number"
                sbatch single_relaxed_anneal_sbatch_job_SL2.sh "$model" "$L" "$swap_move_probability" "$trial_number"
            done
        done
    done
done