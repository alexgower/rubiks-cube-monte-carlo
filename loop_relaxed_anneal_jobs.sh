#!/bin/bash


models_to_test=("clean" "inherent_disorder")
L_values=(3 5 7 9 11)
swap_move_probabilities=(0.0 1.0)
number_of_trials=50

# models_to_test=("inherent_disorder")
# L_values=(3)
# swap_move_probabilities=(0.0)
# number_of_trials=1


for model in "${models_to_test[@]}"; do
    for L in "${L_values[@]}"; do
        for swap_move_probability in ${!swap_move_probabilities[@]}; do
            for trial_number in $(seq 1 $number_of_trials); do
                sbatch single_relaxed_anneal_sbatch_job.sh "$model" "$L" "$swap_move_probability" "$trial_number"
            done
        done
    done
done