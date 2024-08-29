#!/bin/bash

total_jobs=5
start_time=$(date +%s)

for i in {1..5}
do
    job_start_time=$(date +%s)
    echo "Starting job $i of $total_jobs at $(date)"
    julia job_energy_connectivity_cumulative_disorder_average.jl 3 20 "s9_$i"
    job_end_time=$(date +%s)
    echo "Finished job $i of $total_jobs at $(date)"
    job_duration=$((job_end_time - job_start_time))
    echo "Job $i took $job_duration seconds"
    echo "--------------------------"
done

end_time=$(date +%s)
total_duration=$((end_time - start_time))
echo "All jobs completed at $(date)"
echo "Total run time: $total_duration seconds"