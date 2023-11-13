using Distributed
addprocs(10)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/relaxed_anneal_experiment.jl")

L_values = [10,11,12,13,14,15,16,17,18,19]

@sync @distributed for L_value in L_values
    relaxed_anneal_experiment("relaxed_anneal_L=$L_value", L_value, [1.0,0.0], 10.0,10.0,0.1,100; inherent_disorder=true, relaxation_iterations=100000)
end