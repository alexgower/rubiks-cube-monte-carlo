using Distributed
addprocs(3)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/relaxed_anneal_experiment.jl")

L_values = [7,9,11]

@sync @distributed for L_value in L_values
    relaxed_anneal_experiment("emergent_relaxed_anneal_L=$L_value", L_value, [1.0,0.0], 10.0,10.0,0.1,200; inherent_disorder=false, relaxation_iterations=10000)
end