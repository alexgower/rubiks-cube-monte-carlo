include("time_experiment.jl")

function thin_funnel_experiment(N_t::Int64; normalization::String="solved")
    E = 4
    T = 5.0
    p_swap = 0.0

    # TODO redo E=2 for completion?

    # E_3_starting_configurations = [ [[3 4 4 4 4 6 6 5 2; 6 4 4 4 4 4 6 4 5; 4 2 4 4 2 2 5 5 5; 6 5 4 4 1 3 5 4 3; 4 5 4 1 1 1 5 4 4; 4 5 5 4 4 1 3 4 4; 4 5 3 5 5 5 3 3 2; 2 5 3 1 1 1 3 2 2; 5 5 2 1 1 1 2 2 4], [4 4 6 6 5 4 6 6 1; 4 2 4 6 4 6 6 6 4; 4 2 2 6 4 6 6 6 3; 4 4 4 4 4 2 2 3 3; 2 1 1 4 2 2 3 3 3; 2 6 6 5 2 2 2 1 1; 6 6 6 6 4 1 1 4 1; 6 3 1 5 5 1 1 1 1; 6 3 1 3 3 1 1 1 1], [2 5 4 5 5 5 5 2 3; 4 1 1 2 6 6 4 3 3; 4 1 1 1 6 2 4 3 1; 1 1 1 1 1 2 6 3 1; 1 1 1 5 3 2 2 3 3; 5 5 1 1 3 3 3 3 3; 5 3 1 1 3 3 3 3 3; 5 5 5 6 6 5 5 3 3; 5 5 1 2 1 5 2 3 2], [1 3 3 6 4 5 3 5 6; 1 6 6 6 6 5 5 5 5; 5 6 6 6 6 5 5 5 5; 3 3 3 6 6 5 4 4 5; 6 3 3 3 4 4 4 4 6; 3 3 3 3 3 4 4 4 6; 3 3 3 3 3 3 6 1 6; 3 3 3 3 3 3 6 6 6; 3 3 3 3 3 3 1 6 6], [3 1 1 1 5 5 4 1 4; 2 2 1 1 5 5 4 1 1; 2 2 2 2 5 1 4 4 3; 2 2 4 6 3 1 1 1 2; 2 2 2 5 5 5 1 1 5; 2 2 2 5 5 3 6 1 5; 5 2 2 1 1 3 5 5 5; 6 6 6 2 2 2 2 1 3; 5 4 6 2 2 2 2 4 4], [5 2 2 4 2 6 5 6 6; 4 4 4 4 2 3 5 5 6; 3 1 1 4 6 5 5 1 6; 4 2 5 6 6 6 5 6 6; 6 5 5 6 6 6 6 6 6; 6 5 2 2 2 5 6 6 2; 6 2 2 2 2 4 4 4 2; 1 2 2 2 2 2 2 4 2; 1 1 1 1 1 4 4 2 2]],
    #                             [[1 1 1 1 1 1 2 1 1; 6 1 1 1 6 6 1 1 1; 4 2 2 1 6 6 2 1 1; 5 5 2 4 1 1 1 1 1; 5 5 5 5 1 1 1 1 2; 5 5 5 1 1 1 1 2 4; 5 5 1 3 3 1 5 2 2; 5 3 3 3 3 3 3 2 2; 6 1 3 3 3 2 2 2 6], [5 5 5 2 4 4 4 3 4; 5 5 5 2 4 4 4 4 4; 5 5 2 2 4 4 4 4 4; 2 1 1 1 2 2 1 4 4; 3 1 3 3 2 6 6 4 3; 3 3 3 3 5 6 6 6 6; 3 3 5 5 5 4 6 6 6; 6 6 6 4 4 4 6 6 6; 6 6 6 4 4 4 6 4 6], [4 6 1 3 3 3 1 5 5; 6 1 4 4 3 3 3 5 3; 5 5 5 4 3 3 3 2 2; 5 5 5 2 3 3 3 4 5; 2 5 5 5 3 3 3 5 5; 2 6 6 3 3 3 3 1 5; 2 2 6 2 4 3 4 1 1; 2 2 6 3 3 2 4 2 3; 2 4 4 6 6 6 6 6 3], [3 5 6 5 5 3 3 3 3; 3 6 6 6 1 1 4 3 3; 1 6 6 6 1 1 1 3 4; 1 6 6 6 2 6 5 5 4; 1 1 4 4 4 4 4 4 4; 6 1 4 4 4 5 3 4 4; 1 1 4 4 1 5 4 4 2; 1 1 4 2 2 6 2 4 4; 1 1 4 2 2 6 4 2 2], [3 1 3 1 4 6 6 2 2; 1 3 3 3 6 1 1 4 2; 1 3 3 2 2 1 1 4 3; 1 4 2 2 2 2 2 2 2; 1 3 2 2 5 4 2 2 2; 3 3 2 5 1 4 4 2 2; 2 1 3 2 1 4 2 5 3; 3 2 2 2 2 5 5 5 2; 5 2 2 2 6 6 6 4 4], [4 4 4 1 6 6 5 5 5; 4 4 6 6 6 5 5 5 5; 5 2 6 6 2 5 5 2 5; 5 6 6 6 6 5 5 2 4; 6 6 6 6 6 5 5 5 1; 3 1 3 4 6 5 5 5 1; 5 6 3 4 6 6 1 1 6; 3 6 3 3 2 5 5 3 5; 2 4 3 3 5 5 3 6 1]],
    #                             [[5 5 3 3 1 3 3 3 3; 5 5 3 3 3 3 3 1 1; 4 3 3 3 3 3 3 5 5; 1 1 1 1 1 3 3 5 5; 2 1 5 5 1 3 3 2 2; 2 2 5 5 1 3 3 3 3; 2 5 5 3 3 3 3 1 3; 1 4 4 4 3 3 6 1 1; 1 4 4 2 2 1 4 1 1], [4 3 1 3 2 6 6 2 2; 4 3 3 3 2 6 4 4 2; 5 5 4 4 4 6 2 4 1; 5 5 4 4 5 3 3 2 2; 5 5 5 4 2 2 2 2 4; 5 5 5 1 4 4 1 5 1; 2 5 5 5 1 4 1 2 2; 2 2 1 2 1 6 1 2 2; 2 2 2 2 5 6 6 4 4], [2 2 1 2 3 4 1 4 2; 2 2 2 2 5 3 1 1 2; 2 2 2 2 3 3 6 1 2; 2 2 2 2 2 2 2 3 4; 4 2 2 3 3 2 6 6 4; 4 6 6 3 6 6 6 6 2; 5 6 6 2 2 5 3 3 5; 5 6 6 4 4 4 3 3 3; 5 6 6 4 6 3 3 3 3], [6 5 6 6 1 6 6 5 5; 5 5 5 5 5 5 5 5 5; 5 5 5 1 5 5 5 2 5; 5 1 5 1 1 1 4 1 5; 6 4 1 1 4 4 4 4 3; 4 4 1 4 4 4 4 4 4; 4 4 1 4 6 4 2 6 4; 4 5 5 3 6 4 4 3 6; 5 5 5 3 6 6 6 6 6], [3 3 3 1 3 3 3 3 3; 3 3 1 1 3 2 2 2 4; 3 3 1 1 2 2 2 2 1; 5 1 1 2 2 2 2 4 2; 5 5 5 5 5 3 1 4 3; 5 5 5 5 5 5 5 6 6; 5 4 4 2 4 2 4 4 1; 1 4 4 4 1 2 2 4 4; 1 1 1 1 1 3 4 3 4], [6 6 6 6 6 5 2 6 6; 6 6 6 6 6 5 2 6 6; 4 6 6 6 6 6 1 1 2; 4 6 6 6 6 5 1 1 1; 5 6 6 6 6 6 1 1 1; 4 1 6 6 3 6 1 2 6; 4 3 4 4 4 6 6 6 6; 6 1 1 1 3 6 6 6 1; 4 1 1 1 4 1 3 4 1]],
    #                             [[1 2 5 5 4 4 4 2 2; 1 2 2 5 5 3 4 2 4; 2 2 2 2 1 5 4 4 4; 3 5 2 2 1 5 4 4 1; 6 2 2 2 1 4 4 4 1; 5 2 2 5 6 6 5 3 3; 6 4 1 1 3 3 3 6 3; 6 4 1 1 5 6 6 6 3; 6 4 3 1 5 4 6 6 3], [5 1 3 3 3 1 1 5 5; 6 4 4 4 2 2 3 2 1; 6 6 6 5 2 6 2 1 1; 3 1 5 6 2 2 2 2 3; 1 1 1 5 2 2 1 2 3; 1 1 1 1 1 1 1 4 4; 1 1 1 1 1 1 1 4 4; 1 1 1 1 1 1 1 3 3; 2 1 1 1 1 1 1 3 3], [4 6 2 2 2 2 2 2 2; 2 2 2 6 3 3 5 5 2; 2 2 2 6 3 3 3 3 4; 2 2 3 3 3 3 3 3 4; 3 2 5 3 3 3 3 3 3; 3 3 5 3 3 3 3 3 3; 4 4 6 4 4 5 5 3 3; 4 4 4 4 4 5 2 3 3; 4 4 4 4 4 4 2 4 3], [4 5 4 2 4 6 3 3 4; 4 6 3 4 4 2 2 4 2; 4 6 1 6 5 4 2 2 2; 1 6 2 2 4 4 4 2 2; 1 3 5 4 4 1 2 6 2; 3 3 3 4 5 5 5 5 2; 3 3 4 4 5 5 5 3 2; 3 5 2 2 1 5 5 5 2; 3 4 6 5 5 5 5 5 6], [1 3 3 4 2 6 5 5 5; 1 3 3 3 3 1 1 1 5; 5 5 3 3 3 3 3 1 3; 5 4 2 2 1 1 1 1 2; 5 5 4 2 5 5 2 1 4; 6 5 4 4 4 4 2 5 4; 6 5 4 4 4 4 4 4 1; 6 1 5 4 4 4 6 3 1; 6 5 5 2 6 6 6 3 1], [6 2 5 6 6 6 6 6 2; 6 6 6 2 6 6 5 6 4; 6 3 6 2 6 6 6 6 2; 5 6 1 5 6 6 6 6 6; 5 5 6 5 6 6 6 6 2; 5 6 1 1 6 6 6 5 5; 5 5 5 6 6 6 5 1 1; 5 5 6 6 6 1 5 1 6; 5 5 5 6 6 1 1 1 1]],
    #                             [[5 2 2 4 3 4 4 4 4; 5 3 3 2 2 4 4 4 4; 1 3 3 3 6 4 4 4 3; 1 3 4 1 1 1 1 4 3; 1 1 1 1 1 1 1 1 4; 1 1 1 2 2 3 5 5 4; 1 1 1 1 1 5 6 6 6; 1 1 1 1 1 1 6 6 6; 1 1 1 1 1 1 6 6 6], [4 5 5 3 3 3 3 3 3; 4 5 4 3 3 3 3 3 3; 5 4 5 3 3 6 6 3 3; 5 4 2 6 1 6 3 3 3; 4 2 1 5 2 2 3 3 3; 4 6 6 3 2 2 1 3 3; 3 3 3 3 2 2 4 4 3; 3 2 2 2 2 1 1 5 3; 1 1 1 5 2 6 6 3 3], [2 2 6 3 2 2 5 5 6; 6 1 4 5 6 5 5 2 4; 4 5 5 6 6 2 2 2 4; 1 1 2 4 6 2 2 5 2; 1 2 2 2 3 3 3 3 2; 1 2 6 1 3 3 5 3 2; 1 1 1 1 5 5 5 1 2; 3 1 1 1 4 2 2 2 1; 3 2 4 4 4 2 4 2 2], [4 6 4 4 4 4 2 2 2; 4 4 4 3 4 2 3 1 1; 1 4 4 4 4 4 3 3 3; 3 4 4 4 4 1 5 4 2; 3 4 4 4 4 4 2 5 2; 4 4 4 4 4 2 2 1 2; 2 2 2 5 4 2 2 2 2; 6 2 5 5 5 4 5 4 4; 6 6 3 2 5 5 4 4 4], [6 6 2 2 5 6 6 2 2; 6 6 2 2 6 6 6 4 2; 6 6 2 2 2 1 4 6 6; 6 2 3 5 5 5 3 6 6; 6 4 5 5 5 5 6 5 5; 5 5 6 6 3 5 3 6 5; 1 5 6 6 5 5 3 2 5; 1 3 3 6 6 6 6 6 5; 3 3 3 3 6 6 5 5 5], [5 5 5 5 5 5 5 5 5; 1 5 5 5 5 5 5 5 5; 2 1 1 4 5 5 5 5 5; 6 6 6 4 6 5 1 3 6; 6 6 6 6 6 3 3 3 6; 5 6 4 6 6 3 3 2 6; 2 6 6 6 4 1 1 2 6; 4 6 6 4 1 1 1 3 3; 1 2 4 1 1 1 1 1 1]],
    #                             [[5 6 6 4 2 2 2 2 3; 6 5 5 3 2 2 4 4 4; 6 4 3 3 4 4 4 4 4; 6 6 5 5 1 2 4 5 1; 5 3 1 2 1 2 4 2 6; 5 2 2 2 2 2 6 6 6; 1 2 2 6 5 1 6 5 5; 2 2 2 3 3 4 5 6 6; 2 6 2 4 4 6 4 2 4], [5 5 2 3 3 3 3 3 6; 5 5 2 5 5 5 3 3 3; 5 2 5 5 5 3 3 3 2; 2 2 5 1 5 3 3 2 2; 2 2 2 2 2 4 2 1 3; 5 1 6 6 5 4 2 6 4; 5 6 6 2 2 2 2 2 4; 1 4 6 5 5 5 5 6 3; 1 1 1 5 5 5 5 6 4], [5 3 3 3 2 1 1 1 5; 5 3 3 6 6 1 1 1 3; 3 3 3 6 6 1 1 1 3; 3 3 3 3 1 1 1 1 1; 3 3 3 3 3 3 1 1 1; 3 3 3 3 3 6 1 1 1; 1 1 5 2 3 1 1 1 1; 1 1 1 3 3 6 6 1 1; 3 1 6 5 5 6 6 5 1], [3 1 1 1 1 1 3 3 6; 4 1 1 1 1 1 3 3 6; 2 5 1 1 1 1 4 3 2; 2 4 1 1 1 3 4 3 2; 1 4 4 4 4 4 4 5 2; 4 4 4 4 4 4 4 2 2; 3 4 4 4 1 4 4 3 4; 5 4 4 4 1 4 1 6 4; 4 4 4 6 6 4 4 4 3], [4 3 3 3 1 1 1 2 2; 4 4 4 1 4 1 1 2 2; 4 4 3 2 2 3 1 4 4; 4 4 4 4 6 5 2 4 4; 4 4 5 5 5 3 3 2 4; 5 3 5 5 5 2 2 2 4; 5 3 5 5 5 3 6 6 6; 1 3 5 3 6 6 6 6 5; 1 2 2 3 6 6 6 6 6], [6 4 6 6 6 5 3 2 2; 5 5 6 6 6 4 6 2 2; 5 5 6 6 6 3 2 2 2; 5 5 6 6 6 1 5 2 2; 4 4 6 6 6 1 3 5 5; 2 2 6 6 6 5 5 5 6; 1 2 2 6 6 5 5 6 6; 4 5 2 6 6 5 5 2 6; 1 3 5 3 3 1 5 5 2]],
    #                             [[5 6 6 4 4 3 3 3 3; 5 5 6 4 4 4 6 3 3; 1 1 3 6 4 2 2 2 3; 1 1 3 5 5 4 4 5 5; 1 1 3 1 1 2 5 5 5; 2 3 3 1 1 3 3 3 4; 2 4 4 1 1 4 6 6 4; 5 4 4 5 6 6 6 6 6; 5 4 3 1 5 5 6 4 4], [6 6 6 6 6 6 6 5 1; 6 6 6 6 6 6 6 4 4; 4 3 6 2 1 1 4 4 4; 4 1 6 1 1 4 4 4 4; 4 2 6 6 2 4 4 4 4; 2 1 2 2 2 4 4 4 4; 2 5 3 5 5 1 4 4 4; 4 5 1 5 5 4 4 4 4; 1 1 4 3 5 5 5 4 4], [6 1 1 1 2 2 5 5 5; 4 1 1 1 1 2 5 5 5; 5 1 1 1 2 2 1 3 4; 5 5 3 1 2 2 1 2 4; 3 3 3 3 3 5 1 1 1; 3 3 3 3 3 3 4 1 1; 3 3 4 4 3 3 2 3 1; 3 3 4 3 3 3 3 3 4; 3 3 3 3 3 3 1 3 4], [3 3 4 2 2 4 1 6 6; 3 4 5 6 2 4 4 6 6; 2 2 1 5 5 2 3 3 6; 2 2 1 2 4 3 3 3 3; 2 2 4 2 4 3 1 3 3; 2 2 2 2 4 5 6 6 1; 2 2 2 2 2 2 2 1 1; 2 2 2 4 4 2 1 1 1; 2 2 2 6 6 1 1 1 4], [2 2 5 6 6 5 5 1 2; 2 2 5 2 6 2 2 5 2; 2 2 5 5 6 5 1 4 2; 5 1 1 5 5 5 6 6 6; 4 1 2 5 5 4 4 6 6; 1 1 1 1 1 4 4 3 3; 3 1 5 5 2 4 3 3 3; 1 1 3 2 2 4 5 2 2; 1 1 6 6 1 6 5 2 2], [6 6 3 3 1 5 5 5 5; 6 3 2 5 5 5 5 2 5; 6 5 5 5 5 6 6 6 5; 2 3 3 6 6 6 6 5 5; 3 3 3 3 6 6 6 4 5; 1 1 5 6 6 6 6 6 6; 1 1 5 5 6 6 6 5 6; 1 6 6 6 5 5 2 1 3; 1 2 2 2 2 4 4 5 3]],
    #                             [[2 3 3 3 3 1 3 3 3; 3 3 3 1 3 1 3 3 3; 6 6 3 2 1 4 2 2 5; 4 4 1 1 4 5 4 4 4; 4 4 4 1 1 1 1 3 5; 4 4 4 1 1 1 3 3 5; 4 6 1 1 4 4 4 4 4; 5 4 4 4 4 3 3 5 5; 5 5 5 5 5 5 6 5 5], [4 4 1 5 5 3 2 4 4; 1 1 1 2 2 3 2 2 2; 1 1 6 2 2 5 3 2 6; 4 2 2 2 2 2 2 2 3; 2 2 2 2 2 2 2 2 2; 2 2 2 2 2 2 2 2 2; 2 4 2 3 3 1 1 2 4; 2 5 5 3 3 3 1 6 2; 1 4 6 4 6 4 4 4 6], [3 2 2 2 1 3 3 2 1; 3 2 3 3 1 2 2 2 1; 5 5 3 3 3 2 2 2 6; 5 5 3 3 1 4 1 1 6; 1 1 2 3 3 5 1 1 4; 1 1 5 4 3 5 1 1 2; 1 1 5 2 4 5 2 2 2; 1 4 5 2 2 2 5 5 6; 1 6 2 2 1 5 5 5 6], [3 3 3 3 3 1 1 1 1; 3 3 5 5 3 4 4 6 6; 1 3 5 6 3 1 4 4 5; 5 5 6 6 6 1 1 4 2; 5 1 6 4 4 5 5 5 4; 5 5 6 4 4 5 5 1 6; 5 5 6 6 6 5 5 6 5; 6 5 6 6 6 5 6 6 5; 4 5 5 6 6 6 2 1 6], [5 4 4 1 3 3 6 6 6; 4 4 4 5 5 6 1 1 6; 4 5 5 5 5 3 1 5 1; 2 5 5 5 3 3 3 6 6; 2 6 5 5 5 3 3 4 6; 2 6 4 4 5 3 4 6 6; 4 4 4 4 4 3 4 6 6; 5 4 4 4 4 4 2 2 2; 5 4 4 4 4 4 2 2 2], [2 6 6 6 6 6 3 1 3; 2 6 6 6 6 5 3 3 3; 3 3 6 6 6 4 3 3 3; 3 3 6 6 4 3 3 3 3; 2 6 6 6 6 6 5 5 3; 1 6 6 6 6 6 5 1 1; 3 6 6 6 1 1 1 1 2; 4 1 1 6 5 1 1 1 6; 2 1 1 1 1 1 1 1 4]],
    #                             [[2 3 3 6 6 6 4 4 4; 6 6 6 6 6 3 5 4 4; 6 6 1 6 1 1 5 2 4; 6 1 1 1 1 1 1 1 4; 1 1 1 1 1 1 2 2 4; 1 1 1 1 1 1 1 2 4; 1 1 1 4 3 1 1 2 2; 1 1 1 2 2 2 1 2 4; 3 3 5 2 6 2 1 4 4], [4 4 6 2 2 2 1 1 1; 4 4 1 4 2 6 6 1 1; 1 2 1 4 2 6 6 5 6; 1 1 2 2 2 6 6 5 6; 4 6 2 2 2 2 6 3 1; 2 6 6 3 3 2 4 3 1; 4 6 6 3 6 5 3 3 1; 4 4 5 5 3 1 3 1 1; 5 6 5 1 1 1 1 1 1], [2 2 2 3 1 3 3 5 5; 2 3 2 5 1 2 4 6 6; 2 2 2 4 1 1 4 6 6; 2 2 3 3 3 3 3 6 6; 2 1 3 3 3 2 1 6 5; 1 1 5 5 5 2 6 6 6; 2 1 5 5 5 1 6 6 6; 2 2 4 4 2 4 5 5 2; 3 5 3 3 3 5 5 5 2], [2 2 1 1 2 6 6 6 1; 2 5 2 3 6 6 6 6 3; 4 1 2 3 6 6 6 6 3; 4 4 4 4 4 4 6 6 6; 6 4 4 4 4 4 4 1 6; 5 2 2 5 3 6 2 1 1; 5 5 5 2 5 2 2 1 1; 5 3 3 5 5 2 2 2 5; 6 6 6 5 5 4 2 2 4], [6 6 6 5 5 5 5 1 1; 1 1 3 4 5 5 5 5 1; 5 4 3 3 5 5 4 4 4; 5 5 5 5 5 4 4 4 4; 5 5 5 5 5 4 4 4 4; 3 5 5 5 6 4 4 4 4; 3 3 3 3 3 3 4 1 2; 3 3 3 3 3 3 4 5 5; 3 3 3 3 3 3 5 5 5], [6 6 4 2 2 3 2 2 6; 3 2 4 2 3 1 5 6 6; 3 3 2 2 3 3 3 3 3; 5 6 2 2 6 3 4 3 3; 3 4 2 5 6 6 6 5 3; 5 5 5 6 6 6 6 3 2; 5 5 5 5 4 2 4 2 4; 5 3 4 4 4 3 4 4 3; 5 4 4 4 4 4 2 3 3]],
    #                             [[4 3 3 3 2 2 2 2 2; 1 3 3 3 5 2 2 2 2; 3 3 3 3 5 5 5 2 2; 3 3 3 3 1 5 2 2 2; 5 6 6 3 1 5 2 2 2; 6 6 1 1 2 2 2 5 1; 6 6 1 1 2 2 2 1 5; 3 6 2 4 6 6 5 1 5; 6 6 3 5 6 3 3 1 5], [3 3 5 5 6 6 1 5 6; 5 1 1 1 5 6 1 6 6; 4 1 1 1 5 5 2 2 1; 1 1 1 1 5 6 6 2 1; 2 1 1 2 2 2 1 4 1; 1 3 5 2 2 2 1 1 1; 5 1 3 4 4 2 3 3 5; 5 5 5 5 3 3 3 3 2; 5 6 5 6 5 5 6 5 6], [4 4 4 4 5 3 3 3 6; 6 6 6 6 6 1 5 2 6; 6 3 2 6 6 1 1 2 2; 6 3 3 5 6 1 3 2 2; 5 3 3 3 3 1 2 2 3; 5 5 5 5 5 1 3 2 2; 5 5 5 4 3 4 4 4 2; 5 5 5 4 4 4 6 4 4; 5 5 4 4 4 4 2 4 4], [3 4 5 2 2 2 2 2 1; 3 5 6 6 2 6 2 2 6; 6 6 6 6 5 6 6 2 6; 5 3 5 3 1 4 4 5 1; 3 1 1 1 4 4 1 1 1; 3 1 1 3 3 4 6 6 4; 3 1 1 3 3 3 6 6 6; 1 1 1 1 1 2 4 6 2; 1 1 1 1 1 6 6 2 2], [2 2 1 6 1 2 2 3 3; 1 4 5 5 5 2 2 2 2; 1 4 5 5 4 2 4 4 4; 3 5 5 5 5 6 4 4 4; 3 6 6 4 5 4 4 4 4; 1 1 1 4 4 4 4 4 4; 1 3 4 4 4 4 4 4 4; 1 3 4 4 4 4 4 4 4; 5 4 4 5 4 4 4 4 4], [1 6 4 5 6 5 5 1 1; 6 4 4 4 3 3 5 5 5; 6 6 6 6 6 6 5 5 1; 6 6 6 6 6 3 5 5 6; 3 5 5 6 6 6 2 2 6; 2 1 2 2 3 6 2 2 3; 2 1 2 2 3 3 3 3 3; 1 1 6 5 3 3 3 3 3; 2 4 1 4 4 3 3 3 3]]]
            

    starting_configurations = [ [[1 1 6 6 1 1 6 6 6; 1 3 3 3 1 1 6 6 6; 4 3 3 3 3 4 4 6 6; 4 3 3 3 4 4 4 5 5; 4 3 1 1 1 3 3 4 4; 4 4 1 1 1 1 1 1 4; 1 1 2 2 6 1 1 1 4; 1 1 2 2 6 6 6 1 4; 1 1 2 2 2 1 3 3 1], [2 2 2 2 6 2 2 3 1; 2 2 2 2 2 2 2 2 1; 2 2 2 2 2 2 2 2 5; 1 3 3 2 2 2 2 6 6; 2 3 2 2 2 2 2 6 6; 2 2 2 2 2 2 6 6 6; 2 2 2 2 2 2 4 6 6; 2 2 2 2 2 2 4 3 3; 2 2 2 2 2 2 2 3 3], [6 2 4 3 3 3 1 5 4; 6 3 3 3 2 3 1 1 1; 3 3 3 3 3 3 6 1 1; 3 6 3 3 3 1 1 1 1; 6 6 6 3 3 1 1 1 1; 6 6 5 5 1 1 1 1 1; 3 6 5 2 1 1 1 1 1; 3 4 6 1 1 1 1 1 1; 3 1 1 6 1 1 1 4 2], [4 4 3 3 3 4 4 4 4; 3 4 4 4 4 4 4 4 4; 3 3 1 5 4 4 4 4 4; 3 3 5 3 4 4 4 4 4; 3 3 3 3 4 4 4 4 4; 3 3 3 3 4 4 4 4 6; 4 3 3 4 4 4 4 4 6; 3 3 4 4 4 4 4 4 4; 3 3 3 4 4 4 4 4 3], [5 5 5 4 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 1 5 5 5 5 5 6 5 5; 2 2 5 5 5 5 6 5 5; 2 2 5 5 5 6 6 6 5; 2 2 5 5 5 3 3 5 4; 2 2 5 5 5 3 5 5 4; 2 2 5 5 5 5 5 5 5], [4 6 3 3 3 3 1 6 6; 2 6 3 1 3 4 4 6 6; 3 3 6 6 6 6 6 6 6; 2 2 1 6 6 6 6 5 5; 6 1 5 6 6 6 1 5 5; 6 6 6 6 6 4 1 5 5; 6 5 6 6 4 4 1 1 5; 6 6 6 6 6 1 1 5 5; 6 6 6 6 1 1 1 5 5]],
                                [[1 1 1 6 5 1 1 1 1; 1 6 2 1 1 1 1 1 1; 2 2 2 5 5 1 1 1 1; 2 2 2 1 5 1 1 1 1; 3 3 1 1 1 1 1 5 1; 4 5 1 1 1 1 3 3 5; 4 1 1 1 1 3 3 3 3; 5 1 1 1 1 1 3 3 3; 5 6 1 1 1 1 3 3 4], [2 4 6 2 2 5 5 5 4; 2 2 2 2 2 2 2 6 4; 3 6 2 2 2 2 2 4 4; 1 1 6 2 2 2 2 4 4; 1 1 3 2 2 2 2 4 4; 1 1 3 3 1 2 2 4 2; 2 5 5 2 2 2 2 2 5; 2 2 1 2 2 2 2 2 4; 2 2 2 2 2 2 2 2 2], [3 4 6 3 2 3 1 3 3; 2 3 3 3 3 3 3 3 3; 3 3 3 1 3 3 3 3 3; 3 3 3 2 3 3 3 3 6; 3 3 3 3 3 3 1 1 1; 3 3 3 3 3 3 1 1 1; 3 3 6 5 3 6 1 1 6; 3 3 6 5 6 6 1 1 1; 3 3 2 6 6 3 1 1 1], [1 3 6 4 4 4 2 2 3; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 1; 4 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 3; 4 4 4 4 5 4 4 4 4; 2 5 5 6 6 2 5 1 3; 2 5 4 4 3 2 5 5 5], [5 5 5 5 5 5 5 5 5; 5 5 5 5 5 5 5 5 5; 6 2 5 5 5 5 5 5 5; 6 2 5 5 5 5 5 5 5; 6 2 2 2 5 5 5 5 5; 6 6 2 5 5 5 5 5 5; 6 6 6 6 4 5 5 5 5; 6 6 6 5 5 5 1 4 4; 6 6 6 6 5 5 1 1 4], [6 6 3 3 3 3 3 6 6; 6 6 6 3 3 3 3 2 6; 6 6 6 6 6 3 3 6 2; 6 6 6 6 6 6 6 6 6; 2 2 6 6 6 6 6 6 6; 2 2 1 6 6 6 6 6 4; 2 6 1 1 6 6 6 4 4; 2 5 5 6 6 6 2 4 4; 6 1 1 2 6 5 5 6 4]],
                                [[5 5 5 5 5 5 1 1 5; 5 5 5 5 5 5 5 5 5; 5 6 5 5 5 5 5 5 5; 6 6 6 1 1 5 5 5 5; 6 6 6 3 1 5 5 5 5; 6 6 6 3 1 5 5 5 6; 6 6 6 1 1 5 5 5 1; 1 2 6 1 5 5 5 5 5; 5 5 5 5 5 5 5 5 5], [4 6 4 2 2 2 2 2 2; 4 2 2 2 2 2 2 2 2; 5 5 2 2 2 2 2 6 5; 1 1 2 2 2 2 2 3 5; 1 1 2 2 2 2 2 6 5; 1 4 1 2 2 2 5 1 5; 1 6 1 2 2 2 2 3 6; 1 6 2 2 2 2 2 4 1; 1 2 2 2 2 2 2 4 1], [1 1 1 4 4 4 4 1 1; 2 1 1 6 1 3 4 4 3; 4 2 6 6 3 3 3 3 3; 4 4 3 3 3 3 3 3 3; 4 4 4 4 3 3 3 3 3; 4 4 4 4 3 3 3 2 3; 2 4 4 4 3 3 3 3 3; 2 3 1 6 3 3 3 3 3; 3 3 3 3 3 3 3 3 3], [3 3 3 1 1 1 1 3 3; 1 1 1 1 1 1 1 1 3; 1 1 1 1 1 1 1 1 1; 1 3 1 1 1 1 1 1 1; 1 1 1 1 4 5 4 4 1; 1 1 1 1 4 4 4 4 1; 1 1 1 1 5 3 4 6 6; 1 1 1 6 3 6 4 6 6; 4 4 5 5 6 6 3 6 2], [6 6 2 2 2 2 2 5 2; 3 3 3 2 2 2 2 2 2; 3 3 3 2 1 2 2 2 4; 3 3 3 6 6 6 6 2 2; 3 3 3 5 5 6 6 5 2; 3 3 3 5 5 6 6 5 2; 6 3 3 5 5 5 5 5 2; 6 6 4 5 6 5 5 5 2; 6 6 6 6 4 6 6 5 2], [4 2 2 3 3 3 3 6 6; 4 4 4 4 4 3 3 3 4; 4 4 4 4 4 6 6 6 6; 4 4 4 4 4 5 6 6 6; 4 4 4 4 6 6 6 6 6; 4 4 4 4 6 6 6 6 6; 4 4 4 4 6 4 6 6 6; 4 4 4 4 2 1 2 6 6; 4 4 4 4 6 4 4 4 6]]

                                ]


    for (index,starting_configuration) in pairs(starting_configurations)
        cube = RubiksCube(9)
        cube.configuration = starting_configuration

        time_experiment(cube, "funnel_test_E=$(E)_T=$(T)_trial_$index", 9, [p_swap], T, N_t; normalization=normalization)

    end

end