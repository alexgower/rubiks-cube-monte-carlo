using DelimitedFiles
using Plots
using LaTeXStrings
using StatsBase
using CSV
using DataFrames

using LsqFit

using Plots.PlotMeasures


include("../core/rubiks_cube.jl")
include("../core/swap_moves.jl")






function relaxed_anneal_graphs_plotter(simulation_name::String, swap_move_probabilities::Vector{Float64}; inherent_disorder::Bool=false, smooth_window::Int64=1, T_on::Float64=0.0, T_off::Float64=0.0)

    ##### ----- IMPORT RELAXED ANNEAL DATA -----

    ### --- SET UP DEFAULT PARAMETERS ---
    header_line = readlines(joinpath("results/relaxed_anneal_results/"*simulation_name*"_"*string(swap_move_probabilities[end])))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    temperature_vector::Vector{Float64} = []
    slice_rotation_temperature_vector::Vector{Float64} = []
    array_normalised_E_average_by_temperature::Vector{Vector{Float64}} = []
    array_normalised_standard_deviations_by_temperature::Vector{Vector{Float64}} = []
    array_specific_heat_capacities_by_temperature::Vector{Vector{Float64}} = []

    for (index,swap_move_probability) in pairs(swap_move_probabilities)
        simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)

        ## -- READ IN THE DATA --
        data_matrix = readdlm(joinpath("results/relaxed_anneal_results",simulation_name_to_use), ',', Float64, '\n', skipstart=2, comments=true, comment_char='#')

        this_temperature_vector = copy(data_matrix[:,1])
        if swap_move_probability == 1.0
            temperature_vector = this_temperature_vector
        elseif swap_move_probability == 0.0
            slice_rotation_temperature_vector = this_temperature_vector
        end

        E_average_by_temperature = copy(data_matrix[:,2]) 
        E_squared_average_by_temperature = copy(data_matrix[:,4]) 

        normalization_energy = solved_configuration_energy(cube)


        ## -- WINDOW SMOOTHENING PROCESS ---
        E_average_by_temperature = window_average(E_average_by_temperature, smooth_window)
        E_squared_average_by_temperature = window_average(E_squared_average_by_temperature, smooth_window)

        ## -- CALCULATE DESIRED QUANTITIES --
        push!(array_normalised_E_average_by_temperature, -E_average_by_temperature ./ normalization_energy)
        push!(array_normalised_standard_deviations_by_temperature, sqrt.(E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ normalization_energy)
        push!(array_specific_heat_capacities_by_temperature, (E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ this_temperature_vector.^2)

    end

    
    ## -- GENERAL PARAMETERS --
    temperatures = reverse(temperature_vector)
    N_T = length(temperature_vector)-1
    swap_move_one_probabilty_index = findfirst(swap_move_probabilities .== 1.0)
    swap_move_zero_probabilty_index = findfirst(swap_move_probabilities .== 0.0)



    ### --- PLOT RELAXED ANNEAL GRAPHS ---  
    array_renormalised_E_average_by_temperature::Vector{Vector{Float64}} = []
    for vector in array_normalised_E_average_by_temperature
        push!(array_renormalised_E_average_by_temperature, vector .+ 1.0)
    end

    E_star = array_renormalised_E_average_by_temperature[swap_move_zero_probabilty_index][end]
    println("E_star = ", E_star)

    # mean_graph = plot(temperature_vector, array_renormalised_E_average_by_temperature, xlabel="Temperature", ylabel="Average Energy, E", title="L=$L Rubik's Cube Anneal", labels=reshape(["Slice Rotation Cube", "Swap Move Cube"], 1, length(swap_move_probabilities)))
    mean_graph = plot(temperature_vector, array_renormalised_E_average_by_temperature[swap_move_one_probabilty_index], xlabel="Temperature", ylabel="Average Energy, E", title="L=$L Rubik's Cube Anneal", label="Swap Move Cube", color=:blue, margin=3mm)
    plot!(mean_graph, slice_rotation_temperature_vector, array_renormalised_E_average_by_temperature[swap_move_zero_probabilty_index], label="Slice Rotation Cube", color=:red)

    
    # mean_std_graph = plot(temperature_vector, array_renormalised_E_average_by_temperature, yerr=transpose(array_normalised_standard_deviations_by_temperature), markerstrokecolor=:auto, xlabel="Temperature", ylabel="Energy, E", title="L=$L Rubik's Cube Anneal", labels=reshape(["Slice Rotation Cube", "Swap Move Cube"], 1, length(swap_move_probabilities)))
    mean_std_graph = plot(temperature_vector, array_renormalised_E_average_by_temperature[swap_move_one_probabilty_index], yerr=transpose(array_normalised_standard_deviations_by_temperature[swap_move_one_probabilty_index]), xlabel="Temperature", ylabel="Energy, E", title="L=$L Rubik's Cube Anneal", label="Swap Move Cube", color=:blue)
    plot!(mean_std_graph, slice_rotation_temperature_vector, array_renormalised_E_average_by_temperature[swap_move_zero_probabilty_index], yerr=transpose(array_normalised_standard_deviations_by_temperature[swap_move_zero_probabilty_index]), label="Slice Rotation Cube", color=:red)

    if !inherent_disorder
        hline!(mean_std_graph, [-1.0], linestyle=:dash, color=:black, label="")
        hline!(mean_graph, [-1.0], linestyle=:dash, color=:black, label="")
    end

    if T_on != 0.0
        vline!(mean_std_graph, [T_on], linestyle=:dot, color=:orange, label="")
        vline!(mean_graph, [T_on], linestyle=:dot, color=:orange, label="")
        annotate!(mean_std_graph, [(T_on+0.2, ylims(mean_std_graph)[1]-0.025, Plots.text(L"T^*_{on}", 8, :black))])
        annotate!(mean_graph, [(T_on+0.2, ylims(mean_graph)[1]-0.025, Plots.text(L"T^*_{on}", 8, :black))])
    end
    
    if T_off != 0.0
        vline!(mean_std_graph, [T_off], linestyle=:dot, color=:red, label="")
        vline!(mean_graph, [T_off], linestyle=:dot, color=:red, label="")
        annotate!(mean_std_graph, [(T_off-0.2, ylims(mean_std_graph)[1]-0.025, Plots.text(L"T^*_{off}", 8, :black))])
        annotate!(mean_graph, [(T_off-0.2, ylims(mean_graph)[1]-0.025, Plots.text(L"T^*_{off}", 8, :black))])
    end
    
    if E_star != 0.0
        hline!(mean_std_graph, [E_star], linestyle=:dot, color=:green, label="")
        hline!(mean_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(mean_std_graph, [(xlims(mean_std_graph)[1]+0.25, E_star+0.02, Plots.text(L"E^*", 8, :black))])
        annotate!(mean_graph, [(xlims(mean_graph)[1]+0.25, E_star+0.02, Plots.text(L"E^*", 8, :black))])
    end



    # Plot infinite temperature energy
    # hline!(mean_std_graph, [0.16666666666666666], linestyle=:dash, color=:black, label="")
    # hline!(mean_graph, [0.16666666666666666], linestyle=:dash, color=:black, label="")

    savefig(mean_graph, "results/relaxed_anneal_results/$(simulation_name)_mean.png")
    savefig(mean_graph, "results/relaxed_anneal_results/$(simulation_name)_mean.svg")

    savefig(mean_std_graph, "results/relaxed_anneal_results/$(simulation_name)_mean_std.png")
    savefig(mean_std_graph, "results/relaxed_anneal_results/$(simulation_name)_mean_std.svg")


   


    ### --- PLOT SPECIFIC HEAT CAPACITY GRAPH ---
    specific_heat_capacity_graph = plot(temperature_vector, array_specific_heat_capacities_by_temperature[swap_move_one_probabilty_index], xlabel="Temperature", ylabel="Specific Heat Capacity", title="Specific Heat Capacity", label="Swap Move Cube", color=:blue)
    plot!(specific_heat_capacity_graph, slice_rotation_temperature_vector, array_specific_heat_capacities_by_temperature[swap_move_zero_probabilty_index], label="Slice Rotation Cube", color=:red)    
    savefig(specific_heat_capacity_graph, "results/relaxed_anneal_results/$(simulation_name)_specific_heat_capacity.png")
    savefig(specific_heat_capacity_graph, "results/relaxed_anneal_results/$(simulation_name)_specific_heat_capacity.svg")


    ########################################################################################################################################################################################################################################################################



    ### --- PLOT ENTROPY/TEMPERATURE GRAPHS ---



    heat_capacities_by_temperature = reverse(array_specific_heat_capacities_by_temperature[swap_move_one_probabilty_index])
    normalized_energies_by_temperature = reverse(array_normalised_E_average_by_temperature[swap_move_one_probabilty_index])
    absolute_energies_by_temperature = (normalized_energies_by_temperature .+ 1.0)*-solved_configuration_energy(cube)


    ## -- GENERATE S(T=0) FIXED ENTROPY/TEMPERATURE GRAPHS ---

    entropy_by_temperature = zeros(N_T+1)
    alternative_entropy_by_temperature = zeros(N_T+1)

    for temperature_index in 2:N_T+1
        println("Trapezium Width %: ", 0.5 * (temperatures[temperature_index] - temperatures[temperature_index-1])/(temperatures[end] - temperatures[1]))
        entropy_by_temperature[temperature_index] = entropy_by_temperature[temperature_index-1] + 0.5 * (temperatures[temperature_index] - temperatures[temperature_index-1]) * (heat_capacities_by_temperature[temperature_index]/temperatures[temperature_index] + heat_capacities_by_temperature[temperature_index-1]/temperatures[temperature_index-1])
    end

    for temperature_index in 2:N_T+1
        println("Alternative Trapezium Width: % ", 0.5 * (absolute_energies_by_temperature[temperature_index] - absolute_energies_by_temperature[temperature_index-1])/(absolute_energies_by_temperature[end] - absolute_energies_by_temperature[1]))
        alternative_entropy_by_temperature[temperature_index] = alternative_entropy_by_temperature[temperature_index-1] + 0.5 * (absolute_energies_by_temperature[temperature_index] - absolute_energies_by_temperature[temperature_index-1]) * (1/temperatures[temperature_index] + 1/temperatures[temperature_index-1])
    end

    entropy_by_temperature_graph = plot(temperatures, entropy_by_temperature, xlabel="Temperature", ylabel="Entropy, S", title="L=$L Rubik's Cube Anneal", label="Entropy", color=:black)
    plot!(entropy_by_temperature_graph, temperatures, alternative_entropy_by_temperature, label="Alternative Entropy", color=:blue)




    # ## -- GENERATE S(T=INFINITY) FIXED ENTROPY/TEMEPRATURE GRAPHS --

    # upper_fixed_entropy_by_temperature = zeros(N_T+1)
    # upper_fixed_alternative_entropy_by_temperature = zeros(N_T+1)

    # upper_entropies = zeros(11)
    # # From 'Art of Problem Solving' post
    # # upper_entropies[5] = 282870942277741856536180333107150328293127731985672134721536000000000000000
    # # upper_entropies[7] = 19500551183731307835329126754019748794904992692043434567152132912323232706135469180065278712755853360682328551719137311299993600000000000000000000000000000000000
    # # upper_entropies[9] = 14170392390542612915246393916889970752732946384514830589276833655387444667609821068034079045039617216635075219765012566330942990302517903971787699783519265329288048603083134861573075573092224082416866010882486829056000000000000000000000000000000000000000000000000000000000000000
    # # From Salkinder paper
    # upper_entropies[5] = 2582636272886959379162819698174683585918088940054237132144778804568925405184000000000000000
    # upper_entropies[7] = 14841288742494293811131211592812790126449902279909457255335745506781477789026449207045807257330619800316538795293665546204430425928257835670695234259221271924781093824652902400000000000000000000000000000000000
    # upper_entropies[9] = 8207886298937235484582402588758685802310028366167651665616401401542874344822914074802549653732095854271871596416146268508634462875956312770239894329965193325658757918817194281459237701144924933414816036003138370637933027571286472248044977240412593866789491685980731724648262702614573108239398605603042635022336000000000000000000000000000000000000000000000000000000000000000
    
    # physically_indistinguishable_factors = [(24^6/2)^((L^2-4*L+3)/4) for L in 1:11] # Note this factor applies for odd-L cubes swap_move_one_probabilty_index

    # upper_entropies =  log.(upper_entropies ./ physically_indistinguishable_factors)
    
    # upper_entropies[9] = 638.16464041131606054869399362152095974484188834809571212426469616 # This numbed needed WolframAlpha as is too big
    
    # upper_fixed_entropy_by_temperature[end] = upper_entropies[L]
    # upper_fixed_alternative_entropy_by_temperature[end] = upper_entropies[L]


    # for temperature_index in N_T:-1:1
    #     upper_fixed_entropy_by_temperature[temperature_index] += upper_fixed_entropy_by_temperature[temperature_index+1]
    #     upper_fixed_entropy_by_temperature[temperature_index] += 0.5 * (temperatures[temperature_index] - temperatures[temperature_index+1]) * (heat_capacities_by_temperature[temperature_index]/temperatures[temperature_index] + heat_capacities_by_temperature[temperature_index+1]/temperatures[temperature_index+1])
    # end

    # for temperature_index in N_T:-1:1
    #     upper_fixed_alternative_entropy_by_temperature[temperature_index] += upper_fixed_alternative_entropy_by_temperature[temperature_index+1]
    #     upper_fixed_alternative_entropy_by_temperature[temperature_index] += 0.5 * (absolute_energies_by_temperature[temperature_index] - absolute_energies_by_temperature[temperature_index+1]) * (1/temperatures[temperature_index] + 1/temperatures[temperature_index+1]) 
    # end

    # plot!(entropy_by_temperature_graph, temperatures, upper_fixed_entropy_by_temperature, label="Upper Fixed Entropy", color=:red)
    # plot!(entropy_by_temperature_graph, temperatures, upper_fixed_alternative_entropy_by_temperature, label="Alternative Upper Fixed Entropy", color=:green)




    savefig(entropy_by_temperature_graph, "results/relaxed_anneal_results/$(simulation_name)_entropy_temperature.png")
    savefig(entropy_by_temperature_graph, "results/relaxed_anneal_results/$(simulation_name)_entropy_temperature.svg")

    

    ########################################################################################################################################################################################################################################################################








    ### --- PLOT SLICE/SWAP CONNECTIVITIES AND MINIMA/SADDLE CONFIGURAITONS GRAPHS ---

    this_entropy_by_temperature = entropy_by_temperature

    z_slice = isodd(L) ? 6*(L-1) : 6*L
    z_swap::BigInt = total_number_of_swap_moves(cube)

    saddle_connectivity_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), this_entropy_by_temperature, xlabel=L"Energy, E", ylabel="ln(N(E))", title="Number of Saddle Connectivities"; label="Total Configurations, ln(N(E))", color=:black, legend=:topleft)
    saddle_configuration_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), this_entropy_by_temperature, xlabel="Energy, E", ylabel="ln(N(E))", title="Number of Saddle Configurations, L=$L Cube"; label="Total Configurations, ln(N(E))", color=:black, legend=:topleft)
    saddle_configuration_from_even_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), this_entropy_by_temperature, xlabel="Energy, E", ylabel="ln(N(E))", title="Number of Saddle Configurations, L=$L Cube"; label="Total Configurations, ln(N(E))", color=:black, legend=:topleft)
    minima_configuration_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), this_entropy_by_temperature, xlabel="Energy, E", ylabel="ln(N(E))", title="Number of Minima Configurations, L=$L Cube"; label="Total Configurations, ln(N(E))", color=:black, legend=:topleft)
    minima_configuration_from_even_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), this_entropy_by_temperature, xlabel="Energy, E", ylabel="ln(N(E))", title="Number of Minima Configurations, L=$L Cube"; label="Total Configurations, ln(N(E))", color=:black, legend=:topleft)

    # Add analytic model lines
    # slice_rotation_entropy_scaling_gradient = log(6(L-1))*(1/(2*L))
    # plot!(slice_saddle_connectivity_graph, [minimum(absolute_energies_by_temperature[2:end]), maximum(absolute_energies_by_temperature[2:end])], [0, slice_rotation_entropy_scaling_gradient*(maximum(absolute_energies_by_temperature[2:end])-minimum(absolute_energies_by_temperature[2:end]))], line=:dash, color=:blue, lw=2, label=L"\Delta E = 2L"*" (Mode) Layers")
    # slice_rotation_entropy_scaling_gradient = log(6(L-1))*(1/(8*L))
    # plot!(slice_saddle_connectivity_graph, [minimum(absolute_energies_by_temperature[2:end]), maximum(absolute_energies_by_temperature[2:end])], [0, slice_rotation_entropy_scaling_gradient*(maximum(absolute_energies_by_temperature[2:end])-minimum(absolute_energies_by_temperature[2:end]))], line=:dash, color=:green, lw=2, label=L"\Delta E = 8L"*" (Max) Layers")



    # –– RECONSTRUCT (E_0, E_1) CONNECTIVITY HISTOGRAM --
    normalised_E0_E1_histogram_swap = reconstruct_histogram(joinpath("results/relaxed_anneal_results",simulation_name*"_normalised_histogram_data_swap.csv"))
    normalised_E0_E1_histogram_slice = reconstruct_histogram(joinpath("results/relaxed_anneal_results",simulation_name*"_normalised_histogram_data_slice.csv"))


    ## - DO SIMPLE ENERGY INCREASING/DECREASING PROPORTIONS CALCULATIONS AND GRAPHS --
    proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice = zeros(length(normalised_E0_E1_histogram_slice.edges[1])-1)
    proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap = zeros(length(normalised_E0_E1_histogram_slice.edges[1])-1)

    for E0_bin_index in 1:length(normalised_E0_E1_histogram_slice.edges[1])-1
        for E_bin_index in 1:E0_bin_index-1
            proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice[E0_bin_index] += normalised_E0_E1_histogram_slice.weights[E0_bin_index, E_bin_index]
        end
    end

    for E0_bin_index in 1:length(normalised_E0_E1_histogram_swap.edges[1])-1
        for E_bin_index in 1:E0_bin_index-1
            proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap[E0_bin_index] += normalised_E0_E1_histogram_swap.weights[E0_bin_index, E_bin_index]
        end
    end

    E0_values = normalised_E0_E1_histogram_swap.edges[1][1:end-1]

    # Now only keep the proportion values that have non-zero swap move energy decreasing proportion (to exclude energies not reachable in this cube)
    # Also remove all indices where proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap==0.0 and energy is bigger than half of the maximum energy
    valid_E0_indices = filter(i -> ((
        proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap[i] != 0.0) &&
        !(proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice[i] == 0.0 && E0_values[i] > maximum(E0_values) / 2 )
    ), 1:length(E0_values))
    proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap = [proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap[i] for i in valid_E0_indices]
    proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice = [proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice[i] for i in valid_E0_indices]
    E0_values = [E0_values[i] for i in valid_E0_indices]


    energy_decreasing_connections_proportion_graph = plot(E0_values./-solved_configuration_energy(cube), proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice, xlabel="Energy, E", ylabel="Proportion of Energy Decreasing Connections", title="L=$L Cube"; label="Slice Rotation Cube", color=:red, legend=:topleft)
    plot!(energy_decreasing_connections_proportion_graph, E0_values./-solved_configuration_energy(cube), proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap, label="Swap Move Cube", color=:blue)

    energy_increasing_connections_proportion_graph = plot(E0_values./-solved_configuration_energy(cube), 1 .- proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice, xlabel="Energy, E", ylabel="Proportion of Energy Increasing Connections", title="L=$L Cube"; label="Slice Rotation Cube", color=:red, legend=:bottomleft)
    plot!(energy_increasing_connections_proportion_graph, E0_values./-solved_configuration_energy(cube), 1 .- proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap, label="Swap Move Cube", color=:blue)

    # Overplot where slice energy decreasing proportion is =0.0 in green
    entropic_wall_indices = findall(==(0.0), proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice)
    entropic_wall_E0_values = [E0_values[i] for i in entropic_wall_indices]./-solved_configuration_energy(cube)
    plot!(energy_decreasing_connections_proportion_graph, [entropic_wall_E0_values], [proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice[entropic_wall_indices]], color=:green, label="")
    plot!(energy_increasing_connections_proportion_graph, [entropic_wall_E0_values], [1 .- proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice[entropic_wall_indices]], color=:green, label="")

    if E_star != 0.0
        vline!(energy_decreasing_connections_proportion_graph, [E_star], linestyle=:dot, color=:green, label="")
        vline!(energy_increasing_connections_proportion_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(energy_decreasing_connections_proportion_graph, [(E_star+0.02,ylims(energy_decreasing_connections_proportion_graph)[1]-0.02, Plots.text(L"E^*", 8, :black))])
        annotate!(energy_increasing_connections_proportion_graph, [(E_star+0.02, ylims(energy_increasing_connections_proportion_graph)[1]-0.02, Plots.text(L"E^*", 8, :black))])
    end

    savefig(energy_decreasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_decreasing_connections_proportion.png")
    savefig(energy_decreasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_decreasing_connections_proportion.svg")

    savefig(energy_increasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_increasing_connections_proportion.png")
    savefig(energy_increasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_increasing_connections_proportion.svg")


    # Fit exponential to slice energy decreasing proportion
    model(x, p) = p[1] * exp.(p[2] * x)
    fit_exponential = curve_fit(model, E0_values, proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_slice, [0.0, 0.0])
    println("Slice Energy Decreasing Proportion Exponential Fit: ", fit_exponential.param)
    plot!(energy_decreasing_connections_proportion_graph, E0_values./-solved_configuration_energy(cube), fit_exponential.param[1]*exp.(fit_exponential.param[2]*E0_values), label="$(fit_exponential.param[1])exp($(fit_exponential.param[2])E)", color=:black, linestyle=:dash)

    # Fit exponential to swap energy decreasing proportion
    fit_exponential = curve_fit(model, E0_values, proportion_of_energy_decreasing_connections_from_even_by_absolute_energies_swap, [0.0, 0.0])
    println("Swap Energy Decreasing Proportion Exponential Fit: ", fit_exponential.param)
    plot!(energy_decreasing_connections_proportion_graph, E0_values./-solved_configuration_energy(cube), fit_exponential.param[1]*exp.(fit_exponential.param[2]*E0_values), label="$(fit_exponential.param[1])exp($(fit_exponential.param[2])E)", color=:pink, linestyle=:dash)

    # Save the graph
    savefig(energy_decreasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_decreasing_connections_proportion_fit.png")
    savefig(energy_decreasing_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_energy_decreasing_connections_proportion_fit.svg")
































    ## -- DO SLICE/SWAP CONNECTIVITIES CALCULATIONS --
    downwards_swap_connections_from_below_by_absolute_energies = zeros(N_T+1)
    downwards_slice_connections_from_below_by_absolute_energies  = zeros(N_T+1)

    downwards_swap_connections_from_even_by_absolute_energies = zeros(N_T+1)
    downwards_slice_connections_from_even_by_absolute_energies = zeros(N_T+1)

    downwards_swap_connections_proportion_from_below_by_absolute_energies = zeros(N_T+1)
    downwards_slice_connections_proportion_from_below_by_absolute_energies = zeros(N_T+1)

    downwards_swap_connections_proportion_from_even_by_absolute_energies = zeros(N_T+1)
    downwards_slice_connections_proportion_from_even_by_absolute_energies = zeros(N_T+1)


    for (E_absolute_energy_index,absolute_energy) in pairs(absolute_energies_by_temperature)

        E_bin_index = findfirst(>=(absolute_energy), normalised_E0_E1_histogram_swap.edges[2])
        configurations_of_E = exp(this_entropy_by_temperature[E_absolute_energy_index])

        for E0_bin_index in 1:E_bin_index

            E0_absolute_energy_index = findfirst(>=(normalised_E0_E1_histogram_swap.edges[1][E0_bin_index]), absolute_energies_by_temperature)
            configurations_of_E0 = !isnothing(E0_absolute_energy_index) ? exp(this_entropy_by_temperature[E0_absolute_energy_index]) : exp(this_entropy_by_temperature[end])

            # CALCULATIONS FROM BELOW
            proportion_of_E0_swap_connections_to_E = normalised_E0_E1_histogram_swap.weights[E0_bin_index, E_bin_index]
            proportion_of_E0_slice_connections_to_E = normalised_E0_E1_histogram_slice.weights[E0_bin_index, E_bin_index]

            downwards_swap_connections_from_below_by_absolute_energies[E_absolute_energy_index] += proportion_of_E0_swap_connections_to_E * configurations_of_E0 * z_swap
            downwards_slice_connections_from_below_by_absolute_energies[E_absolute_energy_index] += proportion_of_E0_slice_connections_to_E * configurations_of_E0 * z_slice

            # CALCULATIONS FROM EVEN
            proportion_of_E_swap_connections_to_E0 = normalised_E0_E1_histogram_swap.weights[E_bin_index, E0_bin_index]
            proportion_of_E_slice_connections_to_E0 = normalised_E0_E1_histogram_slice.weights[E_bin_index, E0_bin_index]

            downwards_swap_connections_from_even_by_absolute_energies[E_absolute_energy_index] += proportion_of_E_swap_connections_to_E0 * configurations_of_E * z_swap
            downwards_slice_connections_from_even_by_absolute_energies[E_absolute_energy_index] += proportion_of_E_slice_connections_to_E0 * configurations_of_E * z_slice

            downwards_swap_connections_proportion_from_even_by_absolute_energies[E_absolute_energy_index] += proportion_of_E_swap_connections_to_E0
            downwards_slice_connections_proportion_from_even_by_absolute_energies[E_absolute_energy_index] += proportion_of_E_slice_connections_to_E0

        end

    end

    downwards_swap_connections_proportion_from_below_by_absolute_energies = [x/(configurations_of_E*z_swap) for (x, configurations_of_E) in zip(downwards_swap_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]
    downwards_slice_connections_proportion_from_below_by_absolute_energies = [x/(configurations_of_E*z_slice) for (x, configurations_of_E) in zip(downwards_slice_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]


    ## -- DO SADDLE/MINIMA CONFIGURATIONS CALCULATIONS --
    number_swap_minima_configurations_from_below_by_absolute_energies = [x < configurations_of_E ? configurations_of_E - x : 0.0 for (x, configurations_of_E) in zip(downwards_swap_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_slice_minima_configurations_from_below_by_absolute_energies = [x < configurations_of_E ? configurations_of_E - x : 0.0 for (x, configurations_of_E) in zip(downwards_slice_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_swap_minima_configurations_from_even_by_absolute_energies = [x < configurations_of_E ? configurations_of_E - x : 0.0 for (x, configurations_of_E) in zip(downwards_swap_connections_from_even_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_slice_minima_configurations_from_even_by_absolute_energies = [x < configurations_of_E ? configurations_of_E - x : 0.0 for (x, configurations_of_E) in zip(downwards_slice_connections_from_even_by_absolute_energies,exp.(this_entropy_by_temperature))]

    number_swap_saddle_configurations_from_below_by_absolute_energies = [x < configurations_of_E ? x : configurations_of_E for (x, configurations_of_E) in zip(downwards_swap_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_slice_saddle_configurations_from_below_by_absolute_energies = [x < configurations_of_E ? x : configurations_of_E for (x, configurations_of_E) in zip(downwards_slice_connections_from_below_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_swap_saddle_configurations_from_even_by_absolute_energies = [x < configurations_of_E ? x : configurations_of_E for (x, configurations_of_E) in zip(downwards_swap_connections_from_even_by_absolute_energies,exp.(this_entropy_by_temperature))]
    number_slice_saddle_configurations_from_even_by_absolute_energies = [x < configurations_of_E ? x : configurations_of_E for (x, configurations_of_E) in zip(downwards_slice_connections_from_even_by_absolute_energies,exp.(this_entropy_by_temperature))]




    ## -- PRINT RAW DATA RESULTS -- 
    for i in 1:length(downwards_swap_connections_from_below_by_absolute_energies)
        # println("Swap Downward Connections at $(absolute_energies_by_temperature[i]): ", downwards_swap_connections_from_below_by_absolute_energies[i])
        # println("Total Swap Connections at $(absolute_energies_by_temperature[i]): ", exp(this_entropy_by_temperature[i])*z_swap)
        # println("--")
        # println("Slice Saddle Configurations at $(absolute_energies_by_temperature[i]): ", downwards_slice_connections_from_below_by_absolute_energies[i])
        # println("Total Slice Connections at $(absolute_energies_by_temperature[i]): ", exp(this_entropy_by_temperature[i])*z_slice)
        # println("--")
        println("FROM EVEN - Swap Saddle Configurations at $(absolute_energies_by_temperature[i]): ", number_swap_saddle_configurations_from_even_by_absolute_energies[i])
        println("FROM EVEN - Slice Saddle Configurations at $(absolute_energies_by_temperature[i]): ", number_slice_saddle_configurations_from_even_by_absolute_energies[i])
        println("FROM EVEN - Swap Minima Configurations at $(absolute_energies_by_temperature[i]): ", number_swap_minima_configurations_from_even_by_absolute_energies[i])
        println("FROM EVEN - Slice Minima Configurations at $(absolute_energies_by_temperature[i]): ", number_slice_minima_configurations_from_even_by_absolute_energies[i])
        println("FROM BELOW - Swap Saddle Configurations at $(absolute_energies_by_temperature[i]): ", number_swap_saddle_configurations_from_below_by_absolute_energies[i])
        println("FROM BELOW - Slice Saddle Configurations at $(absolute_energies_by_temperature[i]): ", number_slice_saddle_configurations_from_below_by_absolute_energies[i])
        println("FROM BELOW - Swap Minima Configurations at $(absolute_energies_by_temperature[i]): ", number_swap_minima_configurations_from_below_by_absolute_energies[i])
        println("FROM BELOW - Slice Minima Configurations at $(absolute_energies_by_temperature[i]): ", number_slice_minima_configurations_from_below_by_absolute_energies[i])
        println("Total Configurations at $(absolute_energies_by_temperature[i]): ", exp(this_entropy_by_temperature[i]))
        println("")
    end


    ## -- PLOT GRAPHS --
    
    # Saddle Connectivity Graph
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(downwards_swap_connections_from_below_by_absolute_energies), label="Swap Connections From Below", color=:blue, lw=1)   
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(downwards_slice_connections_from_below_by_absolute_energies), label="Slice Connections From Below", color=:red, lw=1)

    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity.png")
    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity.svg")

    # Saddle Connectivity From Below + Total Connectivity Graph
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(exp.(this_entropy_by_temperature)*z_swap), label="Total Swap Connections", color=:purple, lw=1)
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(exp.(this_entropy_by_temperature)*z_slice), label="Total Slice Connections", color=:pink, lw=1)

    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity_total.png")
    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity_total.svg")

    # Saddle Connectivity + Total Connectivity + Even Connections Graph
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(downwards_swap_connections_from_even_by_absolute_energies), label="Slice Connections From Even", color=:green, seriestype=:scatter, lw=1)
    plot!(saddle_connectivity_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(downwards_slice_connections_from_even_by_absolute_energies), label="Swap Connections From Even", color=:orange, seriestype=:scatter, lw=1)
    
    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity_even_too.png")
    savefig(saddle_connectivity_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_connectivity_even_too.svg")




    # Number of Saddle Configurations From Below Graph
    plot!(saddle_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(number_swap_saddle_configurations_from_below_by_absolute_energies), label="Swap Saddle Configurations, ln(N⁻(E))", color=:blue, lw=1)
    plot!(saddle_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(number_slice_saddle_configurations_from_below_by_absolute_energies), label="Slice Saddle Configurations, ln(N⁻(E))", color=:red, lw=1)

    if E_star != 0.0
        vline!(saddle_configuration_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(saddle_configuration_graph, [(E_star+0.02, ylims(saddle_configuration_graph)[1]+40, Plots.text(L"E^*", 8, :black))])
    end

    savefig(saddle_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations.png")
    savefig(saddle_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations.svg")

    # Number of Saddle Configurations + Even Configurations Graph
    plot!(saddle_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(number_swap_saddle_configurations_from_even_by_absolute_energies), label="Swap Saddle Configurations From Even", color=:green, seriestype=:scatter, lw=1)
    plot!(saddle_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), log.(number_slice_saddle_configurations_from_even_by_absolute_energies), label="Slice Saddle Configurations From Even", color=:orange, seriestype=:scatter, lw=1)

    savefig(saddle_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations_even_too.png")
    savefig(saddle_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations_even_too.svg")

    # Number of Saddle Configurations From Even + Total Configurations Graph
    log_number_swap_saddle_configurations_from_even_by_absolute_energies = [isinf(log(x)) ? 0.0 : log(x) for x in number_swap_saddle_configurations_from_even_by_absolute_energies]
    log_number_slice_saddle_configurations_from_even_by_absolute_energies = [isinf(log(x)) ? 0.0 : log(x) for x in number_slice_saddle_configurations_from_even_by_absolute_energies]

    # Find last index where both log_number_swap_saddle_configurations_from_even_by_absolute_energies and log_number_slice_saddle_configurations_from_even_by_absolute_energies are non-zero
    last_non_zero_index = 0
    for i in eachindex(log_number_swap_saddle_configurations_from_even_by_absolute_energies)
        if log_number_swap_saddle_configurations_from_even_by_absolute_energies[i] != 0.0 && log_number_slice_saddle_configurations_from_even_by_absolute_energies[i] != 0.0
            last_non_zero_index = i
        end
    end


    plot!(saddle_configuration_from_even_graph, absolute_energies_by_temperature[1:last_non_zero_index]./-solved_configuration_energy(cube), log_number_swap_saddle_configurations_from_even_by_absolute_energies[1:last_non_zero_index], label="Swap Saddle Configurations, ln(N⁻(E))", color=:blue, lw=1)
    plot!(saddle_configuration_from_even_graph, absolute_energies_by_temperature[1:last_non_zero_index]./-solved_configuration_energy(cube), log_number_slice_saddle_configurations_from_even_by_absolute_energies[1:last_non_zero_index], label="Slice Saddle Configurations, ln(N⁻(E))", color=:red, lw=1)

    if E_star != 0.0
        vline!(saddle_configuration_from_even_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(saddle_configuration_from_even_graph, [(E_star+0.02, ylims(saddle_configuration_from_even_graph)[1]+40, Plots.text(L"E^*", 8, :black))])
    end

    savefig(saddle_configuration_from_even_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations_from_even.png")
    savefig(saddle_configuration_from_even_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_configurations_from_even.svg")



    # Saddle Proportion Graph
    saddle_proportion_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), [isinf(-log(number_swap_saddle_configurations_from_even_by_absolute_energies[i]) + this_entropy_by_temperature[i]) ? 0.0 : -log(number_swap_saddle_configurations_from_even_by_absolute_energies[i])+this_entropy_by_temperature[i] for i in 1:length(number_swap_saddle_configurations_from_even_by_absolute_energies)], xlabel="Energy, E", ylabel="-ln[N⁻(E)/N(E)]", title="Saddle Proportions for L=$L Cube"; label="Swap Move Cube", color=:blue, legend=:topright)
    plot!(saddle_proportion_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), -log.(number_slice_saddle_configurations_from_even_by_absolute_energies) + this_entropy_by_temperature, label="Slice Rotation Cube", color=:red)
    
    if E_star != 0.0
        vline!(saddle_proportion_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(saddle_proportion_graph, [(E_star+0.02, ylims(saddle_proportion_graph)[1]+0.2, Plots.text(L"E^*", 8, :black))])
    end

    for i in eachindex(number_slice_saddle_configurations_from_even_by_absolute_energies)
        println("Slice Saddle Proportion at $(absolute_energies_by_temperature[i]/-solved_configuration_energy(cube)): ", 1/exp(-log(number_slice_saddle_configurations_from_even_by_absolute_energies[i]) + this_entropy_by_temperature[i]))
    end

    savefig(saddle_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_proportion.png")
    savefig(saddle_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_saddle_proportion.svg")


    # Minima proportion graph
    swap_minima_negative_log_likelihood = -log.(number_swap_minima_configurations_from_even_by_absolute_energies) + this_entropy_by_temperature
    slice_minima_negative_log_likelihood = -log.(number_slice_minima_configurations_from_even_by_absolute_energies) + this_entropy_by_temperature
    minima_proportion_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), swap_minima_negative_log_likelihood , xlabel="Energy, E", ylabel="-ln[N⁺(E)/N(E)]", title="Minima Proportions for L=$L Cube"; label="Swap Move Cube", color=:blue, legend=:topleft)
    plot!(minima_proportion_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), slice_minima_negative_log_likelihood, label="Slice Rotation Cube", color=:red)

    for i in eachindex(slice_minima_negative_log_likelihood)
        println("Slice Minima Proportion at $(absolute_energies_by_temperature[i]/-solved_configuration_energy(cube)): ", exp(-slice_minima_negative_log_likelihood[i]))
    end

    if E_star != 0.0
        vline!(minima_proportion_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(minima_proportion_graph, [(E_star+0.02, ylims(minima_proportion_graph)[1]+0.1, Plots.text(L"E^*", 8, :black))])
    end

    savefig(minima_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_proportion.png")
    savefig(minima_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_proportion.svg")



    # Number of Minima Configurations Graph
    swap_minima_configuration_entropy_from_below_by_absolute_energies = [x==0 ? 0 : log(x) for x in number_swap_minima_configurations_from_below_by_absolute_energies]
    slice_minima_configuration_entropy_from_below_by_absolute_energies = [x==0 ? 0 : log(x) for x in number_slice_minima_configurations_from_below_by_absolute_energies]
    plot!(minima_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), swap_minima_configuration_entropy_from_below_by_absolute_energies, label="Swap Minima Configurations, ln(N⁺(E))", color=:blue, lw=1)
    plot!(minima_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), slice_minima_configuration_entropy_from_below_by_absolute_energies, label="Slice Minima Configurations, ln(N⁺(E))", color=:red, lw=1)

    if E_star != 0.0
        vline!(minima_configuration_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(minima_configuration_graph, [(E_star+0.01, ylims(minima_configuration_graph)[1]-0.05, Plots.text(L"E^*", 8, :black))])
    end

    savefig(minima_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations.png")
    savefig(minima_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations.svg")

    # Number of Minima Configurations + Even Configurations Graph
    swap_minima_configuration_entropy_from_even_by_absolute_energies = [x==0 ? 0 : log(x) for x in number_swap_minima_configurations_from_even_by_absolute_energies]
    slice_minima_configuration_entropy_from_even_by_absolute_energies = [x==0 ? 0 : log(x) for x in number_slice_minima_configurations_from_even_by_absolute_energies]
    plot!(minima_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), swap_minima_configuration_entropy_from_even_by_absolute_energies, label="Swap Minima Configurations From Even", color=:green, line=:dash, lw=1)
    plot!(minima_configuration_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), slice_minima_configuration_entropy_from_even_by_absolute_energies, label="Slice Minima Configurations From Even", color=:orange, line=:dash, lw=1)

    savefig(minima_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations_even_too.png")
    savefig(minima_configuration_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations_even_too.svg")


    # Number of Minima Configurations From Even + Total Configurations Graph
    plot!(minima_configuration_from_even_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), swap_minima_configuration_entropy_from_even_by_absolute_energies, label="Swap Minima Configurations, ln(N⁺(E))", color=:blue, lw=1)
    plot!(minima_configuration_from_even_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), slice_minima_configuration_entropy_from_even_by_absolute_energies, label="Slice Minima Configurations, ln(N⁺(E))", color=:red, lw=1)

    if E_star != 0.0
        vline!(minima_configuration_from_even_graph, [E_star], linestyle=:dot, color=:green, label="")
        annotate!(minima_configuration_from_even_graph, [(E_star+0.02, ylims(minima_configuration_from_even_graph)[1]+50, Plots.text(L"E^*", 8, :black))])
    end

    # Print swap minima configuration entropy from even next to absolute energy
    for i in eachindex(slice_minima_configuration_entropy_from_even_by_absolute_energies)
        println("Swap Minima Configurations at $(absolute_energies_by_temperature[i]) = $(1-absolute_energies_by_temperature[i]/solved_configuration_energy(cube)): ", swap_minima_configuration_entropy_from_even_by_absolute_energies[i])
    end

    savefig(minima_configuration_from_even_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations_from_even.png")
    savefig(minima_configuration_from_even_graph, "results/relaxed_anneal_results/$(simulation_name)_minima_configurations_from_even.svg")

    # Downwards Connections Proportions Graph
    downwards_connections_proportion_graph = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), downwards_swap_connections_proportion_from_below_by_absolute_energies, xlabel="Energy, E", ylabel="Proportion of Downward Connections", title="Proportion of Downward Connections"; label="Swap Connections", color=:blue, legend=:topleft)
    plot!(downwards_connections_proportion_graph, absolute_energies_by_temperature./-solved_configuration_energy(cube), downwards_slice_connections_proportion_from_below_by_absolute_energies, label="Slice Connections", color=:red, legend=:topleft)

    savefig(downwards_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_downwards_connections_proportion.png")
    savefig(downwards_connections_proportion_graph, "results/relaxed_anneal_results/$(simulation_name)_downwards_connections_proportion.svg")

    # Downwards Connections Proportions using Even Connections Graph
    downwards_connections_proportion_graph_from_even = plot(absolute_energies_by_temperature./-solved_configuration_energy(cube), downwards_swap_connections_proportion_from_even_by_absolute_energies, xlabel="Energy, E", ylabel="Proportion of Downward Connections", title="Proportion of Downward Connections"; label="Swap Connections", color=:green, lw=1)
    plot!(downwards_connections_proportion_graph_from_even, absolute_energies_by_temperature./-solved_configuration_energy(cube), downwards_slice_connections_proportion_from_even_by_absolute_energies, label="Slice Connections", color=:orange, lw=1)

    # Plot in another colour just the places where the slice downwards connection proportion is zero
    # Get indices where slice connections are zero
    slice_zero_indices = findall(==(0.0), downwards_slice_connections_proportion_from_even_by_absolute_energies)
    plot!(downwards_connections_proportion_graph_from_even, absolute_energies_by_temperature[slice_zero_indices]./-solved_configuration_energy(cube), downwards_slice_connections_proportion_from_even_by_absolute_energies[slice_zero_indices], label="Proportion = Zero", color=:red, seriestype=:scatter, lw=1)

    savefig(downwards_connections_proportion_graph_from_even, "results/relaxed_anneal_results/$(simulation_name)_downwards_connections_proportion_even.png")
    savefig(downwards_connections_proportion_graph_from_even, "results/relaxed_anneal_results/$(simulation_name)_downwards_connections_proportion_even.svg")

    # TODO FOR LOOP DIFFERENT ENTROPIES

end

function reconstruct_histogram(histogram_data_name::String)

    # Step 1: Read the CSV file into a DataFrame
    histogram_data = CSV.read(histogram_data_name, DataFrame)

    # Step 2: Extract the bin edges and weights from the DataFrame
    x_edge_starts = histogram_data.E0
    y_edge_starts = histogram_data.En
    weights_flat = histogram_data.normalised_by_E0_weights

    x_edges_unique = unique(x_edge_starts)
    y_edges_unique = unique(y_edge_starts)

    num_bins_x = length(x_edges_unique)
    num_bins_y = length(y_edges_unique)

    # Step 4: Reshape the weights into a 2D array
    weights_matrix = reshape(weights_flat, num_bins_x, num_bins_y)

    # Step 5: Construct the Histogram object with the inferred edges and reshaped weights
    # We assume the last edge for each dimension is one step beyond the last unique edge found
    # This assumes the bin width is consistent and equal to the difference between consecutive edges
    x_edges_complete = [x_edges_unique; x_edges_unique[end] + (x_edges_unique[2] - x_edges_unique[1])]
    y_edges_complete = [y_edges_unique; y_edges_unique[end] + (y_edges_unique[2] - y_edges_unique[1])]
    reconstructed_edges = (x_edges_complete, y_edges_complete)

    # Create the Histogram object
    hist = Histogram(reconstructed_edges, weights_matrix)

    return hist

end

function window_average(array::Array{Float64, 1}, window_size::Int)
    # Ensure the window size is odd for symmetric averaging
    if window_size % 2 == 0
        error("Window size should be odd.")
    end

    len = length(array)
    half_window = window_size ÷ 2
    result = copy(array)

    for i in (half_window + 1):(len - half_window)
        result[i] = mean(array[(i - half_window):(i + half_window)])
    end

    return result
end


