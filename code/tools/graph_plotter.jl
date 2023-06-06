using DelimitedFiles
using Plots

function plot_graph(simulation_name::String, swap_move_probabilities::Vector{Float64}, L::Int64, T_swap::Float64; normalization::String="solved")

    temperature_vector::Vector{Float64} = []
    normalised_E_average_by_temperature = []
    infinite_temperature_normalised_E_average_by_temperature = []

    for (index,swap_move_probability) in pairs(swap_move_probabilities)
        simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)

        # Read data from file 
        data_matrix = readdlm(joinpath("results",simulation_name_to_use), ',', Float64, '\n', skipstart=2)

        temperature_vector = copy(data_matrix[:,1])
        E_average_by_temperature = data_matrix[:,2]
        this_normalised_E_average_by_temperature = data_matrix[:,3]

        # -------------------------------------


        println("<E>(T)")
        println(E_average_by_temperature)

        push!(normalised_E_average_by_temperature, this_normalised_E_average_by_temperature)
        push!(infinite_temperature_normalised_E_average_by_temperature,this_normalised_E_average_by_temperature .* 6)

        println("<-E/E_0>(T)")
        println(normalised_E_average_by_temperature[index])

        println("<-E/E_infty>(T)")
        println(infinite_temperature_normalised_E_average_by_temperature[index])

        # ------------------------------
    end

    

    # Create Plot --------------
    if normalization == "solved"

        graph = plot(temperature_vector, normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Solved Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
        hline!(graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
        hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")

    elseif normalization == "infinite_temperature"

        graph = plot(temperature_vector, infinite_temperature_normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Infinite Temperature Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
        hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")
        hline!(graph, [-6.0], linestyle=:dash, color=:black, label="")

    end

    vline!([T_swap], linestyle=:dash, label="Swap Moves Enabled")


    # Add other data to graph for comparison if exists ----------
    
    if isfile(joinpath("results","other_data.csv"))
        data_matrix = readdlm(joinpath("results","other_data.csv"), ',', Float64, '\n', skipstart=0)

        other_temperature_vector = copy(data_matrix[:,1])
        other_data = data_matrix[:,2]
        
        plot!(graph, other_temperature_vector, other_data, label="Ollie Results", seriestype=:scatter, color="blue", ms=2, ma=0.5)
    end

    # Save Plot ----------
    savefig("results/new_$simulation_name.png")
end