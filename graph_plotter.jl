using DelimitedFiles
using Plots

function plot_graph(simulation_name::String, swap_move_probabilities::Vector{Float64}, L::Int64, T_swap::Float64; alternative_normalization=false)

    # Cover everything in try/except clause
    # try
        temperature_vector::Vector{Float64} = []
        normalised_E_average_by_temperature = []
        alternative_normalised_E_average_by_temperature = []

        for (index,swap_move_probability) in pairs(swap_move_probabilities)
            simulation_name_to_use = simulation_name * string(swap_move_probability)

            # Read data from file 
            data_matrix = readdlm(joinpath("results",simulation_name_to_use), ',', Float64, '\n', skipstart=2)

            println(data_matrix)

            temperature_vector = copy(data_matrix[:,1])
            E_average_by_temperature = data_matrix[:,2]
            this_normalised_E_average_by_temperature = data_matrix[:,3]
            E_squared_average_by_temperature = data_matrix[:,4]
            relaxation_iterations_by_temperature = data_matrix[:,5]
            accepted_candidates_by_temperature = data_matrix[:,6]
            final_configuration_correlation_function_by_temperature = data_matrix[:,7]


            # -------------------------------------


            println("<E>(T)")
            println(E_average_by_temperature)

            push!(normalised_E_average_by_temperature, this_normalised_E_average_by_temperature)
            push!(alternative_normalised_E_average_by_temperature,this_normalised_E_average_by_temperature .* 6)


            println("<-E/E_0>(T)")
            println(normalised_E_average_by_temperature[index])

            println("<-E/E_infty>(T)")
            println(alternative_normalised_E_average_by_temperature[index])

            # ------------------------------
        end

        

        # Display and Save Plot --------------
        if alternative_normalization

            plot(temperature_vector, alternative_normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Infinite Temperature Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))

        else
            plot(temperature_vector, normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Solved Configuration Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
        end
        vline!([T_swap], linestyle=:dash, label="Swap Moves Enabled")

        savefig("results/new_$simulation_name.png")
end