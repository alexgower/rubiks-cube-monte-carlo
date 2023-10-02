using DelimitedFiles
using Plots
using StatsBase
using LaTeXStrings

function phase_free_energy_plotter(graph_save_name::String, data_simulation_name_to_use::String, ordered_phase_cutoff_energy::Float64; log_scale::Bool=false, plot_disordered_phase::Bool=true, plot_ordered_phase::Bool=true)

    # Read data from file 
    data_matrix = readdlm(joinpath("results/tools/phase_free_energy_plotter/",data_simulation_name_to_use), ',', Any, '\n', skipstart=2)
    temperatures = convert(Vector{Float64},data_matrix[:,1])


    # -------------------------------------

    # New Style
    # Remove '[' and ']' characters at start and end of energy sample array
    start_elements = convert(Vector{String},data_matrix[:,2])
    start_elements = [element[3:end] for element in start_elements]
    start_elements = parse.(Float64,start_elements)
    
    end_elements = convert(Vector{String},data_matrix[:,end])
    end_elements = [element[1:end-2] for element in end_elements]
    end_elements = parse.(Float64,end_elements)

    energy_samples_by_temperature = zeros(length(temperatures), length(data_matrix[1,:])-1) # -1 to remove temperature column

    energy_samples_by_temperature[:,1] .= start_elements
    energy_samples_by_temperature[:,2:end-1] .= convert(Matrix{Float64},data_matrix[:,3:end-1])
    energy_samples_by_temperature[:,end] .= end_elements

    # -------------------------------------



    # Calculate Phase Free Energies

    ordered_phase_free_energy_by_temperature = zeros(length(temperatures))
    disordered_phase_free_energy_by_temperature = zeros(length(temperatures))

    println(temperatures)

    for (index,temperature) in pairs(temperatures)

        # Create a histogram
        histogram = countmap(energy_samples_by_temperature[index,:])

        energy_bins = collect(keys(histogram))
        println("Energies: ", energy_bins)

        probability_densities = collect(values(histogram))
        println("Probability densities: ", probability_densities)

        # Normalize the values
        probability_densities = probability_densities ./ sum(probability_densities)
        println("Normalised Probability densities: ", probability_densities)


        # Calculate ordered phase free energy 
        # Using: \boxed{F_{ord} = \sum_{E_i=0}^{E^*} n(E_i)E_i}
        # where E^* = ordered phase energy cutoff and n(E_i) is the probability density of energy E_i

        ordered_phase_cutoff_energy_index = findfirst(energy_bins .> ordered_phase_cutoff_energy)
        if ordered_phase_cutoff_energy_index == nothing
            ordered_phase_cutoff_energy_index = length(energy_bins)+1
        end

        println("Ordered Phase Cutoff Energy Index: ", ordered_phase_cutoff_energy_index)

        # My definition
        ordered_phase_free_energy_by_temperature[index] = sum(probability_densities[1:ordered_phase_cutoff_energy_index-1] .* energy_bins[1:ordered_phase_cutoff_energy_index-1])
        println("My Ordered Phase Free Energy: ", ordered_phase_free_energy_by_temperature[index])
        # Claudio Definition
        # ordered_phase_free_energy_by_temperature[index] = sum(probability_densities[1:ordered_phase_cutoff_energy_index-1] .* (energy_bins[1:ordered_phase_cutoff_energy_index-1] .- temperature .* log.(probability_densities[1:ordered_phase_cutoff_energy_index-1])))
        # println("Ordered Phase Free Energy: ", ordered_phase_free_energy_by_temperature[index])


        # Calculate disordered phase free energy
        # Using: \boxed{F_{dis} = \sum_{E_i=E^*}^{\infty} n(E_i)E_i}
        # where E^* = ordered phase energy cutoff and n(E_i) is the probability density of energy E_i

        # My definition
        disordered_phase_free_energy_by_temperature[index] = sum(probability_densities[ordered_phase_cutoff_energy_index:end] .* energy_bins[ordered_phase_cutoff_energy_index:end])
        println("My Disordered Phase Free Energy: ", disordered_phase_free_energy_by_temperature[index])
        # Claudio Definition
        # disordered_phase_free_energy_by_temperature[index] = sum(probability_densities[ordered_phase_cutoff_energy_index:end] .* (energy_bins[ordered_phase_cutoff_energy_index:end] .- temperature .* log.(probability_densities[ordered_phase_cutoff_energy_index:end])))
        # println("Disordered Phase Free Energy: ", disordered_phase_free_energy_by_temperature[index])


    end

    println("Ordered Phase Free Energies: ", ordered_phase_free_energy_by_temperature)
    println("Disordered Phase Free Energies: ", disordered_phase_free_energy_by_temperature)
    println("Temperatures: ", temperatures)
    println(temperatures[26])

    # Create graph axis
    graph = log_scale ? plot(xlabel="Temperature", ylabel="Free Energy", title="Ordered and Disordered Phase Free Energies", xscale=:log10, legend=:topright) : plot(xlabel="Temperature", ylabel="Free Energy", title="Ordered and Disordered Phase Free Energies", legend=:bottomleft)
   

    # Strip out 0 Free Energies (where phase doesn't exist in histogram) and corresponding temperatures
    # ordered_phase_temperatures = temperatures[ordered_phase_free_energy_by_temperature .!= 0.0]
    # ordered_phase_free_energy_by_temperature = ordered_phase_free_energy_by_temperature[ordered_phase_free_energy_by_temperature .!= 0.0]
    # disordered_phase_temperatures = temperatures[disordered_phase_free_energy_by_temperature .!= 0.0]
    # disordered_phase_free_energy_by_temperature = disordered_phase_free_energy_by_temperature[disordered_phase_free_energy_by_temperature .!= 0.0]

    # Plot Free Energies
    if plot_disordered_phase
        plot!(graph, temperatures, disordered_phase_free_energy_by_temperature, label=L"F_{dis}", color=:red)
    end
    if plot_ordered_phase
        plot!(graph, temperatures, ordered_phase_free_energy_by_temperature, label=L"F_{ord}", color=:green)
    end

    # Move legend
    graph = plot!(legend=:bottomright, legendfontsize=12)



    ### SAVE GRAPH
    savefig(graph, "results/tools/phase_free_energy_plotter/$(graph_save_name).png")


    display(graph)
end
