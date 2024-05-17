using DelimitedFiles
using LaTeXStrings, Plots; 

using StatsBase
using CSV
# using DataFrames

using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

using Images
using CubicSplines

include("../core/rubiks_cube.jl")



# Main function
function boltzmann_shifted_energy_connectivity_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation"; neighbour_order_to_measure_to::Int64=1, temperature_shift::Float64=0.0)

    E_star = -0.39015151515151514
    

    ### --- READ IN THE DATA ---
    filename = "results/final_paper_results",simulation_name*"_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"

    header_line = readlines(joinpath(filename))[1]
    energy_connections_data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=2)


    ### --- SET UP DEFAULT PARAMETERS ---

    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    # Remove NaNs and overflow 
    energy_connections_data_matrix = remove_bad_rows(energy_connections_data_matrix, L)
    Plots.default(dpi = 300)


    ### --- GET AVERAGE ENERGY AGAINST TEMPERATURE DATA ---
    # TODO GET NICE DATA FOR THIS DISORDER INSTANCE NOT INHERENT DISORDER AVERAGE
    # TODO check spline nicely interpolates curve
    filename = "results/final_paper_results",simulation_name*"_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"

    temperatures = [10.0, 9.440608762859235, 8.912509381337456, 8.413951416451951, 7.943282347242815, 7.498942093324558, 7.079457843841379, 6.6834391756861455, 6.3095734448019325, 5.956621435290105, 5.623413251903491, 5.3088444423098835, 5.011872336272722, 4.7315125896148045, 4.466835921509632, 4.216965034285822, 3.981071705534972, 3.758374042884442, 3.548133892335754, 3.3496543915782766, 3.1622776601683795, 2.9853826189179595, 2.8183829312644537, 2.6607250597988097, 2.51188643150958, 2.371373705661655, 2.2387211385683394, 2.1134890398366464, 1.9952623149688797, 1.8836490894898006, 1.7782794100389228, 1.6788040181225603, 1.5848931924611134, 1.4962356560944334, 1.4125375446227544, 1.333521432163324, 1.2589254117941673, 1.1885022274370183, 1.1220184543019636, 1.0592537251772889, 1.0, 0.9440608762859236, 0.9, 0.8912509381337455, 0.8842105263157894, 0.868421052631579, 0.8526315789473685, 0.8413951416451952, 0.8368421052631578, 0.8210526315789474, 0.8052631578947369, 0.7943282347242814, 0.7894736842105263, 0.7736842105263158, 0.7578947368421052, 0.7498942093324558, 0.7421052631578948, 0.7263157894736842, 0.7105263157894737, 0.707945784384138, 0.6947368421052631, 0.6789473684210526, 0.6683439175686146, 0.6631578947368421, 0.6473684210526315, 0.631578947368421, 0.6309573444801932, 0.6157894736842106, 0.6, 0.5956621435290104, 0.5623413251903491, 0.5308844442309885, 0.5011872336272722, 0.4731512589614806, 0.4466835921509631, 0.4216965034285822, 0.3981071705534973, 0.37583740428844414, 0.3548133892335755, 0.3349654391578276, 0.3162277660168379, 0.298538261891796, 0.28183829312644537, 0.26607250597988097, 0.25118864315095796, 0.23713737056616555, 0.22387211385683398, 0.21134890398366468, 0.19952623149688797, 0.18836490894898, 0.1778279410038923, 0.16788040181225608, 0.15848931924611134, 0.14962356560944334, 0.1333521432163324, 0.12589254117941676, 0.11220184543019636, 0.10592537251772886]
    average_energy_densities_by_temperature = [-0.17839177489177488, -0.17883874458874466, -0.17976443001443004, -0.18036219336219336, -0.18141558441558447, -0.1828080808080808, -0.183728354978355, -0.18547474747474751, -0.18624711399711402, -0.18801767676767678, -0.1880606060606061, -0.19049314574314574, -0.19169480519480522, -0.193234126984127, -0.1955667388167388, -0.1977027417027417, -0.1991569264069264, -0.20138564213564208, -0.203487012987013, -0.2056634199134199, -0.20901623376623377, -0.21251839826839827, -0.21551190476190474, -0.21867352092352088, -0.22276984126984123, -0.22641847041847044, -0.23113780663780664, -0.2357734487734488, -0.24022077922077925, -0.24576803751803755, -0.2515335497835498, -0.25796536796536795, -0.2652431457431457, -0.2731818181818182, -0.2814329004329005, -0.2906125541125541, -0.3006432178932179, -0.31308658008658014, -0.32555158730158734, -0.3393001443001442, -0.3556926406926407, -0.37400360750360756, -0.3929195526695527, -0.3963459595959596, -0.40076623376623377, -0.4074029581529582, -0.41649891774891784, -0.4224444444444444, -0.42484163059163066, -0.43698124098124097, -0.44870779220779233, -0.45488455988456, -0.4589689754689754, -0.4712049062049062, -0.485448051948052, -0.4897723665223664, -0.4982564935064936, -0.5104635642135641, -0.5246349206349207, -0.5262352092352092, -0.5393780663780663, -0.5484971139971141, -0.5559325396825398, -0.5631518759018758, -0.5771948051948053, -0.5896998556998556, -0.5935191197691199, -0.6001255411255413, -0.6055147907647906, -0.612765873015873, -0.62870202020202, -0.6426475468975469, -0.6539480519480521, -0.6623708513708514, -0.6690212842712844, -0.67572113997114, -0.6799285714285713, -0.6834220779220779, -0.6857106782106781, -0.6885407647907649, -0.6905515873015873, -0.6923257575757576, -0.6933246753246752, -0.694488455988456, -0.6954289321789321, -0.695962481962482, -0.696670634920635, -0.6974556277056279, -0.6976601731601731, -0.6978968253968254, -0.6981829004329003, -0.6984693362193363, -0.6986879509379509, -0.6988181818181818, -0.6988398268398269, -0.6988896103896105, -0.6989491341991342, -0.6989949494949496]
    average_energies_by_temperature = average_energy_densities_by_temperature .* -solved_configuration_energy(cube)

    # highlight which energies in average_energies_by_temperature are not in descending order
    for i in 1:size(average_energies_by_temperature, 1)
        if i > 1 && average_energies_by_temperature[i] > average_energies_by_temperature[i-1]
            println("Energy at index $i is not in descending order: $(average_energies_by_temperature[i])")
            println("Corresponding temperature: $(temperatures[i])")
            println("Corresponding energy density: $(average_energy_densities_by_temperature[i])")
        end
    end

    # Create a cubic spline interpolation object
    spline = CubicSpline(sort(average_energies_by_temperature), sort(temperatures), extrapl=[1,], extrapr=[1,])

    # Define the function that uses this spline object to find the equilibrium temperature for a given energy
    function get_equilibrium_temperature_from_energy(energy::Float64)::Float64
        return spline(energy)
    end

    function get_equilibrium_temperature_from_energy_density(energy_density::Float64)::Float64
        return get_equilibrium_temperature_from_energy(energy_density * -solved_configuration_energy(cube))
    end

    function get_equilibrium_temperature_from_absolute_energy(absolute_energy)::Float64
        return get_equilibrium_temperature_from_energy(absolute_energy - abs(solved_configuration_energy(cube)))
    end

    # Test the get_equilibrium_temperature_from_absolute_energy function for 0 to 1320 a few values
    for i in 0:100:1320
        println("E: $i, T: $(get_equilibrium_temperature_from_absolute_energy(i))")
    end

    println("epsilon: -0.3, T: $(get_equilibrium_temperature_from_absolute_energy_density(-0.3))")
    


    ### --- HISTOGRAM GENERAL CALCULATIONS ---

    # Determine the bin edges based on the data
    bin_edges_x = 0:1:-solved_configuration_energy(cube)
    bin_edges_y = 0:1:-solved_configuration_energy(cube)
    edges = (bin_edges_x, bin_edges_y)

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (energy_connections_data_matrix[:,1], energy_connections_data_matrix[:,2]), edges)

    boltzmann_shifted_weights = zeros(size(hist_2d.weights))

    # Apply the Boltzmann factor to the weights
    for i in 1:size(hist_2d.weights, 1)
        for j in 1:size(hist_2d.weights, 2)
            # Note in this basis we have energies from 0 to 1320
            # But since energy discretisation is 1, we can just use the indices i and j
            # Where i = E^{(0)} and j = E^{(1)}
            boltzmann_factor = min(1, exp(-(j-i)/get_equilibrium_temperature_from_absolute_energy(i)+temperature_shift))
            boltzmann_shifted_weights[i, j] = hist_2d.weights[i, j] * boltzmann_factor
        end
    end

    # Normalize the weights by each E^{(0)}
    for i in 1:size(boltzmann_shifted_weights, 1)
        slice_sum = sum(boltzmann_shifted_weights[i, :])
        if slice_sum != 0  # Avoid division by zero
            boltzmann_shifted_weights[i, :] ./= slice_sum
        end
    end



    # Calculate mean j for each i
    # TODO check end-1 calculation
    mean_j = zeros(size(boltzmann_shifted_weights, 1))
    for i in 1:size(boltzmann_shifted_weights, 1)
        mean_j[i] = sum(boltzmann_shifted_weights[i, :] .* bin_edges_y[1:end-1])
    end

    # Calculate mean j-i for each if
    mean_j_minus_i = zeros(size(boltzmann_shifted_weights, 1))
    for i in 1:size(boltzmann_shifted_weights, 1)
        mean_j_minus_i[i] = sum(boltzmann_shifted_weights[i, :] .* (bin_edges_y[1:end-1] .- bin_edges_x[i]))
    end


    # Print each i and each mean j-i
    for i in 1:size(boltzmann_shifted_weights, 1)
        println("i: $(bin_edges_x[i]), mean j-i: $(mean_j_minus_i[i])")
    end


    # Print E_star as absolute energy
    println("E_star: $(E_star)")
    println("E_star as absolute energy: $(E_star*abs(solved_configuration_energy(cube)) + abs(solved_configuration_energy(cube))))")

    
    # Plot mean j-i against i
    # Get smoothened data averaged over n_average energy units each
    mean_j_minus_i_smoothed = [mean(mean_j_minus_i[i:i+20]) for i in 1:20:length(mean_j_minus_i)-20]
    # Dont plot 0.0s
    graph = scatter(bin_edges_x[1:end-1][mean_j_minus_i .!= 0.0], mean_j_minus_i[mean_j_minus_i .!= 0.0], label=L"\langle E^{(1)} - E^{(0)} \rangle", xlabel=L"E^{(0)}", ylabel=L"\langle E^{(1)} - E^{(0)} \rangle", title="", legend=:topright)

    # Do smoothed line
    plot!(graph, bin_edges_x[1:end-1][1:20:end-20][mean_j_minus_i_smoothed .!= 0.0], mean_j_minus_i_smoothed[mean_j_minus_i_smoothed .!= 0.0], line=:dash, color=:orange, lw=2, label="Smoothed")

    # Plot 0.0 line
    hline!(graph, [0.0], line=:dash, color=:red, lw=2, label="")

    # Plot E_star line for it's absolute energy
    E_star_absolute_energy = E_star*abs(solved_configuration_energy(cube)) + abs(solved_configuration_energy(cube))
    vline!(graph, [E_star_absolute_energy], line=:dash, color=:green, lw=2, label=L"E^{*}")

    # Add title as annotated text in top right corner
    annotate!(graph, [(xlims(graph)[2]-200, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])

    # Anotate temperature shift if not 0
    if temperature_shift != 0.0
        annotate!(graph, [(xlims(graph)[1]+250, ylims(graph)[2]-0.25, Plots.text("Temperature Shift: $temperature_shift", 10, :black))])
    end

    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_mean_j_minus_i.png")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_mean_j_minus_i.svg")
    display(graph)

    println("Size of bin_edges_x: ", size(bin_edges_x))
    println("Size of bin_edges_y: ", size(bin_edges_y))
    println("Size of boltzmann_shifted_weights: ", size(boltzmann_shifted_weights))

    indices_where_mean_j_minus_i_is_positive = findall(mean_j_minus_i .> 0.0)

    relative_boltzmann_shifted_weights = zeros(1320, 2*1320)

    for i in 1:1320
        for j in 1:1320
            relative_boltzmann_shifted_weights[i, j-i+1320] = boltzmann_shifted_weights[i, j]
        end
    end

    # Check relative boltzmann shifted weights is still normalised per i
    for i in 1:1320
        slice_sum = sum(relative_boltzmann_shifted_weights[i, :])
        println("Slice sum: ", slice_sum)
    end

    println("Relative boltzmann shifted weights all zero?: ", all(relative_boltzmann_shifted_weights .== 0.0))

    ### --- PLOT RELATIVE BOLTZMANN SHIFTED WEIGHTS ON 2D HISTOGRAM ---
    graph = heatmap(
        [i for i in 1:1320], 
        [i for i in -1320:1320 if i != 0],
        relative_boltzmann_shifted_weights, 
        # color=:bluesreds, 
        show_empty_bins=false, 
        xlabel=L"E^{(0)}", 
        ylabel=L"E^{(1)}", 
        title="", 
        xlims = extrema(indices_where_mean_j_minus_i_is_positive),
        colorbar_title="",
        # xticks=(0:100:1200, 0:100:1200),
        # yticks=(0:100:1200, 0:100:1200),
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=5mm,          # Margin around the plot
        # background_color=:white,
        # grid=false
    )


    # graph = heatmap(
    #     bin_edges_x, 
    #     bin_edges_y,
    #     boltzmann_shifted_weights, 
    #     color=:bluesreds, 
    #     show_empty_bins=false, 
    #     xlabel=L"E^{(0)}", 
    #     ylabel=L"E^{(1)}", 
    #     title="", 
    #     xlims = extrema(indices_where_mean_j_minus_i_is_positive),
    #     ylims = extrema(indices_where_mean_j_minus_i_is_positive),
    #     colorbar_title="",
    #     # xticks=(0:100:1200, 0:100:1200),
    #     # yticks=(0:100:1200, 0:100:1200),
    #     xguidefontsize=12,   # X-axis label font size
    #     yguidefontsize=12,   # Y-axis label font size
    #     margin=5mm,          # Margin around the plot
    #     # background_color=:white,
    #     # grid=false
    # )


    ### -- Plots. jl 2D HISTOGRAM --

    # min_value = minimum([minimum(energy_connections_data_matrix[:,1]), minimum(energy_connections_data_matrix[:,2])])
    # max_value = maximum([maximum(energy_connections_data_matrix[:,1]), maximum(energy_connections_data_matrix[:,2])])


    # # # Initialize an array to hold the biases
    # # biases = zeros(size(energy_connections_data_matrix, 1))

    # # # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
    # # # Calculate the sum of counts in each E_0 slice
    # # for i in 1:size(hist_2d.weights, 1)
    # #     slice_sum = sum(hist_2d.weights[i, :])
    # #     if slice_sum != 0  # Avoid division by zero
    # #         # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_0 slice
    # #         mask = energy_connections_data_matrix[:,1] .== bin_edges_x[i]
    # #         biases[mask] .= 1.0 / slice_sum
    # #     end
    # # end

    # E0_values = energy_connections_data_matrix[:,1]./-solved_configuration_energy(cube)
    # E1_values = energy_connections_data_matrix[:,2]./-solved_configuration_energy(cube)

    # diagonal_ylabel = "Second-Neighbour Energy Density, "*L"\epsilon^{(2)}"
    # diagonal_xlabel = "Neighbour Energy Density, "*L"\epsilon^{(1)}"

    # graph = histogram2d(
    #     E0_values, 
    #     E1_values, 
    #     color=:bluesreds, 
    #     show_empty_bins=false,
    #     normalize=:pdf, 
    #     bins=(bin_edges_x./-solved_configuration_energy(cube), bin_edges_y./-solved_configuration_energy(cube)), 
    #     weights=boltzmann_shifted_weights, 
    #     ylabel=diagonal_ylabel, 
    #     xlabel=diagonal_xlabel, 
    #     # title="$connectivity Cube", 
    #     xlims=(minimum(E0_values), maximum(E0_values)), 
    #     ylims=(minimum(E1_values), maximum(E1_values)),
    #     colorbar_title="",
    #     # titlefontsize=10,   # Title font size
    #     xguidefontsize=12,   # X-axis label font size
    #     yguidefontsize=12,   # Y-axis label font size
    #     margin=5mm,          # Margin around the plot
    #     xticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
    #     yticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
    #     # background_color=:white,
    #     # grid=false
    # )
    # annotate!(1.115 * maximum(E0_values), 0.5 * (minimum(E1_values) + maximum(E1_values)), text("Sampled Frequency", 10, :center, :center, rotation=90))

    # annotation = L"\epsilon^{(1)} = \epsilon^{(0)}"
    # annotate!(graph, [(xlims(graph)[1]+0.1, ylims(graph)[1]+0.03, Plots.text(annotation, 12, :black))])

    # plot!(graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:orange, lw=2, label="")

    # # Add E_star vertical lines to graphs if we are slice rotation cube
    # if connectivity == "Slice-Rotation"
    #     E_star_plot = 1+E_star
    #     vline!(graph, [E_star_plot], line=:dash, color=:green, lw=2, label="")
    #     annotate!(graph, [(E_star_plot+0.02, ylims(graph)[1]+0.05, Plots.text(L"\epsilon^*", 12, :black))])
    # end

    # # Add title as annotated text in top right corner
    # if connectivity=="Slice-Rotation"
    #     annotate!(graph, [(xlims(graph)[2]-0.1, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])
    # else
    #     annotate!(graph, [(xlims(graph)[2]-0.1, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])
    # end

    # ### -- Save and display the graphs --
    # println("Saving diagonal graph...")
    # savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal_boltzmann_shifted.svg")
    # savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal_boltzmann_shifted.png")
    # display(graph)


end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end