using DelimitedFiles
using LaTeXStrings, Plots; 

using StatsBase
using CSV
# using DataFrames

using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

using CubicSplines

include("../core/rubiks_cube.jl")



# Main function
function boltzmann_shifted_energy_connectivity_histogram_figure(connectivity::String="Slice-Rotation"; neighbour_order_to_measure_to::Int64=1, temperature_shift::Float64=0.0)

    # E_star = -0.39015151515151514
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    # E_on = -0.23134348484848488
    E_on = -0.24010378787878783


    ### --- READ IN THE DATA ---
    if connectivity=="Slice-Rotation"
        simulation_name = "combined_L=11_inherent_disorder_slice"
        filename = "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_slice_E0_E1_energy_connections"
    elseif connectivity=="Swap-Move"
        simulation_name = "combined_L=11_inherent_disorder_swap"
        filename = "results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/",simulation_name*"_E0_E1_energy_connections"
    end

    header_line = readlines(joinpath(filename))[1]
    energy_connections_data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=2)


    # # Remove any rows with energy values larger than e10
    # energy_connections_data_matrix = energy_connections_data_matrix[energy_connections_data_matrix[:,1] .< 1e10, :]
    # # or with energy values >=0.0
    # energy_connections_data_matrix = energy_connections_data_matrix[energy_connections_data_matrix[:,1] .< abs(solved_configuration_energy(L)), :]
    # # or with energy values <= e-100
    # energy_connections_data_matrix = energy_connections_data_matrix[energy_connections_data_matrix[:,1] .> -1e-10, :]
    # # or with density values <= e-100
    # energy_connections_data_matrix = energy_connections_data_matrix[energy_connections_data_matrix[:,2] .> -1e-10, :]



    ### --- SET UP DEFAULT PARAMETERS ---

    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    # Remove NaNs and overflow 
    energy_connections_data_matrix = remove_bad_rows(energy_connections_data_matrix, L)
    Plots.default(dpi = 600)


    ### --- GET AVERAGE ENERGY AGAINST TEMPERATURE DATA ---
    filename = "results/final_paper_results",simulation_name*"_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"

    temperatures = [10.0, 9.54992586021436, 9.120108393559098, 8.709635899560805,8.317637711026709, 7.943282347242815, 7.5857757502918375, 7.244359600749901,6.918309709189366, 6.60693448007596, 6.3095734448019325, 6.025595860743578,5.7543993733715695, 5.495408738576245, 5.248074602497725, 5.011872336272722,4.786300923226383, 4.57088189614875, 4.36515832240166, 4.168693834703355,3.981071705534972, 3.801893963205612, 3.6307805477010135, 3.467368504525316,3.311311214825911, 3.1622776601683795, 3.019951720402016, 2.8840315031266055,2.7542287033381663, 2.6302679918953817, 2.51188643150958, 2.3988329190194904,2.2908676527677727, 2.1877616239495525, 2.089296130854039, 1.9952623149688797,1.9054607179632475, 1.8197008586099832, 1.7378008287493754, 1.6595869074375604,1.5848931924611134, 1.5135612484362082, 1.4454397707459277, 1.380384264602885,1.318256738556407, 1.2589254117941673, 1.202264434617413, 1.148153621496883,1.096478196143185, 1.0471285480508996, 1.0, 0.9549925860214359,0.9120108393559097, 0.8709635899560806, 0.831763771102671, 0.7943282347242814,0.7585775750291835, 0.7244359600749903, 0.6918309709189365, 0.6606934480075961,0.6309573444801932, 0.6025595860743578, 0.5754399373371569, 0.5495408738576245,0.5248074602497725, 0.5011872336272722, 0.4786300923226383, 0.45708818961487496,0.4365158322401659, 0.4168693834703355, 0.3981071705534973, 0.38018939632056126,0.36307805477010135, 0.34673685045253166, 0.33113112148259116, 0.3162277660168379,0.3019951720402016, 0.28840315031266056, 0.2754228703338166, 0.26302679918953814,0.25118864315095796, 0.23988329190194896, 0.22908676527677738, 0.2187761623949553,0.20892961308540398, 0.19952623149688797, 0.19054607179632474, 0.18197008586099836,0.17378008287493754, 0.16595869074375605, 0.15848931924611134, 0.1513561248436208,0.1445439770745927, 0.13803842646028847, 0.13182567385564076, 0.12589254117941676,0.12022644346174131, 0.11481536214968828, 0.1096478196143185, 0.10471285480508996]
    average_energy_densities_by_temperature = [-0.1809469696969697, -0.18082575757575758, -0.18192424242424243, -0.17976515151515152,-0.1805378787878788, -0.18525, -0.18692424242424244, -0.18406818181818183,-0.18787121212121213, -0.1863030303030303, -0.18712121212121213, -0.19046969696969698,-0.19040151515151515, -0.19242424242424241, -0.19165151515151516, -0.1923030303030303,-0.19496212121212123, -0.19711363636363635, -0.20083333333333336, -0.20047727272727273,-0.20236363636363636, -0.20400757575757578, -0.20588636363636362, -0.20875000000000002,-0.21309848484848487, -0.2115757575757576, -0.2155, -0.21622727272727274,-0.21995454545454549, -0.2198106060606061, -0.2242121212121212, -0.22921969696969696,-0.2314090909090909, -0.2324848484848485, -0.23884090909090908, -0.24062121212121212,-0.24951515151515152, -0.2516515151515152, -0.258530303030303, -0.2602878787878788,-0.26795454545454545, -0.274, -0.28042424242424246, -0.2888939393939394,-0.2936590909090909, -0.30199242424242423, -0.31156060606060604, -0.3235227272727273,-0.33324242424242423, -0.34381060606060604, -0.3576212121212121, -0.36961363636363637,-0.3893030303030303, -0.4058257575757576, -0.43169696969696975, -0.4528560606060606,-0.48315909090909087, -0.5125, -0.5389015151515152, -0.5718484848484848,-0.5935530303030303, -0.612, -0.6338863636363636, -0.6457954545454546,-0.6553484848484848, -0.6639848484848485, -0.6690681818181818, -0.6777272727272727,-0.686060606060606, -0.6898863636363636, -0.6929015151515151, -0.6951590909090909,-0.6950303030303031, -0.6968106060606061, -0.6981060606060606, -0.6991893939393939,-0.7018560606060606, -0.7027348484848485, -0.7056742424242425, -0.7068560606060607,-0.7070454545454546, -0.707780303030303, -0.7078636363636364, -0.7083181818181818,-0.7086287878787879, -0.7083181818181818, -0.708939393939394, -0.7088181818181818,-0.7088030303030303, -0.7090151515151515, -0.7089772727272727, -0.7090681818181819,-0.7090681818181819, -0.7090757575757576, -0.7090681818181819, -0.7090454545454545,-0.7090757575757576, -0.7090909090909091, -0.7090909090909091, -0.7090909090909091]
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



    # Plot mean j-i against i
    # Get smoothened data averaged over n_average energy units each
    smoothing_window = 60
    mean_j_minus_i_smoothed = [mean(mean_j_minus_i[i:i+smoothing_window]) for i in 1:smoothing_window:length(mean_j_minus_i)-smoothing_window]
    # Dont plot 0.0s

    # INSET PLOT
    yticks = connectivity == "Slice-Rotation" ? [0,0.005,0.01] : [-0.001,0.,0.001]
    # graph = scatter(-1 .+ (bin_edges_x[1:end-1][mean_j_minus_i .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i[mean_j_minus_i .!= 0.0]./abs(solved_configuration_energy(cube)), label=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", title="", legend=false, xtickfontsize=30, ytickfontsize=30, xguidefontsize=30, yguidefontsize=30, margin=5mm, xticks=[-0.3,-0.5,-0.7], yticks=yticks)

    # SM PLOT
    graph = scatter(-1 .+ (bin_edges_x[1:end-1][mean_j_minus_i .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i[mean_j_minus_i .!= 0.0]./abs(solved_configuration_energy(cube)), label=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", title="", legend=false, margin=5mm, yticks=yticks)

    # Do smoothed line
    if connectivity=="Swap-Move"
        plot!(graph, -1 .+ (bin_edges_x[1:end-1][1:smoothing_window:end-smoothing_window][mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i_smoothed[mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube)), line=:solid, color=:orange, lw=3, xlims=(-0.71,-0.25), xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", ylims=(minimum(yticks), maximum(yticks)))
        # plot!(graph, -1 .+ (bin_edges_x[1:end-1][1:smoothing_window:end-smoothing_window][mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i_smoothed[mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube)), line=:solid, color=:orange, lw=3, xlims=(-0.71,-0.25), xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle")

    elseif connectivity=="Slice-Rotation"
        plot!(graph, -1 .+ (bin_edges_x[1:end-1][1:smoothing_window:end-smoothing_window][mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i_smoothed[mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube)), line=:solid, color=:orange, lw=3, xlims=(-0.71,-0.25), xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle", ylims=(minimum(yticks)-0.003, maximum(yticks)+0.001))
        # plot!(graph, -1 .+ (bin_edges_x[1:end-1][1:smoothing_window:end-smoothing_window][mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube))), mean_j_minus_i_smoothed[mean_j_minus_i_smoothed .!= 0.0]./abs(solved_configuration_energy(cube)), line=:solid, color=:orange, lw=3, xlims=(-0.71,-0.25), xlabel=L"\epsilon^{(0)}", ylabel=L"\langle \langle \epsilon^{(1)} - \epsilon^{(0)} \rangle \rangle")

    end

    # Plot 0.0 line
    hline!(graph, [0.0], line=:dash, color=:red, lw=2, label="")

    # Plot E* line
    if connectivity == "Slice-Rotation"
        vline!(graph, [E_star], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star+0.015, ylims(graph)[1]+0.0005, Plots.text(L"\bar\epsilon^*", 10, :black))])
    end

    # Plot E_on line
    if connectivity == "Slice-Rotation"
        vline!(graph, [E_on], line=:dash, color=:red, lw=2, label="")
    end

        

    # Add title as annotated text in top right corner
    # annotate!(graph, [(xlims(graph)[2]-200, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])

    # Anotate temperature shift if not 0
    if temperature_shift != 0.0
        annotate!(graph, [(xlims(graph)[1]+250, ylims(graph)[2]-0.25, Plots.text("Temperature Shift: $temperature_shift", 10, :black))])
    end

    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_boltzmann_shifted.png")
    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_boltzmann_shifted.pdf")
    display(graph)
end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end


boltzmann_shifted_energy_connectivity_histogram_figure("Slice-Rotation")
boltzmann_shifted_energy_connectivity_histogram_figure("Swap-Move")