using DelimitedFiles
using Plots
using LaTeXStrings

function parallel_tempering_figure()
    no_pt_file = "results/final_paper_results/no_PT_swap_moves.csv"
    pt_file = "results/final_paper_results/PT_moves.csv"

    # Read the CSV files into arrays, skipping the header line # i.e. the Ollie Data
    no_pt_data = readdlm(no_pt_file, ',', Float64, skipstart=1)
    pt_data = readdlm(pt_file, ',', Float64, skipstart=1)

    # Extract temperature and energy columns
    no_pt_temp = no_pt_data[:, 1]
    no_pt_energy = no_pt_data[:, 2]
    pt_temp = pt_data[:, 1]
    pt_energy = pt_data[:, 2]

    # Rescale energies from Ollie units to our units by dividing by 6
    no_pt_energy = no_pt_energy / 6
    pt_energy = pt_energy / 6

    # Rescale temperatures from Ollie units to our units by dividing by 2
    no_pt_temp = no_pt_temp / 2
    pt_temp = pt_temp / 2

    # Order the data by temperature
    no_pt_temp, no_pt_energy = sort(no_pt_temp), no_pt_energy[sortperm(no_pt_temp)]
    pt_temp, pt_energy = sort(pt_temp), pt_energy[sortperm(pt_temp)]





    # Get our own L=9 data
    filenames_that_do_not_exist=[]
    model = "clean"
    L = 9
    swap_move_probability = 0.0
    N_T = 100
    temperatures = zeros(N_T)
    running_total_average_energy_densities_by_temperature = zeros(N_T)
    trials = 50
    actual_number_of_trials=0
    for trial in 1:trials
        filename = "results/final_paper_results/relaxed_anneal_results/"* model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
        try
            data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
            
            temperatures .= data_matrix[:,1]
            running_total_average_energy_densities_by_temperature .+= data_matrix[:,3]
            
            actual_number_of_trials += 1
        catch e
            push!(filenames_that_do_not_exist, filename)
        end
    
    end

    alex_slice_temperatures = temperatures
    alex_slice_average_energy_densities = running_total_average_energy_densities_by_temperature/actual_number_of_trials


    filenames_that_do_not_exist=[]
    model = "clean"
    L = 9
    swap_move_probability = 1.0
    N_T = 100
    temperatures = zeros(N_T)
    running_total_average_energy_densities_by_temperature = zeros(N_T)
    trials = 50
    actual_number_of_trials=0
    for trial in 1:trials
        filename = "results/final_paper_results/relaxed_anneal_results/"* model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"
        try
            data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
            
            temperatures .= data_matrix[:,1]
            running_total_average_energy_densities_by_temperature .+= data_matrix[:,3]
            
            actual_number_of_trials += 1
        catch e
            push!(filenames_that_do_not_exist, filename)
        end
    
    end

    alex_swap_temperatures = temperatures
    alex_swap_average_energy_densities = running_total_average_energy_densities_by_temperature/actual_number_of_trials




    println("Filenames that do not exist: ", filenames_that_do_not_exist)





    ### --- COLOURS ---
    Plots.default(dpi = 300)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)

    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Create the plot
    # Add PT data to the same plot
    plot(pt_temp, pt_energy, 
    label="Parallel Tempering", 
    # marker=:square,
    color=alex_red)

    minimum_parallel_tempering_temperature = minimum(pt_temp)

    # ONLY PLOT THIS TO BE CONFIDENT THAT OLLIE'S SWAP MOVE DATA LINES UP WELL WITH OURS BELOW
    # plot!(no_pt_temp, no_pt_energy, 
    #      label="Swap-Moves", 
    #      xlabel="Temperature", 
    #      yaxis="Average Energy Density, "*L"\langle\! \epsilon \rangle = \langle\! E/|\!\!E_s|\!\rangle",
    #      marker=:circle,
    #      color=alex_blue)



    # Only plot data for minimum_parallel_tempering_energy<=T<=5.0
    plot!(alex_slice_temperatures[minimum_parallel_tempering_temperature .<= alex_slice_temperatures .<= 5.0], alex_slice_average_energy_densities[minimum_parallel_tempering_temperature .<= alex_slice_temperatures .<= 5.0],
            label="Slice-Rotations", 
            # marker=:diamond,
            color=alex_green)
            
    # Only plot data for minimum_parallel_tempering_energy<=T<=5.0
    plot!(alex_swap_temperatures[minimum_parallel_tempering_temperature .<= alex_swap_temperatures .<= 5.0], alex_swap_average_energy_densities[minimum_parallel_tempering_temperature .<= alex_slice_temperatures .<= 5.0],
        label="Swap-Moves", 
         xlabel="Temperature", 
         yaxis="Average Energy Density, "*L"\langle\! \epsilon \rangle = \langle\! E/|\!\!E_s|\!\rangle",
        #  marker=:circle,
         color=alex_blue)


    # Annote L=9 Cube
    annotate!(4.3, -0.83, text("L=9 Cube", 12, :right))

    # Display the plot
    display(current())

    # Optionally, save the plot
    savefig("results/final_paper_results/parallel_tempering.png")
    savefig("results/final_paper_results/parallel_tempering.pdf")

end