using Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors, ColorSchemes
using LsqFit
using CubicSplines


include("../core/rubiks_cube.jl")
include("relaxed_anneals_figures.jl")

reduced_model(t, p) = - (t ./ p[1]) .^ p[2]




function relevant_noise_test(print_testing::Bool=false)
    temperatures = [0.9, 0.92, 0.95, 0.97, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.5, 2.0, 3.0]

    tau_10_4 = [1.6262763620132986e10, 6.085131919778271e8, 9.171847145322109e6, 4.399980428543152e6, 612999.7138109296, 40438.257709791935, 14718.74390598317, 4878.230117444054, 2093.8921279151205, 1276.954824869552, 213.61517438882305, 60.92410104070809, 28.04065774699637]
    beta_10_4 = [0.18368635969367528, 0.206711004724558, 0.2906566608784246, 0.26605419278615755, 0.30099532945745455, 0.39372792030453174, 0.4021466028955218, 0.44382798373409393, 0.4746353346496029, 0.49345565734638835, 0.5909050514963158, 0.7470285928529793, 0.8717663684870991]
    mse_10_4 = [9.324046145941833e-5, 0.00017702698347069786, 0.0002770455602538582, 0.0005633934834685305, 0.0010731495430977983, 0.0030765491397208833, 0.005483833064967181, 0.009634126256680073, 0.014263897253959913, 0.017245995215395396, 0.024630246115382677, 0.025302716663662426, 0.025569299489515123]


    tau_10_5 = [2.0414092617524675e8, 1.1423669931788402e7, 3.8529705907567455e6, 769753.5313896592, 199010.32477323478, 39170.48582864105, 15675.891130285305, 4898.924192650086, 2027.481596608481, 1262.6142920660518, 213.617888344709, 60.92410104034641, 28.04065774700506]
    beta_10_5 = [0.265874231537326, 0.323459257293533, 0.33342973830642014, 0.3802334046815314, 0.41318524903449716, 0.3837286076472755, 0.44940001823232334, 0.46529349544000836, 0.44661414693981044, 0.48335119176392227, 0.5909199945820851, 0.7470285928685986, 0.8717663684874589]
    mse_10_5 = [0.0002657514653141614, 0.0006208952515029099, 0.0010704801450796993, 0.0022028141576579258, 0.004542266003450958, 0.010634871982633735, 0.015551259684840132, 0.02144968567004492, 0.024258099781558715, 0.02495189351838861, 0.025584753793608715, 0.02548020865127239, 0.025646283830977053]

    test_MSE_values = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0]

    # First test t=10^4 values
    for i in 1:length(tau_10_4)
        println("-----------------------------------")
        println("Temperature ", temperatures[i])
        println("True Parameters (1e4): ", tau_10_4[i], ", ", beta_10_4[i])

        # Find how much MSE must increase in order for tau to be incorrect by order of magnitude
        for MSE in test_MSE_values
            p0 = [10, 0.9]
            lb = [1e-3, 0.1] # example lower bounds
            ub = [5e11, 2.0]  # example upper bounds

            model(t, p) = exp.(-(t./p[1]).^p[2])
            time_range = 1:1e4

            # Fit the model to the data
            fit = curve_fit(model, time_range, exp.(-(time_range/tau_10_4[i]).^beta_10_4[i]) .+ randn(length(time_range)).*sqrt(MSE), p0, lower=lb, upper=ub)

            # Extract the parameters
            tau_1 = fit.param[1]
            beta_1 = fit.param[2]

            if print_testing
                println("...Testing MSE: ", MSE)
                println("...Testing Estimated Parameters (1e4): ", tau_1, ", ", beta_1)
            end

            if abs(log10(tau_1) - log10(tau_10_4[i])) > 1.0
                println("Broken Estimated Parameters (1e4): ", tau_1, ", ", beta_1)
                println("Data MSE: ", mse_10_4[i])
                println("MSE to Break: ", MSE)
                println("MSE Scale Factor: ", MSE/mse_10_4[i])
                break
            end
        end
    end

    # Now test t=10^5 values
    for i in 1:length(tau_10_5)
        println("-----------------------------------")
        println("Temperature ", temperatures[i])
        println("True Parameters (1e5): ", tau_10_5[i], ", ", beta_10_5[i])

        # Find how much MSE must increase in order for tau to be incorrect by order of magnitude
        for MSE in test_MSE_values
            p0 = [10, 0.9]
            lb = [1e-3, 0.1] # example lower bounds
            ub = [5e11, 2.0]  # example upper bounds

            model(t, p) = exp.(-(t./p[1]).^p[2])
            time_range = 1:1e4

            # Fit the model to the data
            fit = curve_fit(model, time_range, exp.(-(time_range/tau_10_5[i]).^beta_10_5[i]) .+ randn(length(time_range)).*sqrt(MSE), p0, lower=lb, upper=ub)

            # Extract the parameters
            tau_1 = fit.param[1]
            beta_1 = fit.param[2]

            if print_testing
                println("...Testing MSE: ", MSE)
                println("...Testing Estimated Parameters (1e4): ", tau_1, ", ", beta_1)
            end

            if abs(log10(tau_1) - log10(tau_10_5[i])) > 1.0
                println("Broken Estimated Parameters (1e5): ", tau_1, ", ", beta_1)
                println("Data MSE: ", mse_10_5[i])
                println("MSE to Break: ", MSE)
                println("MSE Scale Factor: ", MSE/mse_10_5[i])
                break
            end
        end
    end
end


function noise_test()
    # Make some large tau stretched exponential curves with varying levels of noise and plot
    # taus = [1e7,3e7,1e8,3e8,1e9,3e9,1e10,3e10]
    # taus = [1e7,1e8,1e9,1e10, 1e11]
    taus = [1e6]
    beta = 0.25
    variance_of_noise = 1e-6

    time_range = 1:1e4

    graph = plot(title="", xlabel="Time, $(L"t")", ylabel="Stretched Exponential Functions", legend=:bottomright, ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])

    for tau in taus
        y = exp.(-(time_range./tau).^beta) + randn(length(time_range)).*sqrt(variance_of_noise)
        plot!(graph, time_range, y, label="τ = "*string(tau))

        # Fit the model to the data
        p0 = [10, 0.9]
        lb = [1e-3, 0.1] # example lower bounds
        ub = [5e11, 2.0]  # example upper bounds

        model(t, p) = exp.(-(t./p[1]).^p[2])

        # Fit the model to the data
        fit = curve_fit(model, time_range, y, p0, lower=lb, upper=ub)

        # Extract the parameters
        tau_1 = fit.param[1]
        beta_1 = fit.param[2]


        println("True Parameters: $tau, $beta")
        println("Estimated Parameters: $tau_1, $beta_1")

    end

    annotate!(graph, 1.25e5, 0.5, text("β = "*string(beta), 10))
    annotate!(graph, 1.25e5, 0.4, text("MSE = "*string(variance_of_noise), 10))
    display(graph)




end

function quick_plot()

    ### --- READ IN DATA ---
    # Read this in for all the following filenames for different temperatures
    filename = "L_11_T_0.92_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv"

    data = readdlm(joinpath("results/autocorrelation_anneal_results", filename), ',', skipstart=3)

    sample_temperature = data[1]
    samples_in_average = data[2]
    autocorrelation_functions_by_temperature = data[3:end]

    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])

    plot!(graph, 0:length(autocorrelation_functions_by_temperature)-1, autocorrelation_functions_by_temperature, label="T = "*string(sample_temperature), color=:black, linewidth=2)
    display(graph)

end

function fast_autocorrelation_function_figures(L::Int64)
    
    trials = 50
    autocorrelation_window_length=100000

    temperatures = [1.3,1.4]

    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for temperature in temperatures

        actual_number_of_trials = 0
        running_total_autocorrelation_function = zeros(autocorrelation_window_length+1-100)

        for trial in 1:trials
            println("Importing data for temperature: $(temperature), trial: $(trial)")

            filename = "results/autocorrelation_anneal_results/fast_results/"*"L_$(L)_T_$(temperature)_t_$(autocorrelation_window_length)_trial_$(trial)_slice_1.0_configuration_autocorrelation_averages_by_time.csv"
            try
                data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)

                # Skip if the matrix contains NaNs
                if any(isnan.(data_matrix))
                    continue
                end
            
                running_total_autocorrelation_function .+= data_matrix[3:end]
                
                actual_number_of_trials += 1
            catch e
                push!(filenames_that_do_not_exist, filename)
            end
        
        end

        results_dictionary[temperature] = running_total_autocorrelation_function/actual_number_of_trials
    end

    println("Filenames that do not exist: ", filenames_that_do_not_exist)

    # Plot the data
    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=:topright)

    cutoff_end = Int(3e4)
    for temperature in temperatures
        autocorrelation_function = results_dictionary[temperature]
        plot!(graph, 0:length(autocorrelation_function)-1-cutoff_end, autocorrelation_function[1:end-cutoff_end], label="T = "*string(temperature), linewidth=2)
    end

    display(graph)
    savefig(graph, "results/autocorrelation_anneal_results/fast_results/autocorrelation_averages_by_time.png")
    savefig(graph, "results/autocorrelation_anneal_results/fast_results/autocorrelation_averages_by_time.pdf")

    ## -- SAVE CONFIGURATION AVERAGES BY TIME --
    for temperature in temperatures
        filename = "L_$(L)_T_$(temperature)_t_$(autocorrelation_window_length)_1.0_configuration_autocorrelation_averages_by_time.csv"

        touch(joinpath("results/autocorrelation_anneal_results/fast_results",filename))

        open(joinpath("results/autocorrelation_anneal_results/fast_results",filename), "w") do simulation_file
            write(simulation_file, "Simulation:L=$(L), P_s=1.0, T_1=1.0, T_0=1.0, N_T=100, autocorrelation_sample_size_per_temperature=100 ,autocorrelation_window_length=$autocorrelation_window_length \n")
            write(simulation_file, "Inherent Disorder Average \n")
            write(simulation_file, "Sample Temperature T, Samples in Average N, Configuration autocorrelation average by lag \n")
            write(simulation_file, "$(temperature), 100, $(join(results_dictionary[temperature], ", ")) \n")
        end
    end
    
end


function autocorrelation_function_figures()

    simulation_name = "combined"

    ### --- READ IN DATA ---

    ## -- Define the Data --

    # Read this in for all the following filenames for different temperatures
    filenames = [
    "L_11_T_0.6_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
    # "L_11_T_0.8_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    "L_11_T_0.9_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
    # "L_11_T_0.91_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_0.92_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    # "L_11_T_0.93_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    "L_11_T_0.95_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_0.97_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.05_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.1_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.15_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.2_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.25_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.3_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.4_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_1.5_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_1.75_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_2.0_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_2.25_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_3.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_4.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_5.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_7.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_10.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv"] 

    minimal_graph_main_temperatures = [0.6, 0.9, 0.92, 0.95, 1.0, 1.2]
    minimal_graph_inset_temperatures = [0.9,0.92,0.95]


    ## -- Read in the data --

    autocorrelation_functions_by_temperature = Dict()

    # For each file in filenames, read in the data as a matrix
    for (i, filename) in pairs(filenames)
        data = readdlm(joinpath("results/autocorrelation_anneal_results/final_data", filename), ',', skipstart=3)

        temperature = data[1]
        samples_in_average = data[2]
        autocorrelation_function = data[3:end]
        
        # Store the data in the dictionary
        autocorrelation_functions_by_temperature[temperature] = autocorrelation_function
    end

    all_temperatures = sort(collect(keys(autocorrelation_functions_by_temperature)))


    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Define a fixed table of temperatures and their corresponding colors
    temp_to_color = Dict()
    alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]

    for temperature in all_temperatures
            temp_to_color[temperature] = alex_colors[mod1(findall(x -> x == temperature, all_temperatures)[1], length(alex_colors))]
    end
    temp_to_color[0.6] = alex_red
    temp_to_color[0.9] = alex_pink
    temp_to_color[0.92] = alex_blue
    temp_to_color[0.95] = alex_alt_blue
    temp_to_color[1.0] = alex_orange
    temp_to_color[1.2] = alex_green











    ### --- PLOTTING ---
    
    ## -- COMBINED GRAPHS ---
    long_fit_example_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    short_fit_example_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    minimal_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L" \bar\mathcal{C}(t)", legend=(0.9,0.35), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0], xticks=[0,2.5e4,5.0e4,7.5e4,1.0e5,1.2e5])

    c = 0.201388888888888
    long_window_length = 99900 # Int(1.2e5)
    short_window_length = Int(5e4)

    # -- FIT EXAMPLE GRAPHS --
    for temperature in all_temperatures     
        # TODO CHANGED THESE TWO
        plot!(long_fit_example_graph, 0:long_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:long_window_length] .- c), label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
        plot!(short_fit_example_graph, 0:short_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:short_window_length] .- c), label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
    end

    ## -- MINIMAL GRAPH ---
    main_minimal_graph_window_length = Int(1.2e5)
    inset_minimal_graph_window_length = Int(2.2e5)

    for temperature in minimal_graph_main_temperatures     
        # TODO CHANGED THIS
        plot!(minimal_graph, 0:main_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:main_minimal_graph_window_length] .- c), label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=3)
    end

    ### --- INSET TEMPERATURE GRAPH ON MINIMAL ---
    # Make matrix of the autocorrelation functions for the inset temperatures
    matrix_of_autocorrelation_functions = hcat([autocorrelation_functions_by_temperature[temperature][1:inset_minimal_graph_window_length] for temperature in minimal_graph_inset_temperatures]...)
    
    # plot!(minimal_graph, 0:inset_minimal_graph_window_length-1, (1/(1-c))*(matrix_of_autocorrelation_functions .- c),
    # color=[temp_to_color[inset_temperature] for inset_temperature in minimal_graph_inset_temperatures], 
    # label=["T=$(inset_temperature)" for inset_temperature in minimal_graph_inset_temperatures], inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
    # xlabel=L"t", ylabel=L"\bar\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, linewidth=2)
    
    inset_temperature = 0.9
    # TODO THESE WERE ALREADY DIFFERENT NOW CHANGED TOO
    plot!(minimal_graph, 0:inset_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[inset_temperature][1:inset_minimal_graph_window_length] .- c),
    color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
    xlabel=L"t", ylabel=L"\bar\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, linewidth=2, xticks=[])
    inset_temperature = 0.92
    plot!(minimal_graph, 0:inset_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[inset_temperature][1:inset_minimal_graph_window_length] .- c),
    color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", subplot=2,
    xlabel=L"t", ylabel=L"\bar\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, linewidth=2, xticks=[])
    inset_temperature = 0.95
    plot!(minimal_graph, 0:inset_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[inset_temperature][1:inset_minimal_graph_window_length] .- c),
    color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", subplot=2,
    xlabel=L"t", ylabel=L"\bar\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, linewidth=2, xticks=[0,1e5, 2.2e5], legendfontsize=6)

    # for inset_temperature in minimal_graph_inset_temperatures
        # plot!(minimal_graph, 0:inset_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[inset_temperature][1:inset_minimal_graph_window_length] .- c),
        # color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
        # xlabel=L"t", ylabel=L"\bar\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, linewidth=2)
    # end



    


    ### --- FITTING PARAMETER EXTRACTION ---
    excluded_temperatures_from_fitting = [0.6,0.8]
    fitted_temperatures = [T for T in all_temperatures if T ∉ excluded_temperatures_from_fitting]

    tau_beta_by_temperature = Float64[]
    beta_by_temperature = Float64[]
    tau_errors_by_temperature = Float64[]  # Add this line
    beta_errors_by_temperature = Float64[]  # Add beta error array

    longer_tau_beta_by_temperature = Float64[]
    longer_beta_by_temperature = Float64[]
    longer_tau_errors_by_temperature = Float64[]  # Add this line
    longer_beta_errors_by_temperature = Float64[]  # Add longer beta error array

    for temperature in fitted_temperatures
        if temperature in excluded_temperatures_from_fitting
            continue
        end
        println("Temperature: ", temperature)


        # Fit the data to the autocorrelation function
        p0 = [10, 0.9]
        lb = [1e-3, 0.1] # example lower bounds
        ub = [1e11, 2.0]  # example upper bounds

        beta_relaxation_end = 1e4
        longer_beta_relaxation_end = 99900

        t_values = 1:beta_relaxation_end
        longer_t_values = 1:longer_beta_relaxation_end

        # Define the model function
        c = 0.201388888888888
        # TODO DIDN'T CHANGE THIS BUT NOW THIS IS CORRECT
        autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

        # TODO IS THIS NEEDED NOW WE'VE DONE CHANGES?
        if temperature >= 1.2

            # Fit the model to the data
            autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

            # Extract the parameters for t=10^4 fit
            tau = autocorrelation_fit.param[1]
            beta = autocorrelation_fit.param[2]
            
            # Extract standard errors
            param_errors = stderror(autocorrelation_fit)
            tau_error = param_errors[1]
            # Convert tau error to log10(tau) error: d(log10(tau))/dtau = 1/(tau*ln(10))
            log_tau_error = tau_error / (tau * log(10))

            
            println("t=10^4 Tau = ", tau)
            println("t=10^4 Beta = ", beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [tau, beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^4 Mean Squared Error = ", mse)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
            beta_by_temperature = append!(beta_by_temperature, beta)
            tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
            beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



            # Extract the parameters for t=10^5 fit
            longer_autocorrelation_fit = curve_fit(autocorrelation_model, longer_t_values, autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)], p0, lower=lb, upper=ub)
            longer_tau = longer_autocorrelation_fit.param[1]
            longer_beta = longer_autocorrelation_fit.param[2]
            
            # Extract standard errors for longer fit
            longer_param_errors = stderror(longer_autocorrelation_fit)
            longer_tau_error = longer_param_errors[1]
            # Convert tau error to log10(tau) error
            longer_log_tau_error = longer_tau_error / (longer_tau * log(10))

            
            println("t=10^5 Tau = ", longer_tau)
            println("t=10^5 Beta = ", longer_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_tau, longer_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^5 Mean Squared Error = ", mse)

            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_beta)
            longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
            longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])


        else
            # Now try fitting on log of data with offset already removed
            # TODO CHANGED THIS
            reduced_data = log.((1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- c))


            autocorrelation_fit = curve_fit(reduced_model, t_values, reduced_data, p0, lower=lb, upper=ub)

            reduced_tau = autocorrelation_fit.param[1]
            reduced_beta = autocorrelation_fit.param[2]
            
            # Extract standard errors for reduced model
            param_errors = stderror(autocorrelation_fit)
            tau_error = param_errors[1]
            # For reduced model, tau error is already in appropriate scale
            log_tau_error = tau_error / (reduced_tau * log(10))
            
            println("t=10^4 Tau From Logged Data = ", reduced_tau)
            println("t=10^4 Beta from Logged Data = ", reduced_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [reduced_tau, reduced_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^4 Mean Squared Error = ", mse)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, reduced_tau)
            beta_by_temperature = append!(beta_by_temperature, reduced_beta)
            tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
            beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



            # Extract the parameters for t=10^5 fit
            # TODO CHANGED THIS
            longer_reduced_data = log.((1/(1-c)) * (autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- c))

            longer_autocorrelation_fit = curve_fit(reduced_model, longer_t_values, longer_reduced_data, p0, lower=lb, upper=ub)
            longer_reduced_tau = longer_autocorrelation_fit.param[1]
            longer_reduced_beta = longer_autocorrelation_fit.param[2]
            
            # Extract standard errors for longer reduced model
            longer_param_errors = stderror(longer_autocorrelation_fit)
            longer_tau_error = longer_param_errors[1]
            longer_log_tau_error = longer_tau_error / (longer_reduced_tau * log(10))
            
            println("t=10^5 Tau From Logged Data = ", longer_reduced_tau)
            println("t=10^5 Beta from Logged Data = ", longer_reduced_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_reduced_tau, longer_reduced_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^5 Mean Squared Error = ", mse)
            
            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_reduced_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_reduced_beta)
            longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
            longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])
        end


        # Plot the fit
        dashed_cutoff = beta_relaxation_end
        extrapolated_dashed_cutoff = 5e4

        if temperature in fitted_temperatures
            # label="T = "*string(temperature)*", τ = "*string(tau_beta)*", β = "*string(beta)
            # TODO CHANGED ALL THESE
            plot!(long_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")
            plot!(short_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(extrapolated_dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")

            if temperature in minimal_graph_main_temperatures
                plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="", linewidth=1.5)
            end
        end

        # Plot the longer fit
        longer_dashed_cutoff = longer_beta_relaxation_end
        if temperature in fitted_temperatures
            # TODO CHANGED ALL THESE
            plot!(long_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="", legend=false)
            plot!(short_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(extrapolated_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="", legend=false)
        
            if temperature in minimal_graph_main_temperatures
                plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="", linewidth=1.5)
            end
        end

    end

    # --- SAVE GRAPHS ---
    savefig(long_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.png"))
    savefig(long_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.pdf"))
    display(long_fit_example_graph)

    savefig(short_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.png"))
    savefig(short_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.pdf"))
    display(short_fit_example_graph)

    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.png"))
    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.pdf"))
    display(minimal_graph)






    if simulation_name == "combined"
        # -- RELAXATION TIME FITS ---
        # Make plot with tau_beta_by_temperature against temperature 
        # and then Arrhenius/VFT/parabolic fit curves to data too and output measure
        # Do fits on log of data to avoid numerical errors
        log_tau_beta_by_temperature = log10.(tau_beta_by_temperature) # Note using log10 here
        longer_log_tau_beta_by_temperature = log10.(longer_tau_beta_by_temperature)
        
        combined_fitted_temperatures = vcat(fitted_temperatures, fitted_temperatures+1e-15*randn(length(fitted_temperatures)))
        combined_log_tau_beta_by_temperature = vcat(log_tau_beta_by_temperature, longer_log_tau_beta_by_temperature)
        # Sort combined_fitted_temperatures but sort combined_log_tau_beta_by_temperature in the same way
        combined_fitted_temperatures, combined_log_tau_beta_by_temperature = combined_fitted_temperatures[sortperm(combined_fitted_temperatures)], combined_log_tau_beta_by_temperature[sortperm(combined_fitted_temperatures)]

        fit_temperatures = collect(LinRange(minimum(fitted_temperatures)-0.03, maximum(fitted_temperatures), 100))

        tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3), xticks=[0,1,2,3,4,5,6,7,8,9,10], yticks=[0,1,2,3,4,5,6,7,8,9,10], xlims=(0,10.5), ylims=(0,10))
        squares_only_tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3), xticks=[0,1,2,3,4,5,6,7,8,9,10], yticks=[0,1,2,3,4,5,6,7,8,9,10], ylims=(0,10.5), xlims=(0,10.5))
        circles_only_tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3), xticks=[0,1,2,3,4,5,6,7,8,9,10], yticks=[0,1,2,3,4,5,6,7,8,9,10], ylims=(0,10.5), xlims=(0,10.5))

        scatter!(tau_fits_graph, fitted_temperatures, log_tau_beta_by_temperature, 
                 yerror=tau_errors_by_temperature, color=alex_red, label="", markersize=5,
                 ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], alpha=0.8)
        scatter!(circles_only_tau_fits_graph, fitted_temperatures, log_tau_beta_by_temperature, 
                 yerror=tau_errors_by_temperature, color=alex_red, label="", markersize=5, marker=:circle, 
                 alpha=0.8, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5])

        scatter!(tau_fits_graph, fitted_temperatures, longer_log_tau_beta_by_temperature, 
                 yerror=longer_tau_errors_by_temperature, color=alex_red, label="", markersize=4, 
                 marker=:square, alpha=0.8)
        scatter!(squares_only_tau_fits_graph, fitted_temperatures, longer_log_tau_beta_by_temperature, 
                 yerror=longer_tau_errors_by_temperature, color=alex_red, label="", markersize=4, 
                 marker=:square, alpha=0.8, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5])


        # # Arrhenius fit
        # arrhenius_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T), all_temperatures, log_tau_beta_by_temperature, [0.1, 0.1], lower=[-Inf,-Inf], upper=[Inf, Inf])
        # arrhenius_fit_parameters = arrhenius_fit.param
        # arrhenius_fit_curve = arrhenius_fit_parameters[1] .+ (arrhenius_fit_parameters[2] ./ fit_temperatures)
        # lowest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][1] .+ (confint(arrhenius_fit, level=conf_int)[2][1] ./ fit_temperatures)
        # highest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][2] .+ (confint(arrhenius_fit, level=conf_int)[2][2] ./ fit_temperatures)
        # plot!(tau_fits_graph, fit_temperatures, arrhenius_fit_curve, label="Arrhenius Fit", color=alex_orange, ribbon=(arrhenius_fit_curve .- lowest_confidence_interval_arrhenius_curve, highest_confidence_interval_arrhenius_curve .- arrhenius_fit_curve), fc=:orange, fa=0.05)

        # println("Arrhenius Parameters: ", arrhenius_fit_parameters)
        # println("Arrhenius Residuals", arrhenius_fit.resid)
        # println("Arrhenius Covariance Matrix", vcov(arrhenius_fit))
        # println("Arrhenius Standard Errors", stderror(arrhenius_fit))
        # println("Arrhenius Margin Error: ", margin_error(arrhenius_fit))
        # println("Confdience Intervals", confint(arrhenius_fit))
        # println("")
        # println("")



        # Parabolic fit over all data
        parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        parabolic_fit_parameters = parabolic_fit.param
        parabolic_fit_curve = parabolic_fit_parameters[1] .+ (parabolic_fit_parameters[2] ./ fit_temperatures) .+ (parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(circles_only_tau_fits_graph, fit_temperatures, parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("")
        println("")
        println("Shorter Parabolic Parameters: ", parabolic_fit_parameters)
        # println("Shorter Parabolic Residuals", parabolic_fit.resid)
        # println("Shorter Parabolic Covariance Matrix", vcov(parabolic_fit))
        # println("Shorter Parabolic Standard Errors", stderror(parabolic_fit))
        # println("Shorter Parabolic Margin Error: ", margin_error(parabolic_fit))
        # println("Shorter Confdience Intervals", confint(parabolic_fit))

        longer_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, longer_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        longer_parabolic_fit_parameters = longer_parabolic_fit.param
        longer_parabolic_fit_curve = longer_parabolic_fit_parameters[1] .+ (longer_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(squares_only_tau_fits_graph, fit_temperatures, longer_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("Longer Parabolic Parameters: ", longer_parabolic_fit_parameters)
        # println("Longer Parabolic Residuals", longer_parabolic_fit.resid)
        # println("Longer Parabolic Covariance Matrix", vcov(longer_parabolic_fit))
        # println("Longer Parabolic Standard Errors", stderror(longer_parabolic_fit))
        # println("Longer Parabolic Margin Error: ", margin_error(longer_parabolic_fit))
        # println("Longer Confdience Intervals", confint(longer_parabolic_fit))

        combined_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        combined_parabolic_fit_parameters = combined_parabolic_fit.param
        combined_parabolic_fit_curve = combined_parabolic_fit_parameters[1] .+ (combined_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (combined_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(tau_fits_graph, fit_temperatures, combined_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("Combined Parabolic Parameters: ", combined_parabolic_fit_parameters)
        # println("Combined Parabolic Residuals", combined_parabolic_fit.resid)
        # println("Combined Parabolic Covariance Matrix", vcov(combined_parabolic_fit))
        # println("Combined Parabolic Standard Errors", stderror(combined_parabolic_fit))
        # println("Combined Parabolic Margin Error: ", margin_error(combined_parabolic_fit))
        # println("Combined Confdience Intervals", confint(combined_parabolic_fit))



        # Parabolic data over temperatures < 1.4
        alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures[1:findlast(x -> x < 1.4, fitted_temperatures)], log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, fitted_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        alt_parabolic_fit_parameters = alt_parabolic_fit.param
        alt_parabolic_fit_curve = alt_parabolic_fit_parameters[1] .+ (alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(circles_only_tau_fits_graph, fit_temperatures, alt_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_orange, linewidth=3)

        println("")
        println("")
        println("Alternative Parabolic Parameters: ", alt_parabolic_fit_parameters)
        # println("Alternative Parabolic Residuals", alt_parabolic_fit.resid)
        # println("Alternative Parabolic Covariance Matrix", vcov(alt_parabolic_fit))
        # println("Alternative Parabolic Standard Errors", stderror(alt_parabolic_fit))
        # println("Alternative Parabolic Margin Error: ", margin_error(alt_parabolic_fit))
        # println("Alternative Parabolic Confidence Intervals", confint(alt_parabolic_fit))

        longer_alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures[1:findlast(x -> x < 1.4, fitted_temperatures)], longer_log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, fitted_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        longer_alt_parabolic_fit_parameters = longer_alt_parabolic_fit.param
        longer_alt_parabolic_fit_curve = longer_alt_parabolic_fit_parameters[1] .+ (longer_alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(squares_only_tau_fits_graph, fit_temperatures, longer_alt_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_orange, linewidth=3)

        print("Longer Alternative Parabolic Parameters: ", longer_alt_parabolic_fit_parameters)
        # println("Longer Alternative Parabolic Residuals", longer_alt_parabolic_fit.resid)
        # println("Longer Alternative Parabolic Covariance Matrix", vcov(longer_alt_parabolic_fit))
        # println("Longer Alternative Parabolic Standard Errors", stderror(longer_alt_parabolic_fit))
        # println("Longer Alternative Parabolic Margin Error: ", margin_error(longer_alt_parabolic_fit))
        # println("Longer Alternative Parabolic Confidence Intervals", confint(longer_alt_parabolic_fit))

        combined_alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), combined_fitted_temperatures[1:findlast(x -> x < 1.4, combined_fitted_temperatures)], combined_log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, combined_fitted_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        combined_alt_parabolic_fit_parameters = combined_alt_parabolic_fit.param
        new_fit_temperatures = collect(LinRange(minimum(fitted_temperatures)-0.03, 1.7, 100))
        combined_alt_parabolic_fit_curve = combined_alt_parabolic_fit_parameters[1] .+ (combined_alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (combined_alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        combined_alt_parabolic_fit_curve_new_temperatures = combined_alt_parabolic_fit_parameters[1] .+ (combined_alt_parabolic_fit_parameters[2] ./ new_fit_temperatures) .+ (combined_alt_parabolic_fit_parameters[3] ./ new_fit_temperatures.^2)
        plot!(tau_fits_graph, new_fit_temperatures, combined_alt_parabolic_fit_curve_new_temperatures, label="Parabolic Fit, "*L"T \leq 1.4", linestyle=:dash, color=alex_orange, linewidth=3, legend=(0.75, 0.3))

        print("Combined Alternative Parabolic Parameters: ", combined_alt_parabolic_fit_parameters)
        # println("Combined Alternative Parabolic Residuals", combined_alt_parabolic_fit.resid)
        # println("Combined Alternative Parabolic Covariance Matrix", vcov(combined_alt_parabolic_fit))
        # println("Combined Alternative Parabolic Standard Errors", stderror(combined_alt_parabolic_fit))
        # println("Combined Alternative Parabolic Margin Error: ", margin_error(combined_alt_parabolic_fit))
        # println("Combined Alternative Parabolic Confidence Intervals", confint(combined_alt_parabolic_fit))




        # VFT fit
        vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        vft_fit_parameters = vft_fit.param
        vft_fit_curve = vft_fit_parameters[1] .+ (vft_fit_parameters[2] ./ (fit_temperatures .- vft_fit_parameters[3]))
        plot!(circles_only_tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        println("")
        println("")
        println("VFT Parameters: ", vft_fit_parameters)
        # println("VFT Residuals", vft_fit.resid)
        # println("VFT Covariance Matrix", vcov(vft_fit))
        # println("VFT Standard Errors", stderror(vft_fit))
        # println("VFT Margin Error: ", margin_error(vft_fit))
        # println("Confdience Intervals", confint(vft_fit))


        longer_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, longer_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        longer_vft_fit_parameters = longer_vft_fit.param
        longer_vft_fit_curve = longer_vft_fit_parameters[1] .+ (longer_vft_fit_parameters[2] ./ (fit_temperatures .- longer_vft_fit_parameters[3]))
        plot!(squares_only_tau_fits_graph, fit_temperatures, longer_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        print("Longer VFT Parameters: ", longer_vft_fit_parameters)
        # println("Longer VFT Residuals", longer_vft_fit.resid)
        # println("Longer VFT Covariance Matrix", vcov(longer_vft_fit))
        # println("Longer VFT Standard Errors", stderror(longer_vft_fit))
        # println("Longer VFT Margin Error: ", margin_error(longer_vft_fit))
        # println("Longer Confdience Intervals", confint(longer_vft_fit))

        combined_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        combined_vft_fit_parameters = combined_vft_fit.param
        combined_vft_fit_curve = combined_vft_fit_parameters[1] .+ (combined_vft_fit_parameters[2] ./ (fit_temperatures .- combined_vft_fit_parameters[3]))
        plot!(tau_fits_graph, fit_temperatures, combined_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        print("Combined VFT Parameters: ", combined_vft_fit_parameters)
        # println("Combined VFT Residuals", combined_vft_fit.resid)
        # println("Combined VFT Covariance Matrix", vcov(combined_vft_fit))
        # println("Combined VFT Standard Errors", stderror(combined_vft_fit))
        # println("Combined VFT Margin Error: ", margin_error(combined_vft_fit))
        # println("Combined Confdience Intervals", confint(combined_vft_fit))


        #--- FRAGILITY INSET DATA ---
        scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], log_tau_beta_by_temperature,
        yerror=tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2,
        xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:circle, alpha=0.8)
        scatter!(circles_only_tau_fits_graph, [1/T for T in fitted_temperatures], log_tau_beta_by_temperature,
        yerror=tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2,
        xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:circle, alpha=0.8)

        scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], longer_log_tau_beta_by_temperature,
        yerror=longer_tau_errors_by_temperature, color=alex_red, legend=false, subplot=2, xlabel=L"1/T",
        ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:square, alpha=0.8)
        scatter!(squares_only_tau_fits_graph, [1/T for T in fitted_temperatures], longer_log_tau_beta_by_temperature,
        yerror=longer_tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2, xlabel=L"1/T",
        ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:square, alpha=0.8)


        # Plot fits in inset too
        # plot!(tau_fits_graph, [1/T for T in fit_temperatures], arrhenius_fit_curve,color=alex_orange, label="Arrhenius Fit", linestyle=:dash, subplot=2)
        
        plot!(circles_only_tau_fits_graph, [1/T for T in fit_temperatures], parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(circles_only_tau_fits_graph, [1/T for T in fit_temperatures], alt_parabolic_fit_curve,color=alex_orange, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(circles_only_tau_fits_graph, [1/T for T in fit_temperatures], vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)

        plot!(squares_only_tau_fits_graph, [1/T for T in fit_temperatures], longer_parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(squares_only_tau_fits_graph, [1/T for T in fit_temperatures], longer_alt_parabolic_fit_curve,color=alex_orange, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(squares_only_tau_fits_graph, [1/T for T in fit_temperatures], longer_vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


        # Combined graph
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_alt_parabolic_fit_curve,color=alex_orange, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


        ### EXTRA INSET STUFF

        ## Curve from all energy barrier data
        # temperatures = [10.0, 9.54992586021436, 9.120108393559098, 8.709635899560805,8.317637711026709, 7.943282347242815, 7.5857757502918375, 7.244359600749901,6.918309709189366, 6.60693448007596, 6.3095734448019325, 6.025595860743578,5.7543993733715695, 5.495408738576245, 5.248074602497725, 5.011872336272722,4.786300923226383, 4.57088189614875, 4.36515832240166, 4.168693834703355,3.981071705534972, 3.801893963205612, 3.6307805477010135, 3.467368504525316,3.311311214825911, 3.1622776601683795, 3.019951720402016, 2.8840315031266055,2.7542287033381663, 2.6302679918953817, 2.51188643150958, 2.3988329190194904,2.2908676527677727, 2.1877616239495525, 2.089296130854039, 1.9952623149688797,1.9054607179632475, 1.8197008586099832, 1.7378008287493754, 1.6595869074375604,1.5848931924611134, 1.5135612484362082, 1.4454397707459277, 1.380384264602885,1.318256738556407, 1.2589254117941673, 1.202264434617413, 1.148153621496883,1.096478196143185, 1.0471285480508996, 1.0, 0.9549925860214359,0.9120108393559097, 0.8709635899560806, 0.831763771102671, 0.7943282347242814,0.7585775750291835, 0.7244359600749903, 0.6918309709189365, 0.6606934480075961,0.6309573444801932, 0.6025595860743578, 0.5754399373371569, 0.5495408738576245,0.5248074602497725, 0.5011872336272722, 0.4786300923226383, 0.45708818961487496,0.4365158322401659, 0.4168693834703355, 0.3981071705534973, 0.38018939632056126,0.36307805477010135, 0.34673685045253166, 0.33113112148259116, 0.3162277660168379,0.3019951720402016, 0.28840315031266056, 0.2754228703338166, 0.26302679918953814,0.25118864315095796, 0.23988329190194896, 0.22908676527677738, 0.2187761623949553,0.20892961308540398, 0.19952623149688797, 0.19054607179632474, 0.18197008586099836,0.17378008287493754, 0.16595869074375605, 0.15848931924611134, 0.1513561248436208,0.1445439770745927, 0.13803842646028847, 0.13182567385564076, 0.12589254117941676,0.12022644346174131, 0.11481536214968828, 0.1096478196143185, 0.10471285480508996]
        # average_energy_densities_by_temperature = [-0.1809469696969697, -0.18082575757575758, -0.18192424242424243, -0.17976515151515152,-0.1805378787878788, -0.18525, -0.18692424242424244, -0.18406818181818183,-0.18787121212121213, -0.1863030303030303, -0.18712121212121213, -0.19046969696969698,-0.19040151515151515, -0.19242424242424241, -0.19165151515151516, -0.1923030303030303,-0.19496212121212123, -0.19711363636363635, -0.20083333333333336, -0.20047727272727273,-0.20236363636363636, -0.20400757575757578, -0.20588636363636362, -0.20875000000000002,-0.21309848484848487, -0.2115757575757576, -0.2155, -0.21622727272727274,-0.21995454545454549, -0.2198106060606061, -0.2242121212121212, -0.22921969696969696,-0.2314090909090909, -0.2324848484848485, -0.23884090909090908, -0.24062121212121212,-0.24951515151515152, -0.2516515151515152, -0.258530303030303, -0.2602878787878788,-0.26795454545454545, -0.274, -0.28042424242424246, -0.2888939393939394,-0.2936590909090909, -0.30199242424242423, -0.31156060606060604, -0.3235227272727273,-0.33324242424242423, -0.34381060606060604, -0.3576212121212121, -0.36961363636363637,-0.3893030303030303, -0.4058257575757576, -0.43169696969696975, -0.4528560606060606,-0.48315909090909087, -0.5125, -0.5389015151515152, -0.5718484848484848,-0.5935530303030303, -0.612, -0.6338863636363636, -0.6457954545454546,-0.6553484848484848, -0.6639848484848485, -0.6690681818181818, -0.6777272727272727,-0.686060606060606, -0.6898863636363636, -0.6929015151515151, -0.6951590909090909,-0.6950303030303031, -0.6968106060606061, -0.6981060606060606, -0.6991893939393939,-0.7018560606060606, -0.7027348484848485, -0.7056742424242425, -0.7068560606060607,-0.7070454545454546, -0.707780303030303, -0.7078636363636364, -0.7083181818181818,-0.7086287878787879, -0.7083181818181818, -0.708939393939394, -0.7088181818181818,-0.7088030303030303, -0.7090151515151515, -0.7089772727272727, -0.7090681818181819,-0.7090681818181819, -0.7090757575757576, -0.7090681818181819, -0.7090454545454545,-0.7090757575757576, -0.7090909090909091, -0.7090909090909091, -0.7090909090909091]
        # average_energies_by_temperature = average_energy_densities_by_temperature .* 1320
    
        # # Make spline from temperatures to energy
        # average_energies_by_temperature, temperatures = average_energies_by_temperature[sortperm(average_energies_by_temperature)], temperatures[sortperm(average_energies_by_temperature)]
        # spline_temperature_from_energy = CubicSpline(average_energies_by_temperature, temperatures, extrapl=[1,], extrapr=[1,])
        # println("Temperature at E=-660: ", spline_temperature_from_energy(-660))

        # # Make spline from energy to temperature
        # average_energies_by_temperature, temperatures = average_energies_by_temperature[sortperm(temperatures)], temperatures[sortperm(temperatures)]
        # spline_energy_from_temperature = CubicSpline(temperatures, average_energies_by_temperature, extrapl=[1,], extrapr=[1,])
        # println("Energy at T=1: ", spline_energy_from_temperature(1.0))
    
        # # Make prediction function y(T) = 88((spline_energy_from_temperature(T)/1320) - (1/6))/T)
        # C = -12.5
        # prediction_function = T -> (88 * ((abs(spline_energy_from_temperature(T))/1320) - (1/6))/T) + C

        # # Plot the prediction function on inset®
        # plot!(tau_fits_graph, [1/T for T in fit_temperatures], prediction_function.(fit_temperatures), label="Prediction Function", linestyle=:dash, color=:black, subplot=2)

        
        # # Plot exponential of prediction function on main graph
        # plot!(tau_fits_graph, fit_temperatures, prediction_function.(fit_temperatures), linestyle=:dash, color=:black)

        ## Asymptotic Line from E^* energy barrier
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], (26.4 .* [1/T for T in fit_temperatures]) .- 21, linestyle=:dash, color=:black, linewidth=3, subplot=2)



        

        # --- BETA INSET DATA ---
        # Plot longer beta data
        scatter!(tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
        yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T", 
        ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])
        scatter!(squares_only_tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
        yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T",
        ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])

        # Plot the shorter beta data
        scatter!(tau_fits_graph, fitted_temperatures, beta_by_temperature,
        yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  subplot=3)
        scatter!(circles_only_tau_fits_graph, fitted_temperatures, beta_by_temperature, ylabel=L"\beta",
        yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])






        # STRECHED EXPONENTIAL DEFINITION STUFF --
        # Create spline for all beta data against temperature
        beta_by_temperature = longer_beta_by_temperature
        temperatures_for_beta = fitted_temperatures
        temperature_from_beta_spline = CubicSpline(beta_by_temperature[sortperm(beta_by_temperature)], temperatures_for_beta[sortperm(beta_by_temperature)], extrapl=[1,], extrapr=[1,])
        beta_from_temperature_spline = CubicSpline(temperatures_for_beta[sortperm(temperatures_for_beta)], beta_by_temperature[sortperm(temperatures_for_beta)], extrapl=[1,], extrapr=[1,])
        # Plot spline
        # plot!(tau_fits_graph, temperature_from_beta_spline(beta_by_temperature), beta_by_temperature, color=alex_orange, label="Spline", subplot=3)
        # Onset of stretching definition
        # beta_stretch = 0.8*beta_from_temperature_spline(10.0)
        beta_infinite_temperature_dict = Dict(5 => 1.1812820380521043, 7 => 1.0838099849875709, 11 => 1.0229843232478562, 9 => 1.0370097148780135)
        beta_stretch = 0.75*beta_infinite_temperature_dict[11]

        # Find temperature where beta = beta_stretch
        temperature_at_beta_stretch = temperature_from_beta_spline(beta_stretch)
        println("L=11, Temperature at β = $beta_stretch: ", temperature_at_beta_stretch)
        # Plot hline at \beta=beta_stregch
        hline!(tau_fits_graph, [beta_stretch], color=alex_red, label="β = $(beta_stretch)", linestyle=:dash, linewidth=2, subplot=3)
        # Plot vline at temperature_at_beta_stretch
        vline!(tau_fits_graph, [temperature_at_beta_stretch], color=alex_red, label=L"\bar T^{\rm on}", linestyle=:dash, linewidth=2, subplot=3)
        annotate!(tau_fits_graph, [(temperature_at_beta_stretch+0.2, 0.1, text(L"\bar T^{\rm on}", 10, :left, color=alex_red))], subplot=3)
        

        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.png"))
        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.pdf"))
        display(tau_fits_graph)

        savefig(circles_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_circles_only.png"))
        savefig(circles_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_circles_only.pdf"))
        display(circles_only_tau_fits_graph)

        savefig(squares_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_squares_only.png"))
        savefig(squares_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_squares_only.pdf"))
        display(squares_only_tau_fits_graph)

    end
end







function swap_autocorrelation_function_figures()

    simulation_name = "combined_swap"

    ### --- READ IN DATA ---

    ## -- Define the Data --

    # Read this in for all the following filenames for different temperatures
    filenames = [
    "swap_L_11_T_0.6_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_0.9_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_1.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_3.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",] 


    ## -- Read in the data --

    autocorrelation_functions_by_temperature = Dict()

    # For each file in filenames, read in the data as a matrix
    for (i, filename) in pairs(filenames)
        data = readdlm(joinpath("results/autocorrelation_anneal_results/final_data", filename), ',', skipstart=3)

        temperature = data[1]
        samples_in_average = data[2]
        autocorrelation_function = data[3:end]
        
        # Store the data in the dictionary
        autocorrelation_functions_by_temperature[temperature] = autocorrelation_function
    end

    all_temperatures = sort(collect(keys(autocorrelation_functions_by_temperature)))


    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Define a fixed table of temperatures and their corresponding colors
    temp_to_color = Dict()
    alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]

    for temperature in all_temperatures
            temp_to_color[temperature] = alex_colors[mod1(findall(x -> x == temperature, all_temperatures)[1], length(alex_colors))]
    end
    temp_to_color[0.6] = alex_red
    temp_to_color[0.9] = alex_pink
    temp_to_color[0.92] = alex_blue
    temp_to_color[0.95] = alex_alt_blue
    temp_to_color[1.0] = alex_orange
    temp_to_color[1.2] = alex_green








    ### --- PLOTTING ---
    
    ## -- COMBINED GRAPHS ---
    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.35), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])

    c = 0.201388888888888

    tau_fits_graph = plot(title="Relaxation Time Fits", xlabel="Temperature, $(L"T")", ylabel=L"\log_{10}(\tau)")
    beta_fits_graph = plot(title="Beta Fits", xlabel="Temperature, $(L"T")", ylabel=L"\beta")


    for temperature in all_temperatures
        autocorrelation_function = autocorrelation_functions_by_temperature[temperature]

        plot!(graph, (1/(1-c))*(autocorrelation_function .- c), label="T = $temperature", color=temp_to_color[temperature], linewidth=4, legend=:topright)
    

        # Exponential fit
        p0 = [10, 0.9]
        lb = [1e-3, 0.1] # example lower bounds
        ub = [1e11, 2.0]  # example upper bounds
        autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

        t_values = collect(1:length(autocorrelation_function))
        autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_function, p0, lower=lb, upper=ub)
        autocorrelation_fit_curve = autocorrelation_model(t_values, autocorrelation_fit.param)
        # unrescaled_exponential_fit_curve = (1-c)*exponential_fit_curve .+ c

        plot!(graph, (1/(1-c))*(autocorrelation_fit_curve .-c), label="", color=:black, linestyle=:dash, linewidth=2, margin=5mm)


        println("Temperature: ", temperature)
        println("Stretched Fit Parameters: ")
        println("Tau:: ", autocorrelation_fit.param[1])
        println("Beta:: ", autocorrelation_fit.param[2])

        scatter!(tau_fits_graph, [temperature], [log10(autocorrelation_fit.param[1])], label="", color=temp_to_color[temperature], marker=:circle, markersize=5)
        scatter!(beta_fits_graph, [temperature], [autocorrelation_fit.param[2]], label="", color=temp_to_color[temperature], marker=:circle, markersize=5)

    end

    savefig(graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_functions.png"))
    savefig(graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_functions.pdf"))
    display(graph)

    display(tau_fits_graph)
    display(beta_fits_graph)







end

























function full_swap_autocorrelation_function_figures()

    simulation_name = "combined_swap"

    ### --- READ IN DATA ---

    ## -- Define the Data --

    # Read this in for all the following filenames for different temperatures
    filenames = [
    "swap_L_11_T_0.6_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_0.8_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_0.9_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_0.92_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_0.95_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_1.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_1.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_2.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",  
    "swap_L_11_T_3.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_3.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_4.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_4.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "swap_L_11_T_5.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",] 


    minimal_graph_main_temperatures = [0.6, 0.8, 0.9, 0.92, 0.95, 1.0, 1.5, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]


    ## -- Read in the data --

    autocorrelation_functions_by_temperature = Dict()

    # For each file in filenames, read in the data as a matrix
    for (i, filename) in pairs(filenames)
        data = readdlm(joinpath("results/autocorrelation_anneal_results", filename), ',', skipstart=3)

        temperature = data[1]
        samples_in_average = data[2]
        autocorrelation_function = data[3:end]
        
        # Store the data in the dictionary
        autocorrelation_functions_by_temperature[temperature] = autocorrelation_function
    end

    all_temperatures = sort(collect(keys(autocorrelation_functions_by_temperature)))


    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Define a fixed table of temperatures and their corresponding colors
    temp_to_color = Dict()
    alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]

    for temperature in all_temperatures
            temp_to_color[temperature] = alex_colors[mod1(findall(x -> x == temperature, all_temperatures)[1], length(alex_colors))]
    end
    temp_to_color[0.6] = alex_red
    temp_to_color[0.9] = alex_pink
    temp_to_color[0.92] = alex_blue
    temp_to_color[0.95] = alex_alt_blue
    temp_to_color[1.0] = alex_orange
    temp_to_color[1.2] = alex_green











    ### --- PLOTTING ---
    
    ## -- COMBINED GRAPHS ---
    long_fit_example_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    short_fit_example_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    minimal_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.35), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0], xticks=[0,2.5e4,5.0e4,7.5e4,1.0e5,1.2e5])

    c = 0.201388888888888
    long_window_length = 99900 # Int(1.2e5)
    short_window_length = Int(5e4)

    # -- FIT EXAMPLE GRAPHS --
    for temperature in all_temperatures     
        # TODO CHANGED THESE TWO
        plot!(long_fit_example_graph, 0:long_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:long_window_length] .- c), label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
        plot!(short_fit_example_graph, 0:short_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:short_window_length] .- c), label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
    end

    ## -- MINIMAL GRAPH ---
    main_minimal_graph_window_length = Int(0.999e5)

    for temperature in minimal_graph_main_temperatures     
        # TODO CHANGED THIS
        plot!(minimal_graph, 0:main_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:main_minimal_graph_window_length] .- c), label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
    end


    


    ### --- FITTING PARAMETER EXTRACTION ---
    excluded_temperatures_from_fitting = []
    fitted_temperatures = [T for T in all_temperatures if T ∉ excluded_temperatures_from_fitting]

    tau_beta_by_temperature = []
    beta_by_temperature = []

    longer_tau_beta_by_temperature = []
    longer_beta_by_temperature = []

    for temperature in fitted_temperatures
        if temperature in excluded_temperatures_from_fitting
            continue
        end
        println("Temperature: ", temperature)


        # Fit the data to the autocorrelation function
        p0 = [10, 0.9]
        lb = [1e-3, 0.1] # example lower bounds
        ub = [1e11, 2.0]  # example upper bounds

        beta_relaxation_end = 1e4
        longer_beta_relaxation_end = 99900

        t_values = 1:beta_relaxation_end
        longer_t_values = 1:longer_beta_relaxation_end

        # Define the model function
        c = 0.201388888888888
        # TODO DIDN'T CHANGE THIS BUT NOW THIS IS CORRECT
        autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

        # TODO IS THIS NEEDED NOW WE'VE DONE CHANGES?
        if temperature >= 0.4

            # Fit the model to the data
            autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

            # Extract the parameters for t=10^4 fit
            tau = autocorrelation_fit.param[1]
            beta = autocorrelation_fit.param[2]
            println("t=10^4 Tau = ", tau)
            println("t=10^4 Beta = ", beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [tau, beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^4 Mean Squared Error = ", mse)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
            beta_by_temperature = append!(beta_by_temperature, beta)



            # Extract the parameters for t=10^5 fit
            longer_autocorrelation_fit = curve_fit(autocorrelation_model, longer_t_values, autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)], p0, lower=lb, upper=ub)
            longer_tau = longer_autocorrelation_fit.param[1]
            longer_beta = longer_autocorrelation_fit.param[2]
            println("t=10^5 Tau = ", longer_tau)
            println("t=10^5 Beta = ", longer_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_tau, longer_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^5 Mean Squared Error = ", mse)

            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_beta)


        else
            # Now try fitting on log of data with offset already removed
            # TODO CHANGED THIS
            reduced_data = log.((1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- c))

            autocorrelation_fit = curve_fit(reduced_model, t_values, reduced_data, p0, lower=lb, upper=ub)

            reduced_tau = autocorrelation_fit.param[1]
            reduced_beta = autocorrelation_fit.param[2]
            println("t=10^4 Tau From Logged Data = ", reduced_tau)
            println("t=10^4 Beta from Logged Data = ", reduced_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [reduced_tau, reduced_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^4 Mean Squared Error = ", mse)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, reduced_tau)
            beta_by_temperature = append!(beta_by_temperature, reduced_beta)



            # Extract the parameters for t=10^5 fit
            # TODO CHANGED THIS
            longer_reduced_data = log.((1/(1-c)) * (autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- c))

            longer_autocorrelation_fit = curve_fit(reduced_model, longer_t_values, longer_reduced_data, p0, lower=lb, upper=ub)
            longer_reduced_tau = longer_autocorrelation_fit.param[1]
            longer_reduced_beta = longer_autocorrelation_fit.param[2]
            println("t=10^5 Tau From Logged Data = ", longer_reduced_tau)
            println("t=10^5 Beta from Logged Data = ", longer_reduced_beta)

            # TODO CHANGED THIS
            fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_reduced_tau, longer_reduced_beta]) .- c)
            residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
            mse = mean(residuals.^2)
            println("t=10^5 Mean Squared Error = ", mse)
            
            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_reduced_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_reduced_beta)
        end


        # Plot the fit
        dashed_cutoff = beta_relaxation_end
        extrapolated_dashed_cutoff = 5e4

        if temperature in fitted_temperatures
            # label="T = "*string(temperature)*", τ = "*string(tau_beta)*", β = "*string(beta)
            # TODO CHANGED ALL THESE
            plot!(long_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")
            plot!(short_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(extrapolated_dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")

            if temperature in minimal_graph_main_temperatures
                plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")
            end
        end

        # Plot the longer fit
        longer_dashed_cutoff = longer_beta_relaxation_end
        if temperature in fitted_temperatures
            # TODO CHANGED ALL THESE
            plot!(long_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="")
            plot!(short_fit_example_graph, (1/(1-c))*(autocorrelation_model(1:Int(extrapolated_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="")
        
            if temperature in minimal_graph_main_temperatures
                plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="")
            end
        end

    end

    # --- SAVE GRAPHS ---
    savefig(long_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.png"))
    savefig(long_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.pdf"))
    display(long_fit_example_graph)

    savefig(short_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.png"))
    savefig(short_fit_example_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.pdf"))
    display(short_fit_example_graph)

    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.png"))
    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.pdf"))
    display(minimal_graph)





    



        # -- RELAXATION TIME FITS ---
        # Make plot with tau_beta_by_temperature against temperature 
        # and then Arrhenius/VFT/parabolic fit curves to data too and output measure
        # Do fits on log of data to avoid numerical errors
        log_tau_beta_by_temperature = log10.(tau_beta_by_temperature) # Note using log10 here
        longer_log_tau_beta_by_temperature = log10.(longer_tau_beta_by_temperature)
        
        combined_fitted_temperatures = vcat(fitted_temperatures, fitted_temperatures+1e-15*randn(length(fitted_temperatures)))
        combined_log_tau_beta_by_temperature = vcat(log_tau_beta_by_temperature, longer_log_tau_beta_by_temperature)
        # Sort combined_fitted_temperatures but sort combined_log_tau_beta_by_temperature in the same way
        combined_fitted_temperatures, combined_log_tau_beta_by_temperature = combined_fitted_temperatures[sortperm(combined_fitted_temperatures)], combined_log_tau_beta_by_temperature[sortperm(combined_fitted_temperatures)]

        fit_temperatures = collect(LinRange(minimum(fitted_temperatures)-0.03, maximum(fitted_temperatures), 100))

        tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))
        squares_only_tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))
        circles_only_tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))

        scatter!(tau_fits_graph, fitted_temperatures, log_tau_beta_by_temperature, 
                 yerror=tau_errors_by_temperature, color=alex_red, label="", markersize=5, 
                 ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], alpha=0.8)
        scatter!(circles_only_tau_fits_graph, fitted_temperatures, log_tau_beta_by_temperature, 
                 yerror=tau_errors_by_temperature, color=alex_red, label="", markersize=5, marker=:circle, 
                 alpha=0.8, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5])

        scatter!(tau_fits_graph, fitted_temperatures, longer_log_tau_beta_by_temperature, 
                 yerror=longer_tau_errors_by_temperature, color=alex_red, label="", markersize=4, 
                 marker=:square, alpha=0.8)
        scatter!(squares_only_tau_fits_graph, fitted_temperatures, longer_log_tau_beta_by_temperature, 
                 yerror=longer_tau_errors_by_temperature, color=alex_red, label="", markersize=4, 
                 marker=:square, alpha=0.8, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5])


        # # Arrhenius fit
        # arrhenius_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T), all_temperatures, log_tau_beta_by_temperature, [0.1, 0.1], lower=[-Inf,-Inf], upper=[Inf, Inf])
        # arrhenius_fit_parameters = arrhenius_fit.param
        # arrhenius_fit_curve = arrhenius_fit_parameters[1] .+ (arrhenius_fit_parameters[2] ./ fit_temperatures)
        # lowest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][1] .+ (confint(arrhenius_fit, level=conf_int)[2][1] ./ fit_temperatures)
        # highest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][2] .+ (confint(arrhenius_fit, level=conf_int)[2][2] ./ fit_temperatures)
        # plot!(tau_fits_graph, fit_temperatures, arrhenius_fit_curve, label="Arrhenius Fit", color=alex_orange, ribbon=(arrhenius_fit_curve .- lowest_confidence_interval_arrhenius_curve, highest_confidence_interval_arrhenius_curve .- arrhenius_fit_curve), fc=:orange, fa=0.05)

        # println("Arrhenius Parameters: ", arrhenius_fit_parameters)
        # println("Arrhenius Residuals", arrhenius_fit.resid)
        # println("Arrhenius Covariance Matrix", vcov(arrhenius_fit))
        # println("Arrhenius Standard Errors", stderror(arrhenius_fit))
        # println("Arrhenius Margin Error: ", margin_error(arrhenius_fit))
        # println("Confdience Intervals", confint(arrhenius_fit))
        # println("")
        # println("")



        # Parabolic fit over all data
        parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        parabolic_fit_parameters = parabolic_fit.param
        parabolic_fit_curve = parabolic_fit_parameters[1] .+ (parabolic_fit_parameters[2] ./ fit_temperatures) .+ (parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(circles_only_tau_fits_graph, fit_temperatures, parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("")
        println("")
        println("Shorter Parabolic Parameters: ", parabolic_fit_parameters)
        # println("Shorter Parabolic Residuals", parabolic_fit.resid)
        # println("Shorter Parabolic Covariance Matrix", vcov(parabolic_fit))
        # println("Shorter Parabolic Standard Errors", stderror(parabolic_fit))
        # println("Shorter Parabolic Margin Error: ", margin_error(parabolic_fit))
        # println("Shorter Confdience Intervals", confint(parabolic_fit))

        longer_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, longer_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        longer_parabolic_fit_parameters = longer_parabolic_fit.param
        longer_parabolic_fit_curve = longer_parabolic_fit_parameters[1] .+ (longer_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(squares_only_tau_fits_graph, fit_temperatures, longer_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("Longer Parabolic Parameters: ", longer_parabolic_fit_parameters)
        # println("Longer Parabolic Residuals", longer_parabolic_fit.resid)
        # println("Longer Parabolic Covariance Matrix", vcov(longer_parabolic_fit))
        # println("Longer Parabolic Standard Errors", stderror(longer_parabolic_fit))
        # println("Longer Parabolic Margin Error: ", margin_error(longer_parabolic_fit))
        # println("Longer Confdience Intervals", confint(longer_parabolic_fit))

        combined_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        combined_parabolic_fit_parameters = combined_parabolic_fit.param
        combined_parabolic_fit_curve = combined_parabolic_fit_parameters[1] .+ (combined_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (combined_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        plot!(tau_fits_graph, fit_temperatures, combined_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("Combined Parabolic Parameters: ", combined_parabolic_fit_parameters)
        # println("Combined Parabolic Residuals", combined_parabolic_fit.resid)
        # println("Combined Parabolic Covariance Matrix", vcov(combined_parabolic_fit))
        # println("Combined Parabolic Standard Errors", stderror(combined_parabolic_fit))
        # println("Combined Parabolic Margin Error: ", margin_error(combined_parabolic_fit))
        # println("Combined Confdience Intervals", confint(combined_parabolic_fit))



        # VFT fit
        vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        vft_fit_parameters = vft_fit.param
        vft_fit_curve = vft_fit_parameters[1] .+ (vft_fit_parameters[2] ./ (fit_temperatures .- vft_fit_parameters[3]))
        plot!(circles_only_tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        println("")
        println("")
        println("VFT Parameters: ", vft_fit_parameters)
        # println("VFT Residuals", vft_fit.resid)
        # println("VFT Covariance Matrix", vcov(vft_fit))
        # println("VFT Standard Errors", stderror(vft_fit))
        # println("VFT Margin Error: ", margin_error(vft_fit))
        # println("Confdience Intervals", confint(vft_fit))


        longer_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, longer_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        longer_vft_fit_parameters = longer_vft_fit.param
        longer_vft_fit_curve = longer_vft_fit_parameters[1] .+ (longer_vft_fit_parameters[2] ./ (fit_temperatures .- longer_vft_fit_parameters[3]))
        plot!(squares_only_tau_fits_graph, fit_temperatures, longer_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        print("Longer VFT Parameters: ", longer_vft_fit_parameters)
        # println("Longer VFT Residuals", longer_vft_fit.resid)
        # println("Longer VFT Covariance Matrix", vcov(longer_vft_fit))
        # println("Longer VFT Standard Errors", stderror(longer_vft_fit))
        # println("Longer VFT Margin Error: ", margin_error(longer_vft_fit))
        # println("Longer Confdience Intervals", confint(longer_vft_fit))

        combined_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        combined_vft_fit_parameters = combined_vft_fit.param
        combined_vft_fit_curve = combined_vft_fit_parameters[1] .+ (combined_vft_fit_parameters[2] ./ (fit_temperatures .- combined_vft_fit_parameters[3]))
        plot!(tau_fits_graph, fit_temperatures, combined_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        print("Combined VFT Parameters: ", combined_vft_fit_parameters)
        # println("Combined VFT Residuals", combined_vft_fit.resid)
        # println("Combined VFT Covariance Matrix", vcov(combined_vft_fit))
        # println("Combined VFT Standard Errors", stderror(combined_vft_fit))
        # println("Combined VFT Margin Error: ", margin_error(combined_vft_fit))
        # println("Combined Confdience Intervals", confint(combined_vft_fit))


        #--- FRAGILITY INSET DATA ---
        scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], log_tau_beta_by_temperature,
        yerror=tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2,
        xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:circle, alpha=0.8)
        scatter!(circles_only_tau_fits_graph, [1/T for T in fitted_temperatures], log_tau_beta_by_temperature,
        yerror=tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2,
        xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:circle, alpha=0.8)

        scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], longer_log_tau_beta_by_temperature,
        yerror=longer_tau_errors_by_temperature, color=alex_red, legend=false, subplot=2, xlabel=L"1/T",
        ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:square, alpha=0.8)
        scatter!(squares_only_tau_fits_graph, [1/T for T in fitted_temperatures], longer_log_tau_beta_by_temperature,
        yerror=longer_tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2, xlabel=L"1/T",
        ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:square, alpha=0.8)


        # Plot fits in inset too
        # plot!(tau_fits_graph, [1/T for T in fit_temperatures], arrhenius_fit_curve,color=alex_orange, label="Arrhenius Fit", linestyle=:dash, subplot=2)
        
        plot!(circles_only_tau_fits_graph, [1/T for T in fit_temperatures], parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(circles_only_tau_fits_graph, [1/T for T in fit_temperatures], vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)

        plot!(squares_only_tau_fits_graph, [1/T for T in fit_temperatures], longer_parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(squares_only_tau_fits_graph, [1/T for T in fit_temperatures], longer_vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


        # Combined graph
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


        

        # --- BETA INSET DATA ---
        # Plot longer beta data
        scatter!(tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
        yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T", 
        ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])
        scatter!(squares_only_tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
        yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T",
        ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])

        # Plot the shorter beta data
        scatter!(tau_fits_graph, fitted_temperatures, beta_by_temperature,
        yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  subplot=3)
        scatter!(circles_only_tau_fits_graph, fitted_temperatures, beta_by_temperature, ylabel=L"\beta",
        yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])

        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.png"))
        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.pdf"))
        display(tau_fits_graph)

        savefig(circles_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_circles_only.png"))
        savefig(circles_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_circles_only.pdf"))
        display(circles_only_tau_fits_graph)

        savefig(squares_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_squares_only.png"))
        savefig(squares_only_tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits_squares_only.pdf"))
        display(squares_only_tau_fits_graph)

end














function other_L_autocorrelation_function_figures(L_values::Vector{Int})

    T_on_spline_dict = Dict()

    for L in L_values

        simulation_name = "combined_L=$(L)_slice"

        ### --- READ IN DATA ---

        ## -- Define the Data --

        # Read this in for all the following filenames for different temperatures
        minimal_graph_main_temperatures = [0.8,0.9,0.92,0.95,1.0,1.1,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0,7.5,10.0]

        if L != 5
            minimal_graph_main_temperatures = vcat(minimal_graph_main_temperatures, [1.75,2.25])
        end

        filenames = ["L_$(L)_T_$(T)_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv" for T in minimal_graph_main_temperatures]

        # filenames = [
        #     "L_5_T_0.8_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_0.9_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_0.92_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_0.95_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_1.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_1.1_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_1.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_2.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_2.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        #     "L_5_T_3.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
        # ] 




        ## -- Read in the data --

        autocorrelation_functions_by_temperature = Dict()

        # For each file in filenames, read in the data as a matrix
        for (i, filename) in pairs(filenames)
            data = readdlm(joinpath("results/autocorrelation_anneal_results/final_data", filename), ',', skipstart=3)

            temperature = data[1]
            samples_in_average = data[2]
            autocorrelation_function = data[3:end]
            
            # Store the data in the dictionary
            autocorrelation_functions_by_temperature[temperature] = autocorrelation_function
        end

        all_temperatures = sort(collect(keys(autocorrelation_functions_by_temperature)))


        ### --- COLOURS ---
        Plots.default(dpi = 600)

        alex_red = RGB(227/255, 11/255, 92/255)
        alex_orange = RGB(255/255, 165/255, 0/255)
        alex_pink = RGB(255/255, 105/255, 180/255)
        alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
        alex_blue = RGB(100/255, 149/255, 237/255)
        alex_alt_blue = RGB(4/255, 57/255, 94/255)

        # Define a fixed table of temperatures and their corresponding colors
        temp_to_color = Dict()
        alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]

        for temperature in all_temperatures
                temp_to_color[temperature] = alex_colors[mod1(findall(x -> x == temperature, all_temperatures)[1], length(alex_colors))]
        end
        temp_to_color[0.6] = alex_red
        temp_to_color[0.9] = alex_pink
        temp_to_color[0.92] = alex_blue
        temp_to_color[0.95] = alex_alt_blue
        temp_to_color[1.0] = alex_orange
        temp_to_color[1.2] = alex_green











        ### --- PLOTTING ---
        
        ## -- COMBINED GRAPHS ---
        minimal_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\bar\mathcal{C}(t)", legend=(0.9,0.35), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0], xticks=[0,2.5e4,5.0e4,7.5e4,1.0e5,1.2e5])

        c = 0.201388888888888
        long_window_length = 99900 # Int(1.2e5)
        short_window_length = Int(5e4)


        ## -- MINIMAL GRAPH ---
        main_minimal_graph_window_length = Int(0.999e5)

        for temperature in minimal_graph_main_temperatures     
            # TODO CHANGED THIS
            plot!(minimal_graph, 0:main_minimal_graph_window_length-1, (1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:main_minimal_graph_window_length] .- c), label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
        end


        


        ### --- FITTING PARAMETER EXTRACTION ---
        excluded_temperatures_from_fitting = []
        fitted_temperatures = [T for T in all_temperatures if T ∉ excluded_temperatures_from_fitting]

        tau_beta_by_temperature = Float64[]
        beta_by_temperature = Float64[]
        tau_errors_by_temperature = Float64[]  # Add this line
        beta_errors_by_temperature = Float64[]  # Add beta error array

        longer_tau_beta_by_temperature = Float64[]
        longer_beta_by_temperature = Float64[]
        longer_tau_errors_by_temperature = Float64[]  # Add this line
        longer_beta_errors_by_temperature = Float64[]  # Add longer beta error array

        for temperature in fitted_temperatures
            if temperature in excluded_temperatures_from_fitting
                continue
            end
            println("Temperature: ", temperature)


            # Fit the data to the autocorrelation function
            p0 = [10, 0.9]
            lb = [1e-3, 0.1] # example lower bounds
            ub = [1e11, 2.0]  # example upper bounds

            beta_relaxation_end = 1e4
            longer_beta_relaxation_end = 99900

            t_values = 1:beta_relaxation_end
            longer_t_values = 1:longer_beta_relaxation_end

            # Define the model function
            c = 0.201388888888888
            # TODO DIDN'T CHANGE THIS BUT NOW THIS IS CORRECT
            autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

            # TODO IS THIS NEEDED NOW WE'VE DONE CHANGES?
            if temperature >= 0.4

                # Fit the model to the data
                autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

                # Extract the parameters for t=10^4 fit
                tau = autocorrelation_fit.param[1]
                beta = autocorrelation_fit.param[2]
                
                # Extract standard errors
                param_errors = stderror(autocorrelation_fit)
                tau_error = param_errors[1]
                # Convert tau error to log10(tau) error: d(log10(tau))/dtau = 1/(tau*ln(10))
                log_tau_error = tau_error / (tau * log(10))
                
                println("t=10^4 Tau = ", tau)
                println("t=10^4 Beta = ", beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [tau, beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^4 Mean Squared Error = ", mse)

                tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
                beta_by_temperature = append!(beta_by_temperature, beta)
                tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
                beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



                # Extract the parameters for t=10^5 fit
                longer_autocorrelation_fit = curve_fit(autocorrelation_model, longer_t_values, autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)], p0, lower=lb, upper=ub)
                longer_tau = longer_autocorrelation_fit.param[1]
                longer_beta = longer_autocorrelation_fit.param[2]
                
                # Extract standard errors for longer fit
                longer_param_errors = stderror(longer_autocorrelation_fit)
                longer_tau_error = longer_param_errors[1]
                # Convert tau error to log10(tau) error
                longer_log_tau_error = longer_tau_error / (longer_tau * log(10))


                println("t=10^5 Tau = ", longer_tau)
                println("t=10^5 Beta = ", longer_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_tau, longer_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^5 Mean Squared Error = ", mse)

                longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_tau)
                longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_beta)
                longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
                longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])


            else
                # Now try fitting on log of data with offset already removed
                # TODO CHANGED THIS
                reduced_data = log.((1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- c))


                autocorrelation_fit = curve_fit(reduced_model, t_values, reduced_data, p0, lower=lb, upper=ub)

                reduced_tau = autocorrelation_fit.param[1]
                reduced_beta = autocorrelation_fit.param[2]
                println("t=10^4 Tau From Logged Data = ", reduced_tau)
                println("t=10^4 Beta from Logged Data = ", reduced_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [reduced_tau, reduced_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^4 Mean Squared Error = ", mse)

                tau_beta_by_temperature = append!(tau_beta_by_temperature, reduced_tau)
                beta_by_temperature = append!(beta_by_temperature, reduced_beta)
                tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
                beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



                # Extract the parameters for t=10^5 fit
                # TODO CHANGED THIS
                longer_reduced_data = log.((1/(1-c)) * (autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- c))

                longer_autocorrelation_fit = curve_fit(reduced_model, longer_t_values, longer_reduced_data, p0, lower=lb, upper=ub)
                longer_reduced_tau = longer_autocorrelation_fit.param[1]
                longer_reduced_beta = longer_autocorrelation_fit.param[2]
                println("t=10^5 Tau From Logged Data = ", longer_reduced_tau)
                println("t=10^5 Beta from Logged Data = ", longer_reduced_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_reduced_tau, longer_reduced_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^5 Mean Squared Error = ", mse)
                
                longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_reduced_tau)
                longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_reduced_beta)
                longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
                longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])
            end


            # Plot the fit
            dashed_cutoff = beta_relaxation_end
            extrapolated_dashed_cutoff = 5e4

            if temperature in fitted_temperatures
                # label="T = "*string(temperature)*", τ = "*string(tau_beta)*", β = "*string(beta)
                if temperature in minimal_graph_main_temperatures
                    plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c), color=:black, linestyle=:dash, label="")
                end
            end

            # Plot the longer fit
            longer_dashed_cutoff = longer_beta_relaxation_end
            if temperature in fitted_temperatures
                if temperature in minimal_graph_main_temperatures
                    plot!(minimal_graph, (1/(1-c))*(autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c), color=:green, linestyle=:dash, label="")
                end
            end

        end

        # --- SAVE GRAPHS ---
        # savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.png"))
        # savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.pdf"))
        # display(minimal_graph)





        



            # -- RELAXATION TIME FITS ---
            # Make plot with tau_beta_by_temperature against temperature 
            # and then Arrhenius/VFT/parabolic fit curves to data too and output measure
            # Do fits on log of data to avoid numerical errors
            log_tau_beta_by_temperature = log10.(tau_beta_by_temperature) # Note using log10 here
            longer_log_tau_beta_by_temperature = log10.(longer_tau_beta_by_temperature)
            
            combined_fitted_temperatures = vcat(fitted_temperatures, fitted_temperatures+1e-15*randn(length(fitted_temperatures)))
            combined_log_tau_beta_by_temperature = vcat(log_tau_beta_by_temperature, longer_log_tau_beta_by_temperature)
            # Sort combined_fitted_temperatures but sort combined_log_tau_beta_by_temperature in the same way
            combined_fitted_temperatures, combined_log_tau_beta_by_temperature = combined_fitted_temperatures[sortperm(combined_fitted_temperatures)], combined_log_tau_beta_by_temperature[sortperm(combined_fitted_temperatures)]

            fit_temperatures = collect(LinRange(minimum(fitted_temperatures)-0.03, maximum(fitted_temperatures), 100))

            tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))

            scatter!(tau_fits_graph, fitted_temperatures, log_tau_beta_by_temperature, 
                     yerror=tau_errors_by_temperature, color=alex_red, label="", markersize=5, 
                     ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], alpha=0.8)

            scatter!(tau_fits_graph, fitted_temperatures, longer_log_tau_beta_by_temperature, 
                     yerror=longer_tau_errors_by_temperature, color=alex_red, label="", markersize=4, 
                     marker=:square, alpha=0.8)
    

            # # Arrhenius fit
            # arrhenius_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T), all_temperatures, log_tau_beta_by_temperature, [0.1, 0.1], lower=[-Inf,-Inf], upper=[Inf, Inf])
            # arrhenius_fit_parameters = arrhenius_fit.param
            # arrhenius_fit_curve = arrhenius_fit_parameters[1] .+ (arrhenius_fit_parameters[2] ./ fit_temperatures)
            # lowest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][1] .+ (confint(arrhenius_fit, level=conf_int)[2][1] ./ fit_temperatures)
            # highest_confidence_interval_arrhenius_curve = confint(arrhenius_fit, level=conf_int)[1][2] .+ (confint(arrhenius_fit, level=conf_int)[2][2] ./ fit_temperatures)
            # plot!(tau_fits_graph, fit_temperatures, arrhenius_fit_curve, label="Arrhenius Fit", color=alex_orange, ribbon=(arrhenius_fit_curve .- lowest_confidence_interval_arrhenius_curve, highest_confidence_interval_arrhenius_curve .- arrhenius_fit_curve), fc=:orange, fa=0.05)

            # println("Arrhenius Parameters: ", arrhenius_fit_parameters)
            # println("Arrhenius Residuals", arrhenius_fit.resid)
            # println("Arrhenius Covariance Matrix", vcov(arrhenius_fit))
            # println("Arrhenius Standard Errors", stderror(arrhenius_fit))
            # println("Arrhenius Margin Error: ", margin_error(arrhenius_fit))
            # println("Confdience Intervals", confint(arrhenius_fit))
            # println("")
            # println("")



            # Parabolic fit over all data
            parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
            parabolic_fit_parameters = parabolic_fit.param
            parabolic_fit_curve = parabolic_fit_parameters[1] .+ (parabolic_fit_parameters[2] ./ fit_temperatures) .+ (parabolic_fit_parameters[3] ./ fit_temperatures.^2)

            println("")
            println("")
            println("Shorter Parabolic Parameters: ", parabolic_fit_parameters)
            # println("Shorter Parabolic Residuals", parabolic_fit.resid)
            # println("Shorter Parabolic Covariance Matrix", vcov(parabolic_fit))
            # println("Shorter Parabolic Standard Errors", stderror(parabolic_fit))
            # println("Shorter Parabolic Margin Error: ", margin_error(parabolic_fit))
            # println("Shorter Confdience Intervals", confint(parabolic_fit))

            longer_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), fitted_temperatures, longer_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
            longer_parabolic_fit_parameters = longer_parabolic_fit.param
            longer_parabolic_fit_curve = longer_parabolic_fit_parameters[1] .+ (longer_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_parabolic_fit_parameters[3] ./ fit_temperatures.^2)

            println("Longer Parabolic Parameters: ", longer_parabolic_fit_parameters)
            # println("Longer Parabolic Residuals", longer_parabolic_fit.resid)
            # println("Longer Parabolic Covariance Matrix", vcov(longer_parabolic_fit))
            # println("Longer Parabolic Standard Errors", stderror(longer_parabolic_fit))
            # println("Longer Parabolic Margin Error: ", margin_error(longer_parabolic_fit))
            # println("Longer Confdience Intervals", confint(longer_parabolic_fit))

            combined_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
            combined_parabolic_fit_parameters = combined_parabolic_fit.param
            combined_parabolic_fit_curve = combined_parabolic_fit_parameters[1] .+ (combined_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (combined_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
            plot!(tau_fits_graph, fit_temperatures, combined_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

            println("Combined Parabolic Parameters: ", combined_parabolic_fit_parameters)
            # println("Combined Parabolic Residuals", combined_parabolic_fit.resid)
            # println("Combined Parabolic Covariance Matrix", vcov(combined_parabolic_fit))
            # println("Combined Parabolic Standard Errors", stderror(combined_parabolic_fit))
            # println("Combined Parabolic Margin Error: ", margin_error(combined_parabolic_fit))
            # println("Combined Confdience Intervals", confint(combined_parabolic_fit))



            # VFT fit
            vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
            vft_fit_parameters = vft_fit.param
            vft_fit_curve = vft_fit_parameters[1] .+ (vft_fit_parameters[2] ./ (fit_temperatures .- vft_fit_parameters[3]))

            println("")
            println("")
            println("VFT Parameters: ", vft_fit_parameters)
            # println("VFT Residuals", vft_fit.resid)
            # println("VFT Covariance Matrix", vcov(vft_fit))
            # println("VFT Standard Errors", stderror(vft_fit))
            # println("VFT Margin Error: ", margin_error(vft_fit))
            # println("Confdience Intervals", confint(vft_fit))


            longer_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), fitted_temperatures, longer_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
            longer_vft_fit_parameters = longer_vft_fit.param
            longer_vft_fit_curve = longer_vft_fit_parameters[1] .+ (longer_vft_fit_parameters[2] ./ (fit_temperatures .- longer_vft_fit_parameters[3]))

            print("Longer VFT Parameters: ", longer_vft_fit_parameters)
            # println("Longer VFT Residuals", longer_vft_fit.resid)
            # println("Longer VFT Covariance Matrix", vcov(longer_vft_fit))
            # println("Longer VFT Standard Errors", stderror(longer_vft_fit))
            # println("Longer VFT Margin Error: ", margin_error(longer_vft_fit))
            # println("Longer Confdience Intervals", confint(longer_vft_fit))

            combined_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), combined_fitted_temperatures, combined_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
            combined_vft_fit_parameters = combined_vft_fit.param
            combined_vft_fit_curve = combined_vft_fit_parameters[1] .+ (combined_vft_fit_parameters[2] ./ (fit_temperatures .- combined_vft_fit_parameters[3]))
            plot!(tau_fits_graph, fit_temperatures, combined_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

            print("Combined VFT Parameters: ", combined_vft_fit_parameters)
            # println("Combined VFT Residuals", combined_vft_fit.resid)
            # println("Combined VFT Covariance Matrix", vcov(combined_vft_fit))
            # println("Combined VFT Standard Errors", stderror(combined_vft_fit))
            # println("Combined VFT Margin Error: ", margin_error(combined_vft_fit))
            # println("Combined Confdience Intervals", confint(combined_vft_fit))


            #--- FRAGILITY INSET DATA ---
            # TODO readd
            scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], log_tau_beta_by_temperature,
            yerror=tau_errors_by_temperature, color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35),  subplot=2, 
            # xticks=(0.25,0.75,1.25,1.75),
            xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:circle, alpha=0.8)

            scatter!(tau_fits_graph, [1/T for T in fitted_temperatures], longer_log_tau_beta_by_temperature,
            yerror=longer_tau_errors_by_temperature, color=alex_red, legend=false, subplot=2, xlabel=L"1/T",
            ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(longer_log_tau_beta_by_temperature) .+ [-0.5, 0.5], marker=:square, alpha=0.8)


            # Plot fits in inset too
            # # plot!(tau_fits_graph, [1/T for T in fit_temperatures], arrhenius_fit_curve,color=alex_orange, label="Arrhenius Fit", linestyle=:dash, subplot=2)
            

            # Combined graph
            # TODO readd
            plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
            plot!(tau_fits_graph, [1/T for T in fit_temperatures], combined_vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


            

            # --- BETA INSET DATA ---
            # Plot longer beta data
            scatter!(tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
            yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T", 
            ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])
            scatter!(squares_only_tau_fits_graph, fitted_temperatures, longer_beta_by_temperature,
            yerror=longer_beta_errors_by_temperature, color=alex_pink, marker=:square, legend=false, inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlabel=L"T",
            ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])

            # Plot the shorter beta data
            scatter!(tau_fits_graph, fitted_temperatures, beta_by_temperature,
            yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  subplot=3)
            scatter!(circles_only_tau_fits_graph, fitted_temperatures, beta_by_temperature, ylabel=L"\beta",
            yerror=beta_errors_by_temperature, color=alex_pink, marker=:circle, legend=false,  inset=bbox(0.27,0.05,0.3,0.35), subplot=3, xlims=(0,10.5), ylims=(0,1.05), xticks=[0,2,4,6,8,10], yticks=[0,0.2,0.4,0.6,0.8,1.0])










            # STRECHED EXPONENTIAL DEFINITION STUFF --
            # Create spline for all beta data against temperature
            beta_by_temperature = longer_beta_by_temperature
            temperatures_for_beta = fitted_temperatures
            temperature_from_beta_spline = CubicSpline(beta_by_temperature[sortperm(beta_by_temperature)], temperatures_for_beta[sortperm(beta_by_temperature)], extrapl=[1,], extrapr=[1,])
            beta_from_temperature_spline = CubicSpline(temperatures_for_beta[sortperm(temperatures_for_beta)], beta_by_temperature[sortperm(temperatures_for_beta)], extrapl=[1,], extrapr=[1,])
            # Plot spline
            plot!(tau_fits_graph, temperature_from_beta_spline(beta_by_temperature), beta_by_temperature, color=alex_orange, label="Spline", subplot=3)

            # Onset of stretching definition
            # beta_stretch = 0.8
            # Define beta_stretch as 0.9*beta(T=10.0)
            # beta_stretch = 0.8*beta_from_temperature_spline(10.0)
            beta_infinite_temperature_dict = Dict(5 => 1.1812820380521043, 7 => 1.0838099849875709, 11 => 1.0229843232478562, 9 => 1.0370097148780135)
            beta_stretch = 0.75*beta_infinite_temperature_dict[L]

            # Find temperature where beta = beta_stretch
            temperature_at_beta_stretch = temperature_from_beta_spline(beta_stretch)
            T_on_spline_dict[L] = temperature_at_beta_stretch
            println("L=$(L), Temperature at β = $beta_stretch: ", temperature_at_beta_stretch)
            # Plot temperature_at_beta_stretch
            scatter!(tau_fits_graph, [temperature_at_beta_stretch], [beta_stretch], color=alex_red, label="T at β = 0.$(temperature_at_beta_stretch)", marker=:circle, subplot=3, markersize=4)


            # Plot hline at \beta=beta_stregch
            hline!(tau_fits_graph, [beta_stretch], color=alex_red, label="β = $(beta_stretch)", linestyle=:dash, linewidth=2, subplot=3)
            vline!(tau_fits_graph, [temperature_at_beta_stretch], color=alex_red, label="T = $(temperature_at_beta_stretch)", linestyle=:dash, linewidth=2, subplot=3)
            annotate!(tau_fits_graph, [(temperature_at_beta_stretch+0.2, 0.3, text(L"\bar T^{\rm on}", 8, :left, color=alex_red))], subplot=3)


            savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.png"))
            savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.pdf"))
            display(tau_fits_graph)

    end

    println("T_on_dict: $(T_on_spline_dict)")

end























function combined_L_beta_plot()

    L_values = [5,7,9,11]


    combined_L_beta_graph = plot([], [], xlabel=L"T", ylabel=L"\beta", label="", yguidefontsize=12,xguidefontsize=12, linewidth=3)
    
    combined_L_beta_rescaled_graph = plot([], [], xlabel=L"T", ylabel=L"\beta/\beta(T=\infty)", label="", yguidefontsize=12,xguidefontsize=12, linewidth=3)


            ### --- COLOURS ---
            Plots.default(dpi = 600)

            alex_red = RGB(227/255, 11/255, 92/255)
            alex_orange = RGB(255/255, 165/255, 0/255)
            alex_pink = RGB(255/255, 105/255, 180/255)
            alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
            alex_blue = RGB(100/255, 149/255, 237/255)
            alex_alt_blue = RGB(4/255, 57/255, 94/255)
    
            # Define a fixed table of temperatures and their corresponding colors
            temp_to_color = Dict()
            alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]
    
    


    for L in L_values

        if L != 11
            simulation_name = "combined_L=$(L)_slice"

            # Read this in for all the following filenames for different temperatures
            minimal_graph_main_temperatures = [0.8,0.9,0.92,0.95,1.0,1.1,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0,7.5,10.0]

            if L != 5
                minimal_graph_main_temperatures = vcat(minimal_graph_main_temperatures, [1.75,2.25])
            end

            filenames = ["L_$(L)_T_$(T)_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv" for T in minimal_graph_main_temperatures]


        else

            minimal_graph_main_temperatures = [0.9,0.92,0.95,0.97,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.4,1.5,1.75,2.0,2.25,3.0,4.0,5.0,7.5,10.0]

            filenames = [
                # "L_11_T_0.6_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
                # "L_11_T_0.8_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", #
                "L_11_T_0.9_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
                # "L_11_T_0.91_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_0.92_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv", #
                # "L_11_T_0.93_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", #
                "L_11_T_0.95_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_0.97_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.05_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.1_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.15_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.2_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.25_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.3_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.4_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_1.5_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
                "L_11_T_1.75_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", 
                "L_11_T_2.0_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
                "L_11_T_2.25_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", 
                "L_11_T_3.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_4.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_5.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_7.5_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
                "L_11_T_10.0_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv"] 

        end




        ## -- Read in the data --

        autocorrelation_functions_by_temperature = Dict()

        # For each file in filenames, read in the data as a matrix
        for (i, filename) in pairs(filenames)
            data = readdlm(joinpath("results/autocorrelation_anneal_results/final_data", filename), ',', skipstart=3)

            temperature = data[1]
            samples_in_average = data[2]
            autocorrelation_function = data[3:end]
            
            # Store the data in the dictionary
            autocorrelation_functions_by_temperature[temperature] = autocorrelation_function
        end

        all_temperatures = sort(collect(keys(autocorrelation_functions_by_temperature)))





        ### --- PLOTTING ---

        c = 0.201388888888888
        long_window_length = 99900 # Int(1.2e5)
        short_window_length = Int(5e4)

        fitted_temperatures = [T for T in all_temperatures]

        tau_beta_by_temperature = Float64[]
        beta_by_temperature = Float64[]
        tau_errors_by_temperature = Float64[]  # Add this line
        beta_errors_by_temperature = Float64[]  # Add beta error array

        longer_tau_beta_by_temperature = Float64[]
        longer_beta_by_temperature = Float64[]
        longer_tau_errors_by_temperature = Float64[]  # Add this line
        longer_beta_errors_by_temperature = Float64[]  # Add longer beta error array


        for temperature in fitted_temperatures
            println("Temperature: ", temperature)


            # Fit the data to the autocorrelation function
            p0 = [10, 0.9]
            lb = [1e-3, 0.1] # example lower bounds
            ub = [1e11, 2.0]  # example upper bounds

            beta_relaxation_end = 1e4
            longer_beta_relaxation_end = 99900

            t_values = 1:beta_relaxation_end
            longer_t_values = 1:longer_beta_relaxation_end

            # Define the model function
            c = 0.201388888888888
            # TODO DIDN'T CHANGE THIS BUT NOW THIS IS CORRECT
            autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

            # TODO IS THIS NEEDED NOW WE'VE DONE CHANGES?
            if temperature >= 0.4

                # Fit the model to the data
                autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

                # Extract the parameters for t=10^4 fit
                tau = autocorrelation_fit.param[1]
                beta = autocorrelation_fit.param[2]
                
                # Extract standard errors
                param_errors = stderror(autocorrelation_fit)
                tau_error = param_errors[1]
                # Convert tau error to log10(tau) error: d(log10(tau))/dtau = 1/(tau*ln(10))
                log_tau_error = tau_error / (tau * log(10))
                
                println("t=10^4 Tau = ", tau)
                println("t=10^4 Beta = ", beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [tau, beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^4 Mean Squared Error = ", mse)

                tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
                beta_by_temperature = append!(beta_by_temperature, beta)
                tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
                beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



                # Extract the parameters for t=10^5 fit
                longer_autocorrelation_fit = curve_fit(autocorrelation_model, longer_t_values, autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)], p0, lower=lb, upper=ub)
                longer_tau = longer_autocorrelation_fit.param[1]
                longer_beta = longer_autocorrelation_fit.param[2]
                
                # Extract standard errors for longer fit
                longer_param_errors = stderror(longer_autocorrelation_fit)
                longer_tau_error = longer_param_errors[1]
                # Convert tau error to log10(tau) error
                longer_log_tau_error = longer_tau_error / (longer_tau * log(10))


                println("t=10^5 Tau = ", longer_tau)
                println("t=10^5 Beta = ", longer_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_tau, longer_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^5 Mean Squared Error = ", mse)

                longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_tau)
                longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_beta)
                longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
                longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])


            else
                # Now try fitting on log of data with offset already removed
                # TODO CHANGED THIS
                reduced_data = log.((1/(1-c))*(autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- c))


                autocorrelation_fit = curve_fit(reduced_model, t_values, reduced_data, p0, lower=lb, upper=ub)

                reduced_tau = autocorrelation_fit.param[1]
                reduced_beta = autocorrelation_fit.param[2]
                
                # Extract standard errors for reduced model
                param_errors = stderror(autocorrelation_fit)
                tau_error = param_errors[1]
                # For reduced model, tau error is already in appropriate scale
                log_tau_error = tau_error / (reduced_tau * log(10))
                
                println("t=10^4 Tau From Logged Data = ", reduced_tau)
                println("t=10^4 Beta from Logged Data = ", reduced_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(t_values, [reduced_tau, reduced_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^4 Mean Squared Error = ", mse)

                tau_beta_by_temperature = append!(tau_beta_by_temperature, reduced_tau)
                beta_by_temperature = append!(beta_by_temperature, reduced_beta)
                tau_errors_by_temperature = append!(tau_errors_by_temperature, log_tau_error)
                beta_errors_by_temperature = append!(beta_errors_by_temperature, stderror(autocorrelation_fit)[2])



                # Extract the parameters for t=10^5 fit
                # TODO CHANGED THIS
                longer_reduced_data = log.((1/(1-c)) * (autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- c))

                longer_autocorrelation_fit = curve_fit(reduced_model, longer_t_values, longer_reduced_data, p0, lower=lb, upper=ub)
                longer_reduced_tau = longer_autocorrelation_fit.param[1]
                longer_reduced_beta = longer_autocorrelation_fit.param[2]
                
                # Extract standard errors for longer reduced model
                longer_param_errors = stderror(longer_autocorrelation_fit)
                longer_tau_error = longer_param_errors[1]
                longer_log_tau_error = longer_tau_error / (longer_reduced_tau * log(10))

                println("t=10^5 Tau From Logged Data = ", longer_reduced_tau)
                println("t=10^5 Beta from Logged Data = ", longer_reduced_beta)

                # TODO CHANGED THIS
                fitted_values = (1/(1-c))*(autocorrelation_model(longer_t_values, [longer_reduced_tau, longer_reduced_beta]) .- c)
                residuals = autocorrelation_functions_by_temperature[temperature][1:Int(longer_beta_relaxation_end)] .- fitted_values
                mse = mean(residuals.^2)
                println("t=10^5 Mean Squared Error = ", mse)
                
                longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_reduced_tau)
                longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_reduced_beta)
                longer_tau_errors_by_temperature = append!(longer_tau_errors_by_temperature, longer_log_tau_error)
                longer_beta_errors_by_temperature = append!(longer_beta_errors_by_temperature, stderror(longer_autocorrelation_fit)[2])
            end

        end

        # Plot beta by temperature values for each L on the same graph
        scatter!(combined_L_beta_graph, fitted_temperatures, beta_by_temperature, color=alex_colors[Int(((L-5)/2)+1)], label="L = "*string(L), markersize=5, alpha=0.8)
        scatter!(combined_L_beta_graph, fitted_temperatures, longer_beta_by_temperature, color=alex_colors[Int(((L-5)/2)+1)], label="", markersize=4, marker=:square, alpha=0.8)


        # # Now plot a new graph where they are each rescaled by \Beta(T=10.0)
        # Define the rescaling factor
        # rescale_temperature = 10.0
        # beta_rescale = beta_by_temperature[findall(all_temperatures .== rescale_temperature)[1]]
        beta_infinite_temperature_dict = Dict(5 => 1.1812820380521043, 7 => 1.0838099849875709, 11 => 1.0229843232478562, 9 => 1.0370097148780135)
        beta_rescale = beta_infinite_temperature_dict[L]
        # Plot beta by temperature values for each L on the same graph
        scatter!(combined_L_beta_rescaled_graph, fitted_temperatures, beta_by_temperature./ beta_rescale, color=alex_colors[Int(((L-5)/2)+1)], label="L = "*string(L), markersize=5, alpha=0.8)
        scatter!(combined_L_beta_rescaled_graph, fitted_temperatures, longer_beta_by_temperature./ beta_rescale, color=alex_colors[Int(((L-5)/2)+1)], label="", markersize=4, marker=:square, alpha=0.8)


        # # Now instead of a rescaling factor do a vertical shift alignment to \Beta(T=10.0)
        # # Define the rescaling factor
        # rescale_temperature = 10.0
        # beta_rescale = beta_by_temperature[findall(all_temperatures .== rescale_temperature)[1]]
        # # Plot beta by temperature values for each L on the same graph
        # scatter!(combined_L_beta_rescaled_graph, fitted_temperatures, beta_by_temperature .- beta_rescale, color=alex_colors[Int(((L-5)/2)+1)], label="L = "*string(L), markersize=5, alpha=0.8)
        # scatter!(combined_L_beta_rescaled_graph, fitted_temperatures, longer_beta_by_temperature .- beta_rescale, color=alex_colors[Int(((L-5)/2)+1)], label="", markersize=4, marker=:square, alpha=0.8)

    end



    # Plot T_otimes as vertical lines
    # epsilon_otimes_dict = Dict(5 => -0.2137176988190482, 7 => -0.23134348655947637, 11 => -0.23927293015680562, 9 => -0.23676323523181977, 3 => -0.16330243097243677)
    # T_otimes_dict = Dict(5  => 2.54728, 7  => 2.13763, 11 => 2.06398,9  => 2.06398)
    # Also plot the points (T_otimes, beta(T_otimes)) on the graph
    # for (L, T_otimes) in T_otimes_dict
    #     vline!(combined_L_beta_graph, [T_otimes], color=alex_colors[Int(((L-5)/2)+1)], label="LT_{otimes} = "*string(T_otimes)*" (L=$(L))", linestyle=:dash, linewidth=2)
    #     vline!(combined_L_beta_rescaled_graph, [T_otimes], color=alex_colors[Int(((L-5)/2)+1)], label="LT_{otimes} = "*string(T_otimes)*" (L=$(L))", linestyle=:dash, linewidth=2)

    # end

    savefig(combined_L_beta_graph, joinpath("results/autocorrelation_anneal_results","combined_L_beta_plot.png"))
    savefig(combined_L_beta_graph, joinpath("results/autocorrelation_anneal_results","combined_L_beta_plot.pdf"))
    display(combined_L_beta_graph)

    savefig(combined_L_beta_rescaled_graph, joinpath("results/autocorrelation_anneal_results","combined_L_beta_rescaled_plot.png"))
    savefig(combined_L_beta_rescaled_graph, joinpath("results/autocorrelation_anneal_results","combined_L_beta_rescaled_plot.pdf"))
    display(combined_L_beta_rescaled_graph)

end




function plot_autocorrelation_with_fits(T::Float64, t_final::Int; L::Int = 5)

    # --- Read in the data ---
    filename = "L_$(L)_T_$(T)_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv"
    filepath = joinpath("results/autocorrelation_anneal_results/final_data", filename)
    data = readdlm(filepath, ',', skipstart=3)

    # Extract data
    temperature = data[1]
    samples_in_average = data[2]
    autocorrelation_function = data[3:end]

    # Truncate autocorrelation function to t_final
    t_final = min(t_final, length(autocorrelation_function))
    autocorrelation_function = autocorrelation_function[1:t_final]

    # Adjust the autocorrelation function
    c = 0.201388888888888
    adjusted_autocorrelation_function = (1 / (1 - c)) * (autocorrelation_function .- c)

    # Define time values
    t_values = 1:t_final

    # --- Plot the autocorrelation function ---
    plot(t_values, adjusted_autocorrelation_function,
        label="Data",
        xlabel="Time [MC Steps]",
        ylabel=L"\bar\mathcal{C}(t)",
        title="Autocorrelation Function at T = $(T)",
        legend=:topright)

    # --- Define the model function ---
    autocorrelation_model(t, p) = exp.(- (t ./ p[1]) .^ p[2])

    # Initial guesses and bounds
    p0 = [10.0, 0.9]
    lb = [1e-3, 0.1]
    ub = [1e11, 2.0]

    # --- Shorter fit ---
    beta_relaxation_end = min(1e4, t_final)
    t_fit_values_short = 1:Int(beta_relaxation_end)

    if T >= 0.4
        # Fit the model to the data
        autocorrelation_fit_short = curve_fit(
            autocorrelation_model,
            t_fit_values_short,
            adjusted_autocorrelation_function[1:Int(beta_relaxation_end)],
            p0,
            lower=lb,
            upper=ub)
    else
        # Fit on log of data
        reduced_data = log.(adjusted_autocorrelation_function[1:Int(beta_relaxation_end)])
        autocorrelation_fit_short = curve_fit(
            reduced_model,
            t_fit_values_short,
            reduced_data,
            p0,
            lower=lb,
            upper=ub)
    end

    # Extract parameters
    tau_short = autocorrelation_fit_short.param[1]
    beta_short = autocorrelation_fit_short.param[2]

    # Generate fitted values
    fitted_values_short = autocorrelation_model(t_values, [tau_short, beta_short])

    # Plot the shorter fit
    plot!(t_values, fitted_values_short,
        linestyle=:dash,
        color=:black,
        label="Shorter Fit (τ = $(round(tau_short; digits=2)), β = $(round(beta_short; digits=2)))")

    # --- Longer fit ---
    longer_beta_relaxation_end = min(99900, t_final)
    t_fit_values_long = 1:Int(longer_beta_relaxation_end)

    if T >= 0.4
        autocorrelation_fit_long = curve_fit(
            autocorrelation_model,
            t_fit_values_long,
            adjusted_autocorrelation_function[1:Int(longer_beta_relaxation_end)],
            p0,
            lower=lb,
            upper=ub)
    else
        longer_reduced_data = log.(adjusted_autocorrelation_function[1:Int(longer_beta_relaxation_end)])
        autocorrelation_fit_long = curve_fit(
            reduced_model,
            t_fit_values_long,
            longer_reduced_data,
            p0,
            lower=lb,
            upper=ub)
    end

    # Extract parameters
    tau_long = autocorrelation_fit_long.param[1]
    beta_long = autocorrelation_fit_long.param[2]

    # Generate fitted values
    fitted_values_long = autocorrelation_model(t_values, [tau_long, beta_long])

    # Plot the longer fit
    plot!(t_values, fitted_values_long,
        linestyle=:dash,
        color=:red,
        label="Longer Fit (τ = $(round(tau_long; digits=2)), β = $(round(beta_long; digits=2)))")

    # Display the plot
    display(plot!())

    # Print out the fit parameters
    println("Shorter Fit Parameters: τ = $(tau_short), β = $(beta_short)")
    println("Longer Fit Parameters: τ = $(tau_long), β = $(beta_long)")
end


autocorrelation_function_figures()

combined_L_beta_plot()

