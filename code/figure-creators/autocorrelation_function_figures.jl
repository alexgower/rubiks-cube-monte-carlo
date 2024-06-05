using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors, ColorSchemes
using LsqFit

include("../core/rubiks_cube.jl")


function quick_test()
    # Make some large tau stretched exponential curves with varying levels of noise and plot
    # taus = [1e7,3e7,1e8,3e8,1e9,3e9,1e10,3e10]
    taus = [1e7,1e8,1e9,1e10, 1e11]
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
    filename = "L_11_T_0.9_t_250000_trial_1_1.0_configuration_autocorrelation_averages_by_time.csv"

    data = readdlm(joinpath("results/autocorrelation_anneal_results/fast_results", filename), ',', skipstart=3)

    sample_temperature = data[1]
    samples_in_average = data[2]
    autocorrelation_functions_by_temperature = data[3:end]

    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])

    plot!(graph, 0:length(autocorrelation_functions_by_temperature)-1, autocorrelation_functions_by_temperature, label="T = "*string(sample_temperature), color=:black, linewidth=2)
    display(graph)

end

function fast_autocorrelation_function_figures()

    trials = 100
    autocorrelation_window_length=250000

    temperatures = [0.9,0.92,0.95]

    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for temperature in temperatures

        actual_number_of_trials = 0
        running_total_autocorrelation_function = zeros(autocorrelation_window_length+1-100)

        for trial in 1:trials
            println("Importing data for temperature: $(temperature), trial: $(trial)")

            filename = "results/autocorrelation_anneal_results/fast_results/"*"L_11_T_$(temperature)_t_$(autocorrelation_window_length)_trial_$(trial)_1.0_configuration_autocorrelation_averages_by_time.csv"
            try
                data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
            
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
    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=:topright)

    cutoff_end = Int(3e4)
    for temperature in temperatures
        autocorrelation_function = results_dictionary[temperature]
        plot!(graph, 0:length(autocorrelation_function)-1-cutoff_end, autocorrelation_function[1:end-cutoff_end], label="T = "*string(temperature), linewidth=2)
    end

    display(graph)
    savefig(graph, "results/autocorrelation_anneal_results/fast_results/autocorrelation_averages_by_time.png")
    savefig(graph, "results/autocorrelation_anneal_results/fast_results/autocorrelation_averages_by_time.svg")

    ## -- SAVE CONFIGURATION AVERAGES BY TIME --
    for temperature in temperatures
        filename = "L_11_T_$(temperature)_t_$(autocorrelation_window_length)_1.0_configuration_autocorrelation_averages_by_time.csv"

        touch(joinpath("results/autocorrelation_anneal_results/fast_results",filename))

        open(joinpath("results/autocorrelation_anneal_results/fast_results",filename), "w") do simulation_file
            write(simulation_file, "Simulation:L=11, P_s=1.0, T_1=1.0, T_0=1.0, N_T=100, autocorrelation_sample_size_per_temperature=100 ,autocorrelation_window_length=$autocorrelation_window_length \n")
            write(simulation_file, "Inherent Disorder Average \n")
            write(simulation_file, "Sample Temperature T, Samples in Average N, Configuration autocorrelation average by lag \n")
            write(simulation_file, "$(temperature), 100, $(join(results_dictionary[temperature], ", ")) \n")
        end
    end
    
end


function autocorrelation_function_figures()

    # TODO check if other autocorrelation stuff starts on time step 0
    
    simulation_name = "combined"

    ### --- READ IN DATA ---
    # Read this in for all the following filenames for different temperatures
    filenames = [
    # "L_11_T_0.6_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv", #
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
    "L_11_T_1.5_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_2.0_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", 
    "L_11_T_3.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv"] 

    all_temperatures = [0.6, 0.8, 0.9, 0.91, 0.92, 0.93, 0.95, 0.97, 1.0, 1.02, 1.05, 1.1, 1.15, 1.2, 1.25, 1.5, 2.0, 3.0]
    sample_temperatures_to_overlay_streched_exponential_fits_onto = all_temperatures



    minimal_filenames = [
        "L_11_T_0.6_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv", #
        "L_11_T_0.9_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_0.92_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_0.95_t_250000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_1.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_1.2_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv"] 

    minimal_graph_main_temperatures = [0.6, 0.9, 0.92, 1.0, 1.2]
    inset_temperatures = [0.9,0.92,0.95]





    sample_temperatures = []
    samples_in_average = []
    autocorrelation_functions_by_temperature = []

    # For each file in filenames, read in the data as a matrix
    for (i, filename) in pairs(filenames)
        data = readdlm(joinpath("results/autocorrelation_anneal_results", filename), ',', skipstart=3)

        sample_temperatures = append!(sample_temperatures, data[1])
        samples_in_average = append!(samples_in_average, data[2])
        autocorrelation_functions_by_temperature = append!(autocorrelation_functions_by_temperature, [data[3:end]])
    end

    minimal_sample_temperatures = []
    minimal_samples_in_average = []
    minimal_autocorrelation_functions_by_temperature = []

    # For each file in filenames, read in the data as a matrix
    for (i, filename) in pairs(minimal_filenames)
        data = readdlm(joinpath("results/autocorrelation_anneal_results", filename), ',', skipstart=3)

        minimal_sample_temperatures = append!(minimal_sample_temperatures, data[1])
        minimal_samples_in_average = append!(minimal_samples_in_average, data[2])
        minimal_autocorrelation_functions_by_temperature = append!(minimal_autocorrelation_functions_by_temperature, [data[3:end]])
    end

    ### --- COLOURS ---
    Plots.default(dpi = 300)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    # Define a fixed table of temperatures and their corresponding colors
    temp_to_color = Dict()
    alex_colors = [alex_red, alex_orange, alex_pink, alex_green, alex_blue]

    for temperature in all_temperatures
            temp_to_color[temperature] = alex_colors[mod1(findall(x -> x == temperature, all_temperatures)[1], length(alex_colors))]
    end
    temp_to_color[0.6] = alex_red
    temp_to_color[0.9] = alex_pink
    temp_to_color[1.0] = alex_orange
    temp_to_color[1.2] = alex_green


    ### --- PLOTTING ---
    ## -- COMBINED GRAPHS ---
    graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    alternative_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    minimal_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=(0.9,0.7), ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])

    c = 0.201388888888888
    lag_limit_relative_subtracted_cutoff = 1e4
    altenative_lag_limit_cutoff = 5e4

    for (i, temperature) in pairs(sample_temperatures)     
        plot!(graph, 0:length(autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)])-1, (1+c)*autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)] .- c, label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
        plot!(alternative_graph, 0:length(autocorrelation_functions_by_temperature[i][1:Int(altenative_lag_limit_cutoff)])-1, (1+c)*autocorrelation_functions_by_temperature[i][1:Int(altenative_lag_limit_cutoff)] .- c, label="T = $(temperature)", color=temp_to_color[temperature], linewidth=2)
    end

    ## -- MINIMAL GRAPH ---
    main_minimal_graph_window = Int(1.2e5)
    for (i, temperature) in pairs(minimal_graph_main_temperatures)     
        plot!(minimal_graph, 0:Int(main_minimal_graph_window-1), (1+c)*minimal_autocorrelation_functions_by_temperature[i][1:minimal_autocorrelation_functions_by_temperature] .- c, label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
    end

    ### --- INSET TEMPERATURE GRAPH ON MINIMAL ---
    for inset_temperature in inset_temperatures
        temperature_index = findall(x -> x == inset_temperature, sample_temperatures)[1]

        # Plot same data on an inset graph within the minimal graph
        inset_graph_window_size = Int(2.2e5)
        # inset_graph = plot(title="", xlabel="Time, $(L"t") [MC Steps]", ylabel="Autocorrelation Function, "*L"\mathcal{C}(t)", legend=:topright, ylims = (0.0,1.05), yticks=[0.0,0.2,0.4,0.6,0.8,1.0])
        # plot!(inset_graph, 0:inset_graph_window_size-1, (1+c)*autocorrelation_functions_by_temperature[temperature_index][1:inset_graph_window_size] .- c, color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", linewidth=2)
        plot!(minimal_graph, 0:main_minimal_graph_window-1, (1+c)*minimal_autocorrelation_functions_by_temperature[temperature_index][1:main_minimal_graph_window] .- c, color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
        xlabel=L"t", ylabel=L"\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10,  xticks=[0,6e4,1.2e5], linewidth=2)        
        
        # plot!(minimal_graph, (1+c)*(minimal_autocorrelation_functions_by_temperature[temperature_index][1:end-Int(lag_limit_relative_subtracted_cutoff)] .- c),
        # color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
        # xlabel=L"t", ylabel=L"\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10,  xticks=[0,6e4,1.2e5], linewidth=2)
    end


    ### --- FITTING PARAMETER EXTRACTION ---
    tau_beta_by_temperature = []
    beta_by_temperature = []

    longer_tau_beta_by_temperature = []
    longer_beta_by_temperature = []

    for (i, temperature) in pairs(sample_temperatures) 
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
        autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

        if temperature >= 1.2

            # Fit the model to the data
            autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[i][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

            # Extract the parameters
            tau = autocorrelation_fit.param[1]
            beta = autocorrelation_fit.param[2]
            println("Tau = ", tau)
            println("Beta = ", beta)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
            beta_by_temperature = append!(beta_by_temperature, beta)


            longer_autocorrelation_fit = curve_fit(autocorrelation_model, longer_t_values, autocorrelation_functions_by_temperature[i][1:Int(longer_beta_relaxation_end)], p0, lower=lb, upper=ub)
            longer_tau = longer_autocorrelation_fit.param[1]
            longer_beta = longer_autocorrelation_fit.param[2]
            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_beta)


        else
            # Now try fitting on log of data with offset already removed
            reduced_data = log.((1-c)^(-1) * (autocorrelation_functions_by_temperature[i][1:Int(beta_relaxation_end)] .- c))

            reduced_model(t,p) = -(t./p[1]).^p[2]

            autocorrelation_fit = curve_fit(reduced_model, t_values, reduced_data, p0, lower=lb, upper=ub)

            reduced_tau = autocorrelation_fit.param[1]
            reduced_beta = autocorrelation_fit.param[2]
            println("Tau From Logged Data = ", reduced_tau)
            println("Beta from Logged Data = ", reduced_beta)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, reduced_tau)
            beta_by_temperature = append!(beta_by_temperature, reduced_beta)


            longer_reduced_data = log.((1-c)^(-1) * (autocorrelation_functions_by_temperature[i][1:Int(longer_beta_relaxation_end)] .- c))

            longer_autocorrelation_fit = curve_fit(reduced_model, longer_t_values, longer_reduced_data, p0, lower=lb, upper=ub)
            longer_reduced_tau = longer_autocorrelation_fit.param[1]
            longer_reduced_beta = longer_autocorrelation_fit.param[2]
            longer_tau_beta_by_temperature = append!(longer_tau_beta_by_temperature, longer_reduced_tau)
            longer_beta_by_temperature = append!(longer_beta_by_temperature, longer_reduced_beta)


        end

        # Plot the fit
        dashed_cutoff = beta_relaxation_end
        alternative_dashed_cutoff = 5e4

        if temperature in sample_temperatures_to_overlay_streched_exponential_fits_onto
            # label="T = "*string(temperature)*", τ = "*string(tau_beta)*", β = "*string(beta)
            plot!(graph, (1+c)*autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c, color=:black, linestyle=:dash, label="")
            plot!(alternative_graph, (1+c)*autocorrelation_model(1:Int(alternative_dashed_cutoff), autocorrelation_fit.param) .- c, color=:black, linestyle=:dash, label="")

            if temperature in minimal_sample_temperatures
                plot!(minimal_graph, (1+c)*autocorrelation_model(1:Int(dashed_cutoff), autocorrelation_fit.param) .- c, color=:black, linestyle=:dash, label="")
            end
        end

        # Plot the longer fit
        longer_dashed_cutoff = longer_beta_relaxation_end
        if temperature in sample_temperatures_to_overlay_streched_exponential_fits_onto
            plot!(graph, (1+c)*autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c, color=:green, linestyle=:dash, label="")
            plot!(alternative_graph, (1+c)*autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c, color=:green, linestyle=:dash, label="")
        
            if temperature in minimal_sample_temperatures
                plot!(minimal_graph, (1+c)*autocorrelation_model(1:Int(longer_dashed_cutoff), longer_autocorrelation_fit.param) .- c, color=:green, linestyle=:dash, label="")
            end
        end

    end

    # --- SAVE GRAPHS ---
    savefig(graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.png"))
    savefig(graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time.svg"))
    display(graph)

    savefig(alternative_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.png"))
    savefig(alternative_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_alternative.svg"))
    display(alternative_graph)

    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.png"))
    savefig(minimal_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_autocorrelation_averages_by_time_minimal.svg"))
    display(minimal_graph)






    if simulation_name == "combined"
        # -- RELAXATION TIME FITS ---
        # Make plot with tau_beta_by_temperature against temperature 
        # and then Arrhenius/VFT/parabolic fit curves to data too and output measure
        # Do fits on log of data to avoid numerical errors
        log_tau_beta_by_temperature = log10.(tau_beta_by_temperature) # Note using log10 here
        longer_log_tau_beta_by_temperature = log10.(longer_tau_beta_by_temperature)

        fit_temperatures = collect(LinRange(minimum(sample_temperatures)-0.03, maximum(sample_temperatures), 100))

        tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))
        scatter!(tau_fits_graph, sample_temperatures, log_tau_beta_by_temperature, color=alex_red, label="", markersize=5, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5], alpha=0.8)
        scatter!(tau_fits_graph, sample_temperatures, longer_log_tau_beta_by_temperature, color=alex_red, label="", markersize=4, marker=:square, alpha=0.8)

        conf_int = 0.05


        # # Arrhenius fit
        # arrhenius_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T), sample_temperatures, log_tau_beta_by_temperature, [0.1, 0.1], lower=[-Inf,-Inf], upper=[Inf, Inf])
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
        parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), sample_temperatures, log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        parabolic_fit_parameters = parabolic_fit.param
        parabolic_fit_curve = parabolic_fit_parameters[1] .+ (parabolic_fit_parameters[2] ./ fit_temperatures) .+ (parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        lowest_confidence_interval_parabolic_curve = confint(parabolic_fit, level=conf_int)[1][1] .+ (confint(parabolic_fit, level=conf_int)[2][1] ./ fit_temperatures) .+ (confint(parabolic_fit, level=conf_int)[3][1] ./ fit_temperatures.^2)
        highest_confidence_interval_parabolic_curve = confint(parabolic_fit, level=conf_int)[1][2] .+ (confint(parabolic_fit, level=conf_int)[2][2] ./ fit_temperatures) .+ (confint(parabolic_fit, level=conf_int)[3][2] ./ fit_temperatures.^2)
        # plot!(tau_fits_graph, fit_temperatures, parabolic_fit_curve, label="Parabolic Fit", color=alex_blue, ribbon=(parabolic_fit_curve .- lowest_confidence_interval_parabolic_curve, highest_confidence_interval_parabolic_curve .- parabolic_fit_curve), fc=:green, fa=0.05) # with confidence interval
        plot!(tau_fits_graph, fit_temperatures, parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        println("")
        println("")
        println("Parabolic Parameters: ", parabolic_fit_parameters)
        println("Parabolic Residuals", parabolic_fit.resid)
        println("Parabolic Covariance Matrix", vcov(parabolic_fit))
        println("Parabolic Standard Errors", stderror(parabolic_fit))
        println("Parabolic Margin Error: ", margin_error(parabolic_fit))
        println("Confdience Intervals", confint(parabolic_fit))

        longer_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), sample_temperatures, longer_log_tau_beta_by_temperature, [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        longer_parabolic_fit_parameters = longer_parabolic_fit.param
        longer_parabolic_fit_curve = longer_parabolic_fit_parameters[1] .+ (longer_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        # plot!(tau_fits_graph, fit_temperatures, longer_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_blue, linewidth=3)

        print("Longer Parabolic Parameters: ", longer_parabolic_fit_parameters)



        # Parabolic data over temperatures < 1.4
        alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), sample_temperatures[1:findlast(x -> x < 1.4, sample_temperatures)], log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, sample_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        alt_parabolic_fit_parameters = alt_parabolic_fit.param
        alt_parabolic_fit_curve = alt_parabolic_fit_parameters[1] .+ (alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        lowest_confidence_interval_alt_parabolic_curve = confint(alt_parabolic_fit, level=conf_int)[1][1] .+ (confint(alt_parabolic_fit, level=conf_int)[2][1] ./ fit_temperatures) .+ (confint(alt_parabolic_fit, level=conf_int)[3][1] ./ fit_temperatures.^2)
        highest_confidence_interval_alt_parabolic_curve = confint(alt_parabolic_fit, level=conf_int)[1][2] .+ (confint(alt_parabolic_fit, level=conf_int)[2][2] ./ fit_temperatures) .+ (confint(alt_parabolic_fit, level=conf_int)[3][2] ./ fit_temperatures.^2)
        # plot!(tau_fits_graph, fit_temperatures, alt_parabolic_fit_curve, label="Parabolic Fit", color=alex_green, ribbon=(alt_parabolic_fit_curve .- lowest_confidence_interval_alt_parabolic_curve, highest_confidence_interval_alt_parabolic_curve .- alt_parabolic_fit_curve), fc=:green, fa=0.05)
        plot!(tau_fits_graph, fit_temperatures, alt_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_orange, linewidth=3)

        println("")
        println("")
        println("Alternative Parabolic Parameters: ", alt_parabolic_fit_parameters)
        println("Alternative Parabolic Residuals", alt_parabolic_fit.resid)
        println("Alternative Parabolic Covariance Matrix", vcov(alt_parabolic_fit))
        println("Alternative Parabolic Standard Errors", stderror(alt_parabolic_fit))
        println("Alternative Parabolic Margin Error: ", margin_error(alt_parabolic_fit))
        println("Alternative Parabolic Confidence Intervals", confint(alt_parabolic_fit))

        longer_alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), sample_temperatures[1:findlast(x -> x < 1.4, sample_temperatures)], longer_log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, sample_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        longer_alt_parabolic_fit_parameters = longer_alt_parabolic_fit.param
        longer_alt_parabolic_fit_curve = longer_alt_parabolic_fit_parameters[1] .+ (longer_alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (longer_alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        # plot!(tau_fits_graph, fit_temperatures, longer_alt_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_orange, linewidth=3)

        print("Longer Alternative Parabolic Parameters: ", longer_alt_parabolic_fit_parameters)


        # VFT fit
        vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), sample_temperatures, log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        vft_fit_parameters = vft_fit.param
        vft_fit_curve = vft_fit_parameters[1] .+ (vft_fit_parameters[2] ./ (fit_temperatures .- vft_fit_parameters[3]))
        lowest_confidence_interval_vft_curve = confint(vft_fit, level=conf_int)[1][1] .+ (confint(vft_fit, level=conf_int)[2][1] ./ (fit_temperatures .- confint(vft_fit, level=conf_int)[3][1]))
        highest_confidence_interval_vft_curve = confint(vft_fit, level=conf_int)[1][2] .+ (confint(vft_fit, level=conf_int)[2][2] ./ (fit_temperatures .- confint(vft_fit, level=conf_int)[3][2]))
        # plot!(tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_orange, ribbon=(vft_fit_curve .- lowest_confidence_interval_vft_curve, highest_confidence_interval_vft_curve .- vft_fit_curve), fc=:blue, fa=0.05) # with confidence intervals
        plot!(tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        println("")
        println("")
        println("VFT Parameters: ", vft_fit_parameters)
        println("VFT Residuals", vft_fit.resid)
        println("VFT Covariance Matrix", vcov(vft_fit))
        println("VFT Standard Errors", stderror(vft_fit))
        println("VFT Margin Error: ", margin_error(vft_fit))
        println("Confdience Intervals", confint(vft_fit))


        longer_vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), sample_temperatures, longer_log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        longer_vft_fit_parameters = longer_vft_fit.param
        longer_vft_fit_curve = longer_vft_fit_parameters[1] .+ (longer_vft_fit_parameters[2] ./ (fit_temperatures .- longer_vft_fit_parameters[3]))
        # plot!(tau_fits_graph, fit_temperatures, longer_vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)

        print("Longer VFT Parameters: ", longer_vft_fit_parameters)


        #--- FRAGILITY INSET DATA ---
        scatter!(tau_fits_graph, [1/T for T in sample_temperatures], log_tau_beta_by_temperature,
        color=alex_red, legend=false, inset=bbox(0.65,0.05,0.3,0.35), subplot=2,
        xlabel=L"1/T", ylabel=L"\log_{10}(\tau)", yguidefontsize=12,xguidefontsize=12, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5])

        # Plot fits in inset too
        # plot!(tau_fits_graph, [1/T for T in fit_temperatures], arrhenius_fit_curve,color=alex_orange, label="Arrhenius Fit", linestyle=:dash, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], parabolic_fit_curve,color=alex_blue, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], alt_parabolic_fit_curve,color=alex_orange, label="Parabolic Fit", linestyle=:dash, linewidth=3, subplot=2)
        plot!(tau_fits_graph, [1/T for T in fit_temperatures], vft_fit_curve,color=alex_green, label="VFT Fit", linewidth=3, subplot=2)


        # --- BETA INSET DATA ---
        plot!(tau_fits_graph, sample_temperatures, beta_by_temperature,
        color=alex_pink, marker=:circle, legend=false, inset=bbox(0.23,0.05,0.3,0.35), subplot=3,
        xlabel=L"T", ylabel=L"\beta", yguidefontsize=12,xguidefontsize=12, linewidth=3)

        # Plot longer beta data
        # plot!(tau_fits_graph, sample_temperatures, longer_beta_by_temperature,
        # color=alex_pink, marker=:square, legend=false, subplot=3)

        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.png"))
        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.svg"))
        display(tau_fits_graph)

    end
end
