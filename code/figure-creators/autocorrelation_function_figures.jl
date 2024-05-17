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

function autocorrelation_function_figures()
    
    simulation_name = "combined"

    ### --- READ IN DATA ---
    # Read this in for all the following filenames for different temperatures
    filenames = [
    # "L_11_T_0.6_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    # "L_11_T_0.8_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    "L_11_T_0.9_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
    "L_11_T_0.91_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
    # "L_11_T_0.92_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    # "L_11_T_0.93_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv", #
    "L_11_T_0.95_t_100000_1.0_configuration_autocorrelation_averages_by_time.csv",
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
        "L_11_T_0.9_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_1.0_t_130000_1.0_configuration_autocorrelation_averages_by_time.csv",
        "L_11_T_1.2_t_140000_1.0_configuration_autocorrelation_averages_by_time.csv"] 

    inset_temperatures = [0.9]





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
        plot!(graph, 0:length(autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)])-1, (1+c)*autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)] .- c, label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
        plot!(alternative_graph, 0:length(autocorrelation_functions_by_temperature[i][1:Int(altenative_lag_limit_cutoff)])-1, (1+c)*autocorrelation_functions_by_temperature[i][1:Int(altenative_lag_limit_cutoff)] .- c, label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
    end

    ## -- MINIMAL GRAPH ---
    for (i, temperature) in pairs(minimal_sample_temperatures)     
        plot!(minimal_graph, 0:length(minimal_autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)])-1, (1+c)*minimal_autocorrelation_functions_by_temperature[i][1:end-Int(lag_limit_relative_subtracted_cutoff)] .- c, label="T = "*string(temperature), color=temp_to_color[temperature], linewidth=2)
    end

    ### --- INSET TEMPERATURE GRAPH ON MINIMAL ---
    for inset_temperature in inset_temperatures
        temperature_index = findall(x -> x == inset_temperature, sample_temperatures)[1]

        plot!(minimal_graph, (1+c)*(autocorrelation_functions_by_temperature[temperature_index][1:end-Int(lag_limit_relative_subtracted_cutoff)] .- c),
        color=temp_to_color[inset_temperature], label="T=$(inset_temperature)", inset=bbox(0.25,0.44,0.28,0.28), subplot=2,
        xlabel=L"t", ylabel=L"\mathcal{C}(t)", yguidefontsize=10,xguidefontsize=10, yticks=[0.95,0.9,0.85], xticks=[0,6e4, 1.2e5], linewidth=2)
    end


    ### --- FITTING PARAMETER EXTRACTION ---
    tau_beta_by_temperature = []
    beta_by_temperature = []

    for (i, temperature) in pairs(sample_temperatures) 
        println("Temperature: ", temperature)


        # Fit the data to the autocorrelation function
        p0 = [10, 0.9]
        lb = [1e-3, 0.1] # example lower bounds
        ub = [1e11, 2.0]  # example upper bounds

        beta_relaxation_end = 1e4
        t_values = 1:beta_relaxation_end

        # Define the model function
        c = 0.201388888888888
        autocorrelation_model(t, p) = (1-c).* exp.(-(t./p[1]).^p[2]) .+ c

        if temperature >= 1.5

            # Fit the model to the data
            autocorrelation_fit = curve_fit(autocorrelation_model, t_values, autocorrelation_functions_by_temperature[i][1:Int(beta_relaxation_end)], p0, lower=lb, upper=ub)

            # Extract the parameters
            tau = autocorrelation_fit.param[1]
            beta = autocorrelation_fit.param[2]
            println("Tau = ", tau)
            println("Beta = ", beta)

            tau_beta_by_temperature = append!(tau_beta_by_temperature, tau)
            beta_by_temperature = append!(beta_by_temperature, beta)

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

        fit_temperatures = collect(LinRange(minimum(sample_temperatures), maximum(sample_temperatures), 100))

        tau_fits_graph = plot(title="", xlabel="Temperature, "*L"T", ylabel="Relaxation Time, "*L"\log_{10}(\tau)", legend=(0.85,0.3))
        scatter!(tau_fits_graph, sample_temperatures, log_tau_beta_by_temperature, color=alex_red, label="", markersize=5, ylims=extrema(log_tau_beta_by_temperature) .+ [-0.5, 0.5])

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

        println("Parabolic Parameters: ", parabolic_fit_parameters)
        println("Parabolic Residuals", parabolic_fit.resid)
        println("Parabolic Covariance Matrix", vcov(parabolic_fit))
        println("Parabolic Standard Errors", stderror(parabolic_fit))
        println("Parabolic Margin Error: ", margin_error(parabolic_fit))
        println("Confdience Intervals", confint(parabolic_fit))
        println("")
        println("")

        # Parabolic data over temperatures < 1.4
        alt_parabolic_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ T) .+ (p[3] ./ T.^2), sample_temperatures[1:findlast(x -> x < 1.4, sample_temperatures)], log_tau_beta_by_temperature[1:findlast(x -> x < 1.4, sample_temperatures)], [0.1, 0.1, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf])
        alt_parabolic_fit_parameters = alt_parabolic_fit.param
        alt_parabolic_fit_curve = alt_parabolic_fit_parameters[1] .+ (alt_parabolic_fit_parameters[2] ./ fit_temperatures) .+ (alt_parabolic_fit_parameters[3] ./ fit_temperatures.^2)
        lowest_confidence_interval_alt_parabolic_curve = confint(alt_parabolic_fit, level=conf_int)[1][1] .+ (confint(alt_parabolic_fit, level=conf_int)[2][1] ./ fit_temperatures) .+ (confint(alt_parabolic_fit, level=conf_int)[3][1] ./ fit_temperatures.^2)
        highest_confidence_interval_alt_parabolic_curve = confint(alt_parabolic_fit, level=conf_int)[1][2] .+ (confint(alt_parabolic_fit, level=conf_int)[2][2] ./ fit_temperatures) .+ (confint(alt_parabolic_fit, level=conf_int)[3][2] ./ fit_temperatures.^2)
        # plot!(tau_fits_graph, fit_temperatures, alt_parabolic_fit_curve, label="Parabolic Fit", color=alex_green, ribbon=(alt_parabolic_fit_curve .- lowest_confidence_interval_alt_parabolic_curve, highest_confidence_interval_alt_parabolic_curve .- alt_parabolic_fit_curve), fc=:green, fa=0.05)
        plot!(tau_fits_graph, fit_temperatures, alt_parabolic_fit_curve, label="Parabolic Fit", linestyle=:dash, color=alex_orange, linewidth=3)

        println("Alternative Parabolic Parameters: ", alt_parabolic_fit_parameters)
        println("Alternative Parabolic Residuals", alt_parabolic_fit.resid)
        println("Alternative Parabolic Covariance Matrix", vcov(alt_parabolic_fit))
        println("Alternative Parabolic Standard Errors", stderror(alt_parabolic_fit))
        println("Alternative Parabolic Margin Error: ", margin_error(alt_parabolic_fit))
        println("Alternative Parabolic Confidence Intervals", confint(alt_parabolic_fit))



        # VFT fit
        vft_fit = curve_fit((T, p) -> p[1] .+ (p[2] ./ (T .- p[3])), sample_temperatures, log_tau_beta_by_temperature, [1.0, 1.0, 0.1], lower=[-Inf, -Inf, -Inf], upper=[Inf, Inf, Inf], show_trace=true, maxIter=100)
        vft_fit_parameters = vft_fit.param
        vft_fit_curve = vft_fit_parameters[1] .+ (vft_fit_parameters[2] ./ (fit_temperatures .- vft_fit_parameters[3]))
        lowest_confidence_interval_vft_curve = confint(vft_fit, level=conf_int)[1][1] .+ (confint(vft_fit, level=conf_int)[2][1] ./ (fit_temperatures .- confint(vft_fit, level=conf_int)[3][1]))
        highest_confidence_interval_vft_curve = confint(vft_fit, level=conf_int)[1][2] .+ (confint(vft_fit, level=conf_int)[2][2] ./ (fit_temperatures .- confint(vft_fit, level=conf_int)[3][2]))

        # plot!(tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_orange, ribbon=(vft_fit_curve .- lowest_confidence_interval_vft_curve, highest_confidence_interval_vft_curve .- vft_fit_curve), fc=:blue, fa=0.05) # with confidence intervals
        plot!(tau_fits_graph, fit_temperatures, vft_fit_curve, label="VFT Fit", color=alex_green, linewidth=3)



        println("VFT Parameters: ", vft_fit_parameters)
        println("VFT Residuals", vft_fit.resid)
        println("VFT Covariance Matrix", vcov(vft_fit))
        println("VFT Standard Errors", stderror(vft_fit))
        println("VFT Margin Error: ", margin_error(vft_fit))
        println("Confdience Intervals", confint(vft_fit))
        println("")
        println("")


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

        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.png"))
        savefig(tau_fits_graph, joinpath("results/autocorrelation_anneal_results",simulation_name*"_relaxation_time_fits.svg"))
        display(tau_fits_graph)

    end
end
