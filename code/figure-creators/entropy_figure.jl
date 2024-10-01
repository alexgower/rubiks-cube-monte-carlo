using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

using Images

using CubicSplines
using QuadGK
using CurveFit
using Statistics


include("../core/rubiks_cube.jl")

function entropy_figure()

    # --- All L Figure ---
    # models = ["clean", "inherent_disorder"]
    models = ["clean", "inherent_disorder"]
    Ls = [11]
    swap_move_probabilities = [1.0]
    trials = 50
    N_T = 100


    ### --- COLOURS ---
    Plots.default(dpi = 300)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)

    alex_alt_blue = RGB(4/255, 57/255, 94/255)


    ### --- READ IN DATA ---
    filenames_that_do_not_exist=[]

    results_dictionary = Dict()
    for model in models
        for L in Ls
            for swap_move_probability in swap_move_probabilities
                N_T = (model=="clean" && L==11) ? 200 : 100

                temperatures = zeros(N_T)
                running_total_specific_heat_capacities = zeros(N_T)
                running_total_energies = zeros(N_T)

                actual_number_of_trials=0
                for trial in 1:trials
                    filename = "results/relaxed_anneal_results/data/" * model * "_L_" * string(L) * "_trial_" * string(trial) * "_$(swap_move_probability)"

                    # try
                        data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=3)
                        
                        temperatures .= data_matrix[:,1]
                        running_total_energies .+= data_matrix[:,2]
                        running_total_specific_heat_capacities .+= data_matrix[:,5] 

                        actual_number_of_trials += 1
                    # catch e
                    #     push!(filenames_that_do_not_exist, filename)
                    # end
                
                end

                results_dictionary[(model, L, swap_move_probability, "temperatures")] = temperatures
                results_dictionary[(model, L, swap_move_probability, "average_energies")] = running_total_energies/actual_number_of_trials
                results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")] = running_total_specific_heat_capacities/actual_number_of_trials
            end
        end
    end











    color_index = 1
    for model in models
        if model=="clean"
            colors =  [alex_alt_blue, alex_green, alex_blue, alex_grey]
        elseif model=="inherent_disorder" 
            colors = [alex_pink, alex_red, alex_orange]
        elseif model=="custom"
            colors = [alex_blue]
        end

        for L in Ls
            for swap_move_probability in swap_move_probabilities

                # Gather read in data
                temperatures = results_dictionary[(model, L, swap_move_probability, "temperatures")]
                average_specific_heat_capacities = results_dictionary[(model, L, swap_move_probability, "average_specific_heat_capacities")]
                average_energies = results_dictionary[(model, L, swap_move_probability, "average_energies")]

                # SPLINES
                temperature_sorted_indices = sortperm(temperatures)
                sorted_temperatures = temperatures[temperature_sorted_indices]
                temperature_sorted_specific_heat_capacities = average_specific_heat_capacities[temperature_sorted_indices]
                temperature_sorted_energies = average_energies[temperature_sorted_indices]

                spline_specific_heat_capacity_from_temperature = CubicSpline(sorted_temperatures, temperature_sorted_specific_heat_capacities, extrapl=[1,], extrapr=[1,])
                println("Spline specific heat capacity from temperature at 0.1: $(spline_specific_heat_capacity_from_temperature(0.1))")
                spline_energy_from_temperature = CubicSpline(sorted_temperatures, temperature_sorted_energies, extrapl=[1,], extrapr=[1,])
                println("Spline energy from temperature at 0.1: $(spline_energy_from_temperature(0.1))")
                spline_absolute_energy_from_temperature = CubicSpline(sorted_temperatures, temperature_sorted_energies .+ abs(minimum(temperature_sorted_energies)), extrapl=[1,], extrapr=[1,])
                println("Spline absolute energy from temperature at 0.1: $(spline_absolute_energy_from_temperature(0.1))")

                energy_sorted_indices = sortperm(average_energies)
                sorted_absolute_energies = average_energies[energy_sorted_indices] .+ abs(minimum(average_energies))
                energy_sorted_temperatures = temperatures[energy_sorted_indices]

                spline_temperature_from_absolute_energy = CubicSpline(sorted_absolute_energies, energy_sorted_temperatures, extrapl=[1,], extrapr=[1,])





                # SPLINE GRAPHS
                energy_graph = plot(temperatures, average_energies, label="Energy", color=:red, linestyle=:solid)
                heat_capacity_graph = plot(temperatures, average_specific_heat_capacities, label="Specific Heat Capacity", color=:red, linestyle=:solid)
                
                test_temperatures_start_point = model == "clean" ? 0.2 : 0.01
                test_temperatures = range(sorted_temperatures[1]+test_temperatures_start_point, sorted_temperatures[end], length=1000)
                
                plot!(energy_graph, test_temperatures, spline_energy_from_temperature(test_temperatures), label="Energy Spline", color=:blue, linestyle=:dash, xlabel="Temperature", ylabel="Energy")
                plot!(heat_capacity_graph, test_temperatures, spline_specific_heat_capacity_from_temperature(test_temperatures), label="Specific Heat Capacity Spline", color=:blue, linestyle=:dash, xlabel="Temperature", ylabel="Specific Heat Capacity")


                display(energy_graph)
                display(heat_capacity_graph)





                # ENTROPY CALCULATIONS
                # First using specific heat capacity spline over 1000 temperatures \int_0^T \frac{C_V(T')}{T'} dT'
                entropy_ala_heat_capacity_by_test_temperature = zeros(1000)
                for i in 1:1000
                    entropy_ala_heat_capacity_by_test_temperature[i] = quadgk(x->spline_specific_heat_capacity_from_temperature(x)/x, 0.001, test_temperatures[i])[1]
                end
                spline_entropy_ala_heat_capacity_from_temperature = CubicSpline(test_temperatures, entropy_ala_heat_capacity_by_test_temperature, extrapl=[1,], extrapr=[1,])


                # Second using energy spline over 1000 temperatures \int_0^E \frac{1}{T(E')} dE'
                entropy_ala_energy_by_energy = zeros(1000)
                for i in 1:1000
                    entropy_ala_energy_by_energy[i] = quadgk(x->1/spline_temperature_from_absolute_energy(x), 0.1, spline_absolute_energy_from_temperature(test_temperatures[i]))[1]
                end
                sorted_absolute_energies_by_test_temperature_indices = sortperm(spline_absolute_energy_from_temperature(test_temperatures))
                sorted_absolute_energies_by_test_temperature = spline_absolute_energy_from_temperature(test_temperatures)[sorted_absolute_energies_by_test_temperature_indices]
                sorted_entropy_ala_energy_by_energy = entropy_ala_energy_by_energy[sorted_absolute_energies_by_test_temperature_indices]
                spline_entropy_ala_energy_from_absolute_energy = CubicSpline(sorted_absolute_energies_by_test_temperature, sorted_entropy_ala_energy_by_energy, extrapl=[1,], extrapr=[1,])


                # PLOT S(T) Combined
                # entropy_temperature_graph = plot(test_temperatures, entropy_ala_heat_capacity_by_test_temperature, label="Specific Heat Capacity Entropy", color=:red, linestyle=:solid, xlabel="Temperature", ylabel="Entropy")
                # plot!(entropy_temperature_graph, test_temperatures, spline_entropy_ala_energy_from_absolute_energy(spline_absolute_energy_from_temperature(test_temperatures)), label="Energy Entropy", color=:blue, linestyle=:solid, xlabel="Temperature", ylabel="Entropy")
                # display(entropy_temperature_graph)

                # PLOT S(E) Combined
                entropy_energy_graph = plot(spline_absolute_energy_from_temperature(test_temperatures), spline_entropy_ala_heat_capacity_from_temperature(test_temperatures), label="Specific Heat Capacity Entropy", color=:red, linestyle=:solid, xlabel="Energy", ylabel="Entropy")
                plot!(entropy_energy_graph, spline_absolute_energy_from_temperature(test_temperatures), entropy_ala_energy_by_energy, label="Entropy", color=:blue, linestyle=:solid)

                display(entropy_energy_graph)


                # S(E) from energy is the better one so use that
                legend_label = model == "clean" ? "L = 11, Original RC" : "L = 11, Randomised RC"

                final_graph = plot(spline_absolute_energy_from_temperature(test_temperatures), entropy_ala_energy_by_energy, label=legend_label, color=:blue, linestyle=:solid, xlabel="Absolute Energy, "*L"E - E_0", ylabel="Entropy, "*L"S \sim \ln(\mathcal{N}(E))")


                if model == "clean"
                    maximum_energy_to_fit_to = 600
                else
                    maximum_energy_to_fit_to = 600
                end

                indices_to_fit_to = findall(x -> x < maximum_energy_to_fit_to, spline_absolute_energy_from_temperature(test_temperatures))
                x = spline_absolute_energy_from_temperature(test_temperatures)[indices_to_fit_to]
                y = entropy_ala_energy_by_energy[indices_to_fit_to]

                # Perform linear fit
                intercept, gradient = linear_fit(x, y)

                # Calculate residuals and standard errors
                y_fit = gradient .* x .+ intercept
                residuals = y .- y_fit
                n = length(x)
                mse = sum(residuals.^2) / (n - 2)
                x_mean = mean(x)
                gradient_se = sqrt(mse / sum((x .- x_mean).^2))
                intercept_se = sqrt(mse * (1/n + x_mean^2 / sum((x .- x_mean).^2)))

                println("Gradient: $(round(gradient, digits=6)) ± $(round(gradient_se, digits=6))")
                println("Intercept: $(round(intercept, digits=6)) ± $(round(intercept_se, digits=6))")

                # Add linear fit to graph
                plot!(final_graph, x, y_fit, label="Linear Fit, S = $(round(gradient, digits=2))E + $(round(intercept, digits=2))", color=:red, linestyle=:dash)

                display(final_graph)

                savefig(final_graph, "results/relaxed_anneal_results/entropy_figure_L_$(L)_$(model).png")
                savefig(final_graph, "results/relaxed_anneal_results/entropy_figure_L_$(L)_$(model).pdf")


            end
            color_index += 1
        end

    end


    ### --- PRINT ERRORS ---
    println("The following files do not exist:")
    for filename in filenames_that_do_not_exist
        println(filename)
    end

end
