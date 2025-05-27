using Plots
using LsqFit
using Printf
using LaTeXStrings
using Plots.PlotMeasures

function swap_moves_growth_figure()
    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)

    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Data
    L = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    N_s = [168, 2588, 2149196, 281635556, 44845673155896, 5940848633430244,
           260219508822162622045800, 34350321361181406176576276,
           416759185539018736945142423356183016, 55012222938456767359870599949634752100]

    odd_L = [ 3, 5, 7, 9, 11]
    N_s_odd = [2588, 281635556, 5940848633430244, 34350321361181406176576276, 55012222938456767359870599949634752100]

    even_L = [2, 4, 6, 8, 10]
    N_s_even = [168, 2149196, 44845673155896, 260219508822162622045800, 416759185539018736945142423356183016]
    
    x = L.^2
    y = log.(N_s)
    
    # Combined Linear fit
    model(x, p) = p[1] .+ p[2]*x
    p0 = [0.0, 0.0]  # Initial guess
    fit = curve_fit(model, x, y, p0)
    
    # Generate points for the fitted curve
    x_fit = range(minimum(x), maximum(x), length=100)
    y_fit = model(x_fit, fit.param)



    # Odd Linear fit
    x_odd = odd_L.^2
    y_odd = log.(N_s_odd)
    fit_odd = curve_fit(model, x_odd, y_odd, p0)
    x_fit_odd = range(minimum(x_odd), maximum(x_odd), length=100)
    y_fit_odd = model(x_fit_odd, fit_odd.param)

    # Even Linear fit
    x_even = even_L.^2
    y_even = log.(N_s_even)
    fit_even = curve_fit(model, x_even, y_even, p0)
    x_fit_even = range(minimum(x_even), maximum(x_even), length=100)
    y_fit_even = model(x_fit_even, fit_even.param)


    
    # Plot
    graph = plot(x_even, y_even, seriestype=:scatter, label=L"N_s, \rm{even \ L}", markersize=5, color=alex_red)
    plot!(graph, x_odd, y_odd, seriestype=:scatter, label=L"N_s, \rm{odd \ L}", markersize=5, color=alex_orange)
    
    plot!(graph,x_fit, y_fit, label=L"\ln(N_s) \approx 2.53 + 0.73 L^2", lw=2, color=alex_green, xlabel=L"L^2", ylabel=L"\ln(N_s)", legend=:bottomright, margin=5mm)
    plot!(graph, x_fit_odd, y_fit_odd, label=L"\ln(N_s) \approx 1.70 + 0.70 L^2, \rm{odd \ L}", lw=2, color=alex_orange)
    plot!(graph, x_fit_even, y_fit_even, label=L"\ln(N_s) \approx 2.07 + 0.80 L^2, \rm{even \ L}", lw=2, color=alex_red)


    # Plot odd L graph only 
    graph_odd = plot(x_odd, y_odd, seriestype=:scatter, label=L"N_s", markersize=5, color=alex_orange)
    plot!(graph_odd, x_fit_odd, y_fit_odd, label=L"\ln(N_s) \approx 1.70 + 0.70 L^2", lw=2, color=alex_red, xlabel=L"L^2", ylabel=L"\ln(N_s)", legend=:bottomright, margin=5mm)
    
    display(graph)
    savefig(graph, "results/miscellaneous/swap_moves_growth_figure.png")
    savefig(graph, "results/miscellaneous/swap_moves_growth_figure.pdf")

    display(graph_odd)
    savefig(graph_odd, "results/miscellaneous/swap_moves_growth_figure_odd.png")
    savefig(graph_odd, "results/miscellaneous/swap_moves_growth_figure_odd.pdf")
    
    # Print the fitted parameters
    println("Fitted parameters:")
    println("a = ", fit.param[1])
    println("b = ", fit.param[2])

    println("Fitted parameters (odd):")
    println("a = ", fit_odd.param[1])
    println("b = ", fit_odd.param[2])

    println("Fitted parameters (even):")
    println("a = ", fit_even.param[1])
    println("b = ", fit_even.param[2])
end

function configurations_growth_figure()
    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)

    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # Data
    L = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    N_tot = [
        BigInt(1),
        BigInt(3674160),
        BigInt(43252003274489856000),
        BigInt(7.4e45),
        BigInt(2.8e74),
        BigInt(1.6e116),
        BigInt(2.0e160),
        BigInt(3.5e217),
        BigInt(1.4e277),
    ]
    
    # Remove L=1 as log(0) is undefined
    L = L[2:end]
    N_tot = N_tot[2:end]
    
    x = L.^2
    y = log.(N_tot)
    
    # Linear fit
    model(x, p) = p[1] .+ p[2]*x
    p0 = [0.0, 0.0]  # Initial guess
    fit = curve_fit(model, x, y, p0)
    x_fit = range(minimum(x), maximum(x), length=100)
    y_fit = model(x_fit, fit.param)

    # Odd L Linear Fit
    odd_L = [3, 5, 7, 9]
    N_tot_odd = [BigInt(43252003274489856000), BigInt(2.8e74), BigInt(2.0e160), BigInt(1.4e277)]
    x_odd = odd_L.^2
    y_odd = log.(N_tot_odd)
    fit_odd = curve_fit(model, x_odd, y_odd, p0)
    x_fit_odd = range(minimum(x_odd), maximum(x_odd), length=100)
    y_fit_odd = model(x_fit_odd, fit_odd.param)
    
    # Plot
    graph = plot(x, y, seriestype=:scatter, label=L"N_{tot}", xlabel=L"L^2", ylabel=L"\ln(N_{tot})", legend=:bottomright, margin=6mm, markersize=5, color=alex_blue)
    plot!(graph, x_fit, y_fit, label=L"\ln(N_{tot}) \approx -25.9 + 8.17 L^2", lw=2, color=alex_grey)

    # Plot odd L graph only
    graph_odd = plot(x_odd, y_odd, seriestype=:scatter, label=L"N_{tot}", xlabel=L"L^2", ylabel=L"\ln(N_{tot})", legend=:bottomright, margin=6mm, markersize=5, color=alex_blue)
    plot!(graph_odd, x_fit_odd, y_fit_odd, label=L"\ln(N_{tot}) \approx -32.3 + 8.25 L^2", lw=2, color=alex_grey)


    display(graph)
    savefig(graph, "results/miscellaneous/configurations_growth_figure.png")
    savefig(graph, "results/miscellaneous/configurations_growth_figure.pdf")

    display(graph_odd)
    savefig(graph_odd, "results/miscellaneous/configurations_growth_figure_odd.png")
    savefig(graph_odd, "results/miscellaneous/configurations_growth_figure_odd.pdf")

    # Print the fitted parameters
    println("Fitted parameters for N_tot:")
    println("a = ", fit.param[1])
    println("b = ", fit.param[2])

    println("Fitted parameters for N_tot (odd):")
    println("a = ", fit_odd.param[1])
    println("b = ", fit_odd.param[2])


end

swap_moves_growth_figure()
configurations_growth_figure()