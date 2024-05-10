# Illustrating the Newton-Raphson iteration scheme to find the minimum of a given f(x).
using Pkg
Pkg.activate(".")
using CairoMakie, MathTeXEngine, Zygote

function newtonRaphson()
    
    # Function definitions
    f(x) =  x.^2 .- 2.0                                         # Function to minimize
    df(x) = 2.0 * x                                             # Its analytical derivative
    tangent_line(f, x, x0) = f(x0) .+ f'(x0) .* (x .- x0)       # Function returning points that ly on the tangent line at any point in f(x)
    x0   = 10.0                                                 # Initial guess
    
    # Initialisation of variables
    xn_z     = copy(x0)                                         # Zygote
    xn_n     = copy(x0)                                         # Numerical
    xn_a     = copy(x0)                                         # Analytical copy of initial guess
    xn_z_all = [x0]                                             # History arrays for the above variables
    xn_a_all = [x0]
    xn_n_all = [x0]
    xplot    = [0.0 + i*(x0 / 100) for i in 1:100]              # Array for plotting
    α        = [0.1 + i*(1.0 / 8.0) for i in 1:8]               # Transparency of markers and lines

    # Solver settings
    iter_max = 1000                                             # Maximum no. iterations
    ϵ_tol    = 1e-11                                            # Solver tolerance

    # Visualize
    fg1 = Figure(fontsize = 18)
    ax1 = Axis(fg1[1,1], xlabel = L"x_n", ylabel = L"f(x)", limits = (0.0, x0, -10.0, f(x0) + 2.0), aspect = 1)
    ax2 = Axis(fg1[1,2], xlabel = L"n", ylabel = L"| x_{n+1} - x_n |", yscale = log10, limits = (0, 10, ϵ_tol / 10.0, 100.0), aspect = 1)
    ln3 = lines!(ax1, xplot, (f.(xplot)), label = "f(x)", color = :steelblue, linewidth = 3)

    # Iteration loop
    err_a = 2.0 * ϵ_tol; err_z = 2.0 * ϵ_tol; err_n = 2.0 * ϵ_tol; iters = 0;
    while err_z > ϵ_tol && err_n > ϵ_tol && iters <= iter_max
        # Update counter
        iters += 1

        # Calculate xn using Zygote derivative
        xn_z -= f(xn_z) / f'(xn_z)
        push!(xn_z_all, xn_z)

        # Calculate xn using the numerical dervicative
        xn_p   = xn_n + xn_n * 1e-11
        df_dxn = (f(xn_n) - f(xn_p)) / (xn_n - xn_p)
        xn_n  -= f(xn_n) / df_dxn
        push!(xn_n_all, xn_n)

        # Calculate xn using the analytical derivative
        xn_a -= f(xn_a) / df(xn_a)
        push!(xn_a_all, xn_a)

        # Calculate errs
        err_z = abs(xn_z_all[iters+1] - xn_z_all[iters])
        err_n = abs(xn_n_all[iters+1] - xn_n_all[iters])
        err_a = abs(xn_a_all[iters+1] - xn_a_all[iters])

        # Visualize
        ln4 = lines!(ax1, xplot, tangent_line.(f, xplot, xn_z), color = (:gold, α[iters]), linewidth = 3, label = "f'(x)")
        sc2 = scatter!(ax1, xn_z, f(xn_z), color = (:purple4, α[iters]), markersize = 20.0, label = "Zyg.")
        sc3 = scatter!(ax1, xn_n, f(xn_n), color = (:purple4, α[iters]), markersize = 20.0, marker = :utriangle, label = "Num.")
        sc3 = scatter!(ax1, xn_a, f(xn_a), color = (:purple4, α[iters]), markersize = 20.0, marker = :rect, label = "Ana.")
        sc4 = scatter!(ax2, iters, err_z, label = "Zyg.", color = (:magenta1, 0.5), markersize = 20.0)
        sc4 = scatter!(ax2, iters, err_n, label = "Num.", color = (:cyan, 0.5), marker = :utriangle, markersize = 20.0)
        sc4 = scatter!(ax2, iters, err_a, label = "Ana.", color = (:purple4, 0.5), marker = :rect, markersize = 20.0)
        if iters == 1   # Add legends
            axislegend(ax1, position = :lt)
            axislegend(ax2, position = :lb)
        end
        display(fg1)

        # Save Figure
        save("/Users/lcandiot/Developer/DevelopersCard/doc/png/NewtonRaphson_$(iters).png", fg1, px_per_unit = 2)

        # Print progress to screen
        println("xn_z = $(xn_z); xn_n = $xn_n; xn_a = $xn_a ; err_z = $err_z; err_n = $err_n; err_a = $err_a; iters = $iters")
    end

    # Return
    return nothing
end

# Run
newtonRaphson()