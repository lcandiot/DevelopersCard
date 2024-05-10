# Illustrating the Newton-Raphson iteration scheme to find the minimum of a given f(x).
using Pkg
Pkg.activate(".")
using CairoMakie, MathTeXEngine, Zygote

function newtonRaphson()
    
    # Initial guess and function definition
    f(x) =  x.^2 .- 2.0
    df(x) = 2.0 * x
    tangent_line(f, x, x0) = f(x0) .+ f'(x0) .* (x .- x0)
    x0   = 10.0
    xn_z   = copy(x0)
    xn_n   = copy(x0)
    xn_a   = copy(x0)
    xn_z_all = [x0]
    xn_a_all = [x0]
    xn_n_all = [x0]

    # Solver settings
    iter_max = 1000
    ϵ_tol    = 1e-5

    # Visualize
    xplot = [0.0 + i*(x0 / 100) for i in 1:100]
    α = [0.1 + i*(1.0 / 11.0) for i in 1:11]
    fg1 = Figure(fontsize = 18)
    ax1 = Axis(fg1[1,1], xlabel = L"x_n", ylabel = L"f(x)", limits = (0.0, x0, -10.0, f(x0) + 2.0), aspect = 1)
    ax2 = Axis(fg1[1,2], xlabel = "iters", ylabel = "errs", yscale = log10, limits = (0, 7.0, ϵ_tol / 10.0, 100.0), aspect = 1)
    ln3 = lines!(ax1, xplot, (f.(xplot)), label = "f(x)", color = :steelblue, linewidth = 3)

    # Iteration loop
    err_a = 2.0 * ϵ_tol; err_z = 2.0 * ϵ_tol; err_n = 2.0 * ϵ_tol; iters = 0;
    while err_z > ϵ_tol && err_n > ϵ_tol && iters <= iter_max
        # Update counter
        iters += 1

        # Calculate xn Zygote
        xn_z -= f(xn_z) / f'(xn_z)
        push!(xn_z_all, xn_z)

        # Calculate xn numerically
        xn_p   = xn_n + xn_n * 1e-11
        df_dxn = (f(xn_n) - f(xn_p)) / (xn_n - xn_p)
        xn_n  -= f(xn_n) / df_dxn
        push!(xn_n_all, xn_n)

        # Calculate xn analytically
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
        sc4 = scatter!(ax2, iters, err_z, label = "Zyg.", color = (:magenta1, 0.5), marker = :circ, markersize = 20.0)
        sc4 = scatter!(ax2, iters, err_n, label = "Num.", color = (:cyan, 0.5), marker = :utriangle, markersize = 20.0)
        sc4 = scatter!(ax2, iters, err_a, label = "Ana.", color = (:purple4, 0.5), marker = :rect, markersize = 20.0)
        if iters == 1
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