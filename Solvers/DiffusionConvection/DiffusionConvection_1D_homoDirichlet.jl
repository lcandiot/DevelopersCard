# Solving 1D diffusion-convection equation numerically and comparing results to analytical solutions. Initial condition is a Gaussian and homogeneous Dirichlet BC are used.
using CairoMakie, QuadGK

# main
@views function DiffusionConvection_1D()
    # Physics
    Lx   = 10.0             # Length [m]
    xc   = Lx / 2.0         # Location of amplitude [m]
    D    = 1.0              # Diffusion coefficient [m^2/s]
    time = 0.0              # Physical time [s]
    v    = 1.5              # Velocity solid [m/s]
    w    = 0.5              # Perturbation width [m]
    TA   = 1.0              # Perturbation amplitude [K]

    # Numerics
    nx        = 501                  # No. of grid points
    nt        = 2000                 # No. of time steps
    nEig      = 70                   # No. of Eigenvalues
    nviz      = 100                  # Visualization increment
    cfl       = 1.0 + sqrt(2.1)      # CFL criterion
    ϵtol      = 1.0e-6               # Solver tolerance
    max_iter  = 1000 * nx^2          # Iteration cap
    ncheck    = 1                    # Convergence check increment
    plot_res  = true                 # Residual plotting switch
    plot_leg  = true                 # Legend plotting switch
    print_fig = true
    
    # Derived numerics
    dx     = Lx / (nx - 1)                                  # Grid spacing [m]

    # Initialisation
    T      = Vector{Float64}(undef, nx - 1)                 # Temperature [K]
    dTdt   = Vector{Float64}(undef, nx - 1)                 # Temperature time derivative [K/s]
    qT     = Vector{Float64}(undef, nx    )                 # Heat flux (W/m2)
    sumEig = Vector{Float64}(undef, nx - 1)                 # Sum of Eigenvalues
    ζ_n    = Vector{Float64}(undef, nEig  )                 # Eigenvalues
    x      = [0.0 - dx / 2.0 + ix * dx for ix = 1:nx-1]     # Coordinate array
    T     .= TA .* exp.( .-( (x .- xc) / w) .^2)            # Initial condition
    Ti     = copy(T)
    Ta     = copy(Ti)                                       # Analytical T array
    T_old  = copy(Ti)                                       # Old physical T array
    dTdt  .= 0.0
    qT    .= 0.0

    # Compute coefficients for analytical solution - One has to take the integral between 0, Lx, which is why we don't use x array
    for iEig in eachindex(ζ_n)
        int, int_err = quadgk(xint -> exp(-(v / 2.0 / D) * xint) * TA * exp( -((xint - xc) / w) ^2) * sin((iEig * π / Lx) * xint), 0, Lx, rtol= ϵtol)
        ζ_n[iEig] = 2.0 / Lx * int
    end

    # Visualize initial configuration
    f     = Figure(fontsize = 16)
    ax1   = Axis(f[1,1], title="Temperature at time = $(time)", xlabel = "x [m]", ylabel = "T [K]")
    ax2   = Axis(f[2,1], yscale = log10, xlabel = "i []", ylabel = "err []")
    lines!(ax1, x, Ti)
    display(f)

    # Time loop
    for iTime = 1:nt
        # Time stepping
        dt = if iTime == 1
            min(dx^2 / D / 2.1, dx / abs(v) / 2.1)
        else
            min(dx^2 / D / 2.1, dx / abs(v) / 2.1)
        end

        # Update convergence parameters
        re_T    = π + sqrt(π^2 + Lx^2 / D / dt)
        θ_dτ_T  = Lx / re_T / cfl / dx
        β_dτ_T  = (re_T * D) / (cfl * dx * Lx)

        # Update physical parameters
        time += dt
        T_old .= T

        # Iteration loop
        err = 1.0; iter = 0; iter_all = []; err_all = []
        while err > ϵtol && iter < max_iter
            iter += 1
            # Calculate heat flux
            for iX in eachindex(T)
                if iX < nx - 1
                    qT[iX + 1, 1] -= (qT[iX + 1, 1] + D * (T[iX + 1, 1] - T[iX, 1]) / dx) * 1.0 / (1.0 + θ_dτ_T)
                end
            end

            # Calculate temperature change 
            for iX in eachindex(dTdt)
                dTdt[iX, 1] = (qT[iX + 1, 1] - qT[iX, 1]) / dx                          # Diffusion
            end

            for iX in eachindex(dTdt)
                if iX > 1 && iX < nx - 1 
                    dTdt[iX  , 1] += max(0.0, v ) * (T[iX, 1] - T[iX - 1, 1]) / dx      # Advection (upwind)
                    dTdt[iX-1, 1] += min(v , 0.0) * (T[iX, 1] - T[iX - 1, 1]) / dx
                end
            end

            # Update temperature
            for iX in eachindex(T)
                T[iX, 1] -= 1.0 / (1.0 / dt + β_dτ_T) * ((T[iX, 1] - T_old[iX, 1]) / dt + dTdt[iX, 1])
            end

            # Homogeneous Dirichlet BC
            T[[1 end], 1] .= 0.0

            # Check the error
            if iter % ncheck == 0
                err = maximum(abs.(-(T[2:end-1, 1] .- T_old[2:end-1, 1]) ./ dt .- dTdt[2:end-1, 1] )) * dt / TA
                push!(iter_all, iter)
                push!(err_all, err)
                println("err = $(err); iter = $(iter)")
                println("min(T ) = $(minimum(T));  max(T ) = $(maximum(T))")
                println("min(qT) = $(minimum(qT)); max(qT) = $(maximum(qT))")
            end
        end

        # Calculate analytical solution
        sumEig      .= 0.0                                                                              # Reset
        for iEig in eachindex(ζ_n)
            sumEig .+= ζ_n[iEig] .* sin.((iEig * π / Lx) .* x) .* exp(-D * (iEig * π / Lx)^2 * time)    # Calculate sum of Eigenvalues
        end
        Ta          .= exp(-(v^2 / 4.0 / D) * time) .* exp.((v / 2.0 / D) .* x) .* sumEig               # Compute analytical solution

        # Visualize
        if iTime % nviz == 0
            empty!(ax1)
            li1 = lines!(ax1, x, Ti, color=:purple4, linewidth = 2, label = "Initial")
            li2 = lines!(ax1, x, T , color=:steelblue, linewidth = 2, label = "Num. solution")
            sc1 = scatter!(ax1, x[1:15:end, 1], Ta[1:15:end, 1], color=:steelblue, markersize = 10.0, label = "Ana. solution")
            ax1.title = "Temperature at time = $(time)"
            if plot_leg
                axislegend(ax1, position = :rt)
            end
            plot_leg = false
            if plot_res
                scatter!(ax2, Float64.(iter_all), Float64.(err_all), color=:black)
                ax2.title = "Residual at $(iter) iterations"
            end
            display(f)
            if print_fig
                # Save results
                save("/Users/lcandiot/Developer/DevelopersCard/doc/png/DiffusionConvection/DiffusionConvection_1D_homoDirichlet_$(iTime).png", f, px_per_unit=3)
            end
        end

        # Clear plot
        empty!(ax2)
    end

    # Return
    return nothing

end

DiffusionConvection_1D()