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
    nx       = 501                  # No. of grid points
    nt       = 2000                  # No. of time steps
    nEig     = 70                   # No. of Eigenvalues
    nviz     = 100                   # Visualization increment
    cfl      = 1.0 + sqrt(2.1)      # CFL criterion
    ϵtol     = 1.0e-6               # Solver tolerance
    max_iter = 1000 * nx^2          # Iteration cap
    ncheck   = 1                    # Convergence check increment
    plot_res = true                 # Residual plotting switch

    # Initialisation
    dx     = Lx / (nx - 1)           # Grid spacing [m]
    T      = Vector{Float64}(undef, nx - 1)
    dTdt   = Vector{Float64}(undef, nx - 1)
    qT     = Vector{Float64}(undef, nx    )
    sumEig = Vector{Float64}(undef, nx - 1)
    ζ_n    = Vector{Float64}(undef, nEig  )
    x      = [0.0 - dx / 2.0 + ix * dx for ix = 1:nx-1]
    T     .= TA .* exp.( .-( (x .- xc) / w) .^2)
    Ti     = copy(T)
    Ta     = copy(Ti)
    T_old  = copy(Ti)
    dTdt  .= 0.0
    qT    .= 0.0

    # Compute coefficients for analytical solution - One has to take the integral between 0, Lx, which is why we don't use x array
    for iEig in eachindex(ζ_n)
        int, int_err = quadgk(xint -> exp(-(v / 2.0 / D) * xint) * TA * exp( -((xint - xc) / w) ^2) * sin((iEig * π / Lx) * xint), 0, Lx, rtol= ϵtol)
        ζ_n[iEig] = 2.0 / Lx * int
    end

    # Visualize initial configuration
    f     = Figure()
    ax1   = Axis(f[1,1], title="Temperature at time = $(time)")
    ax2   = Axis(f[2,1], yscale = log10)
    lines!(ax1, x, Ti)
    display(f)

    # Time loop
    for iTime = 1:nt
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
        err = 1.0; iter = 0
        while err > ϵtol && iter < max_iter
            iter += 1
            # Calculate heat flux
            for iX in eachindex(T)
                if iX < nx - 1
                    qT[iX + 1, 1] -= (qT[iX + 1, 1] + D * (T[iX + 1, 1] - T[iX, 1]) / dx) * 1.0 / (1.0 + θ_dτ_T)
                end
            end
            # qT[end, 1] = 0.0# qT[end - 1, 1] # BCs
            # qT[1  , 1] = 0.0# qT[end    , 1]

            # Calculate temperature change 
            for iX in eachindex(dTdt)
                dTdt[iX, 1] = (qT[iX + 1, 1] - qT[iX, 1]) / dx             # Diffusion
            end

            for iX in eachindex(dTdt)
                if iX > 1 && iX < nx - 1 
                    dTdt[iX  , 1] += max(0.0, v ) * (T[iX, 1] - T[iX - 1, 1]) / dx     # Advection (upwind)
                    dTdt[iX-1, 1] += min(v , 0.0) * (T[iX, 1] - T[iX - 1, 1]) / dx
                end
            end

            # Update temperature
            for iX in eachindex(T)
                T[iX, 1] -= 1.0 / (1.0 / dt + β_dτ_T) * ((T[iX, 1] - T_old[iX, 1]) / dt + dTdt[iX, 1])
            end

            # Dirichlet BC
            T[[1 end], 1] .= 0.0

            # Check the error
            if iter % ncheck == 0
                err = maximum(abs.(-(T[2:end-1, 1] .- T_old[2:end-1, 1]) ./ dt .- dTdt[2:end-1, 1] )) * dt / TA
                println("err = $(err); iter = $(iter)")
                println("min(T ) = $(minimum(T));  max(T ) = $(maximum(T))")
                println("min(qT) = $(minimum(qT)); max(qT) = $(maximum(qT))")
                if plot_res && iTime == nt
                    scatter!(ax2, iter, err, color=:black)
                    ax2.title = "Residual at $(iter) iterations"
                    display(f)
                end
            end
        end

        # Calculate analytical solution
        sumEig .= 0.0
        for iEig in eachindex(ζ_n)
            sumEig .+= ζ_n[iEig] .* sin.((iEig * π / Lx) .* x) .* exp(-D * (iEig * π / Lx)^2 * time)
        end
        Ta .= exp(-(v^2 / 4.0 / D) * time) .* exp.((v / 2.0 / D) .* x) .* sumEig

        # Visualize
        if iTime % nviz == 0
            empty!(ax1)
            li1 = lines!(ax1, x, Ti, color=:tomato, label = "Initial")
            li2 = lines!(ax1, x, T , color=:blue, label = "Num. solution")
            sc1 = scatter!(ax1, x[1:10:end, 1], Ta[1:10:end, 1], color=:blue, label = "Ana. solution")
            ax1.title = "Temperature at time = $(time)"
            axislegend(ax1, position = :rt)
            display(f)

        end

        # Save final result
        if iTime == nt
            save("./doc/png/DiffusionConvection_1D_homoDirichlet.png", f, px_per_unit=3)
        end

        # Clear plot
        empty!(ax2)
    end

    # Return
    return nothing

end

DiffusionConvection_1D()