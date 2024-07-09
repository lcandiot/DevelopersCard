# Solving 2D stokes flow of an incompressible highly viscous fluid
using CairoMakie

@views function incompressible_viscous_stokes2D()
    # Real values
    ρ_r    = 1.0
    gy     = 1.0
    gx     = 0.0
    yr     = 3600.0 * 24.0 * 365.25
    ηs_inc = 1e-3
    Ṽ      = 6000.0

    # Indepentent physics
    ρg   = ρ_r * gy
    ε̇_bg = 1.0
    η_bg = 1.0

    # Nondimensional number
    Lx_Ly = 1.0
    Lx_Lc = 10.0

    # Dependent physics
    Lc   = 1.0 / ( ρg / (η_bg * ε̇_bg) )
    tc   = 1.0 / ε̇_bg
    dt   = 1e-3 * tc
    Lx   = Lc * Lx_Lc
    Ly   = Lx / Lx_Ly
    Lx_r = 10.0
    ρgx  = ρ_r * gx

    # Dependent scales
    Psc = η_bg * ε̇_bg
    vsc = Lc   * ε̇_bg
    tsc = 1.0  / ε̇_bg

    # Numerics
    nx, ny   = 151, 151
    ϵ_tol    = 1e-6
    max_iter = 1000
    dx, dy   = Lx / (nx - 1), Ly / (ny - 1)
    Re       = 3 * √10 / 2.0 * π
    r̃        = 0.5
    ρ̃        = Re * ηs_bg / Ṽ / max(Lx, Ly)
    G̃        = ρ̃ * Ṽ^2 / (r̃ + 2.0)
    K̃        = r̃ * G̃
    
    # Mechanics switches
    is_buoyant = 0

    # Initialization
    vx     = zeros(Float64, nx + 1, ny    )
    vy     = zeros(Float64, nx    , ny + 1)
    τxx    = zeros(Float64, nx    , ny    )
    τxy    = zeros(Float64, nx - 1, ny - 1)
    Res_vx = zeros(Float64, nx - 1, ny - 2)
    Res_vy = zeros(Float64, nx - 2, ny - 1)
    p      = zeros(Float64, nx    , ny    )
    ηs     = zeros(Float64, nx    , ny    ) .+ η_bg
    x      = LinRange(0, Lx, nx)
    y      = LinRange(0, Ly, ny)
    xc     = Lx / 2.0
    yc     = Ly / 2.0
    r      = Lx / Lx_r
    x2D    = reshape(repeat(x , ny), nx, ny)
    y2D    = reshape(repeat(y', nx), nx, ny)

    for idx_y in eachindex(y)
        for idx_x in eachindex(x)
            if ( (x2D[idx_x, idx_y] - xc)^2 + (y2D[idx_x, idx_y] - yc)^2 < r^2 )
                ηs[idx_x, idx_y] = ηs_inc
            end
        end
    end

    smooth_2DArray_diffusion(ηs, 10, dx, dy)

    # Visualize
    fg1 = Figure(size = (600, 600))
    ax1 = Axis(fg1[1, 1], xlabel = "x", ylabel = "y")
    hm1 = heatmap!(ax1, x, y, ηs)
    display(fg1)

    # Pseudo - transient solver loop
    for iter in 1:max_iter
        
        # Conservation of mass - pressure update
        Res_p = - K̃ * dτ * (diff(vx, dims = 1) / dx + diff(vy, dims = 2) / dy)
        p   .+= dτ * Res_p

        # Constitutive equation - stress update
        τxx   .-= 2.0 .* G̃ .* dτ .* (τxx ./ 2.0 ./ ηs .- diff(vx, dims = 1) ./ dx)
        τyy   .-= 2.0 .* G̃ .* dτ .* (τyy ./ 2.0 ./ ηs .- diff(vy, dims = 2) ./ dy)
        τxy   .-= 2.0 .* G̃ .* dτ .* (τxy ./ 2.0 ./ ηs .- 0.5 .* (diff(vx, dims = 2) ./ dy + diff(vy, dims = 1) ./ dx))

        # Conservation of linear momentum - velocity update
        Res_vx .= - 1.0 ./ ρ̃ .* ( diff(p, dims = 1) / dx - diff(τxx, dims = 1) / dx - diff(τxy, dims = 2) / dy) + Float64(is_buoyant) * ρgx 
        Res_vy .= - 1.0 ./ ρ̃ .* ( diff(p, dims = 2) / dy - diff(τyy, dims = 2) / dy - diff(τxy, dims = 1) / dx) + Float64(is_buoyant) * ρgy
        vx[2:end-1, 2:end-1]    .+= dτ .* Res_vx 
        vy[2:end-1, 2:end-1]    .+= dτ .* Res_vy 

    end

    # Return
    return nothing
end

# Diffuse array
function smooth_2DArray_diffusion(A :: Matrix{Float64}, nsteps :: Int64, dx :: Float64, dy :: Float64)
    for _ in 1:nsteps
        A[2:end-1, 2:end-1] .+= (min.(dx.^2, dy.^2) ./ 1.0 ./ 4.1) .* ( diff( 1.0 .* diff(A[:, 2:end-1], dims = 1), dims = 1) ./ dx^2 + diff(1.0 .* diff(A[2:end-1, :], dims = 2), dims = 2) ./ dy^2)
    end

    # Return
    return nothing
end

# Run main
incompressible_viscous_stokes2D()

