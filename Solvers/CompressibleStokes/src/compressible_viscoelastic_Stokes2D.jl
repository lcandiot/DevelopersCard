# Solving 2D stokes flow of an incompressible highly viscous fluid
using CairoMakie, Printf, ColorSchemes

macro av(A)  esc(:( ($A[1:end-1, 1:end-1] .+ $A[2:end, 1:end-1] .+ $A[1:end-1, 2:end] .+ $A[2:end, 2:end]) .* 0.25 )) end
macro av_xi(A) esc(:( ($A[1:end-1, 2:end-1] .+ $A[2:end, 2:end-1]) .* 0.5 )) end
macro av_yi(A) esc(:( ($A[2:end-1, 1:end-1] .+ $A[2:end-1, 2:end]) .* 0.5 )) end

@views function compressible_viscoelastic_stokes2D()
    # Real values
    ρ_r      = 2800.0
    gy       = 9.81
    gx       = 0.0
    sec_year = 3600.0 * 24.0 * 365.25
    ηs_inc   = 1e20
    T0       = 500.0 + 273.15
    G_r      = 1e10

    # Indepentent physics
    ρg       = ρ_r * gy
    ε̇_bg     = 1e-15
    η_bg     = 1e23
    TA       = 200.0

    # Nondimensional number
    Lx_Ly    = 1.0
    Lx_Lc    = 10.0

    # Dependent scales
    Psc      = η_bg * ε̇_bg
    Lsc      = 1.0 / ( ρg / (η_bg * ε̇_bg) )
    tsc      = 1.0 / ε̇_bg
    vsc      = Lsc / tsc

    # Dependent physics
    dt       = 1e-3 * tsc
    Lx       = Lsc * Lx_Lc
    Ly       = Lx / Lx_Ly
    Lx_r     = 10.0
    ρgx      = ρ_r * gx
    ρgy      = ρ_r * gy
    max_LxLy = max(Lx, Ly)


    # Numerics
    nx, ny   = 255, 255
    ϵ_tol    = 1e-10
    max_iter = 1e5
    ncheck   = 1000
    CFL      = 0.9 / sqrt(2.0)
    dx, dy   = Lx / nx, Ly / ny
    Ṽ        = min(dx,dy) * CFL
    Re       = 3.0 * sqrt(10.0) / 2.0 * π
    r̃        = 0.5

    # Mechanics switches
    is_buoyant = 0

    # Initialization
    vx       = zeros(Float64, nx + 1, ny    )
    vy       = zeros(Float64, nx    , ny + 1)
    τxx      = zeros(Float64, nx    , ny    )
    τyy      = zeros(Float64, nx    , ny    )
    τxy      = zeros(Float64, nx - 1, ny - 1)
    Res_vx   = zeros(Float64, nx - 1, ny - 2)
    Res_vy   = zeros(Float64, nx - 2, ny - 1)
    Res_p    = zeros(Float64, nx    , ny    )
    G̃dτ      = zeros(Float64, nx    , ny    )
    dτ_ρ     = zeros(Float64, nx    , ny    )
    p        = zeros(Float64, nx    , ny    )
    ηs       = zeros(Float64, nx    , ny    ) .+ η_bg
    ηve      = zeros(Float64, nx    , ny    )
    G        = zeros(Float64, nx    , ny    ) .+ G_r
    x        = LinRange(0 + dx / 2.0, Lx - dx / 2.0, nx)
    y        = LinRange(0 + dy / 2.0, Ly - dy / 2.0, ny)
    xv       = LinRange(0, Lx, nx + 1)
    yv       = LinRange(0, Ly, ny + 1)
    x_inc    = Lx / 2.0
    y_inc    = Ly / 2.0
    r        = Lx / Lx_r
    x2D      = reshape(repeat(x , ny), nx, ny)
    y2D      = reshape(repeat(y', nx), nx, ny)

    # Boundary conditions
    for idy in axes(vx, 2)
        vx[:, idy] .= -ε̇_bg .* LinRange(-0.5*Lx, 0.5*Lx, nx+1)
    end
    for idx in axes(vy, 1)
        vy[idx, :] .= ε̇_bg .* LinRange(-0.5*Ly, 0.5*Ly, ny+1)
    end

    # Set circular inclusion
    for idx_y in eachindex(y)
        for idx_x in eachindex(x)
            if ( (x2D[idx_x, idx_y] - x_inc)^2 + (y2D[idx_x, idx_y] - y_inc)^2 < r^2 )
                ηs[idx_x, idx_y] = ηs_inc
            end
        end
    end

    smooth_2DArray_diffusion!(ηs, 10, dx, dy)

    # Visualize
    fg1   = Figure(size = (600, 600))
    ax1   = Axis(fg1[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Pt")
    ax2   = Axis(fg1[2, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vx")
    ax3   = Axis(fg1[2, 2], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vy")
    hm1   = heatmap!(ax1, x, y, ηs)
    hm2   = heatmap!(ax2, xv, yv, vx)
    hm3   = heatmap!(ax3, xv, yv, vy)
    display(fg1)

    # Calculate pseudo-transient params
    dτ_ρ .= Ṽ .* max_LxLy ./ Re ./ ηs
    G̃dτ  .= Ṽ^2 ./ dτ_ρ ./ (r̃ + 2.0)

    # Pseudo - transient solver loop
    max_Res_all = []; err = 2ϵ_tol
    for iter in 1:max_iter

        # Conservation of mass - pressure update
        Res_p .= - r̃ .* G̃dτ .* (diff(vx, dims = 1) ./ dx .+ diff(vy, dims = 2) ./ dy)
        p    .+= Res_p

        # Constitutive equation - stress update
        τxx   .= ( τxx .+ 2.0 .* G̃dτ .* diff(vx, dims = 1) ./ dx ) ./ (G̃dτ ./ ηs .+ 1.0)
        τyy   .= ( τyy .+ 2.0 .* G̃dτ .* diff(vy, dims = 2) ./ dy ) ./ (G̃dτ ./ ηs .+ 1.0)
        τxy   .= ( τxy .+ 2.0 .* @av(G̃dτ) .* (0.5 .* (diff(vx[2:end-1,:], dims = 2) ./ dy + diff(vy[:, 2:end-1], dims = 1) ./ dx)) ) ./ (@av(G̃dτ) ./ @av(ηs) .+ 1.0)

        # Conservation of linear momentum - velocity update
        Res_vx .= @av_xi(dτ_ρ) .* ( diff(τxx[:, 2:end-1], dims = 1) ./ dx .+ diff(τxy, dims = 2) ./ dy .- diff(p[:, 2:end-1], dims = 1) ./ dx ) 
        Res_vy .= @av_yi(dτ_ρ) .* ( diff(τyy[2:end-1, :], dims = 2) ./ dy .+ diff(τxy, dims = 1) ./ dx .- diff(p[2:end-1, :], dims = 2) ./ dy )
        vx[2:end-1, 2:end-1]    .+= Res_vx
        vy[2:end-1, 2:end-1]    .+= Res_vy

        # Boundary conditions
        vx[:,   1] .= vx[:,       2]
        vx[:, end] .= vx[:, end - 1]
        vy[1,   :] .= vy[2,       :]
        vy[end, :] .= vy[end - 1, :]

        # Monitor residuals 
        if iter % ncheck == 0
            @printf("\n")
            @printf("Iteration = %d \n", iter)
            @printf("Residual p  = %.2e \n", maximum(Res_p ) / Psc)
            @printf("Residual vx = %.2e \n", maximum(Res_vx) / vsc)
            @printf("Residual vy = %.2e \n", maximum(Res_vy) / vsc)
            err = maximum([maximum(Res_p ) / Psc maximum(Res_vx) / vsc maximum(Res_vy) / vsc])
        end

        if err < ϵ_tol
            break
        end
    end

    # Visualize
    empty!(ax1)
    empty!(ax2)
    empty!(ax3)
    hm1 = heatmap!(ax1, x, y, p, colormap = :viridis)
    hm2 = heatmap!(ax2, xv, yv, vx, colormap = :roma)
    hm3 = heatmap!(ax3, xv, yv, vy, colormap = :roma)

    Colorbar(fg1[1,2][1,1], hm1, label = "P [Pa]")
    Colorbar(fg1[1,2][1,2], hm2, label = "Vx [m.s⁻¹]")

    display(fg1)

    # Return
    return nothing
end

# Diffuse array
function smooth_2DArray_diffusion!(A :: Matrix{Float64}, nsteps :: Int64, dx :: Float64, dy :: Float64)
    for _ in 1:nsteps
        A[2:end-1, 2:end-1] .+= (min.(dx.^2, dy.^2) ./ 1.0 ./ 4.1) .* ( diff( 1.0 .* diff(A[:, 2:end-1], dims = 1), dims = 1) ./ dx^2 + diff(1.0 .* diff(A[2:end-1, :], dims = 2), dims = 2) ./ dy^2)
    end

    # Return
    return nothing
end

# Run main
compressible_viscoelastic_stokes2D();

