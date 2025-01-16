# Solving 2D stokes flow of a compressible highly viscous fluid
using CairoMakie, Printf, ColorSchemes

macro av(A)  esc(:( ($A[1:end-1, 1:end-1] .+ $A[2:end, 1:end-1] .+ $A[1:end-1, 2:end] .+ $A[2:end, 2:end]) .* 0.25 )) end
macro av_xi(A) esc(:( ($A[1:end-1, 2:end-1] .+ $A[2:end,   2:end-1]) .* 0.5 )) end
macro av_yi(A) esc(:( ($A[2:end-1, 1:end-1] .+ $A[2:end-1, 2:end  ]) .* 0.5 )) end
macro av_xa(A) esc(:( ($A[1:end-1,       :] .+ $A[2:end,         :]) .* 0.5 )) end
macro av_ya(A) esc(:( ($A[:,       1:end-1] .+ $A[:,         2:end]) .* 0.5 )) end
macro ex_x(A)  esc(:( vcat($A[1, :]', $A) |> B -> vcat(B, $A[end, :]'))) end
macro ex_y(A)  esc(:( hcat($A[:, 1], $A) |> B -> hcat(B, $A[:, end]))) end

@views function compressible_viscous_stokes2D()
    # Real values
    ρ_r      = 2800.0                   # Density [kg.m⁻³]
    gy       = 9.81                     # Grav. acc. y [m.s⁻²]
    gx       = 0.0                      # Grav. acc. x [m.s⁻²]
    sec_year = 3600.0 * 24.0 * 365.25   # Seconds in a year
    ηs_inc   = 1e21                     # Viscosity of the inclusion [Pa.s]
    K_r      = 1e11                     # Bulk modulus [Pa]
    K_inc    = 1e9                     # Bulk modulus [Pa]
    p0       = 101300.0                 # Reference pressure [Pa]

    # Indepentent physics
    ρg       = ρ_r * gy                 # Buoyancy term
    ε̇_bg     = 1e-15                    # Background strain rate [s⁻¹]
    η_bg     = 1e23                     # Viscosity of the host rock [Pa.s]

    # Nondimensional numbers
    Lx_Ly    = 1.0                      # Domain ratio
    Lx_Lc    = 10.0                     # Ratio between char. length and domain length

    # Dependent scales
    Psc      = η_bg * ε̇_bg              # Pressure scale
    Lsc      = 1.0 / ( ρg / (η_bg * ε̇_bg) ) # Length scale
    tsc      = 1.0 / ε̇_bg               # Time scale
    vsc      = Lsc / tsc                # Velocity scale

    # Dependent physics
    dt       = 1e-3 * tsc               # Time step [s]
    Lx       = Lsc * Lx_Lc              # Domain length [m]
    Ly       = Lx / Lx_Ly               # Domain height [m]
    Lx_r     = 10.0                     # Ratio of domain length over inc. radius []
    ρgx      = ρ_r * gx                 # Buoyancy in x
    ρgy      = ρ_r * gy                 # Buoyancy in y
    max_LxLy = max(Lx, Ly)              # Maximum domain extend [m]


    # Numerics
    nx, ny   = 155, 155                 # Resolution
    nt       = 1                      # No. time steps
    ϵ_tol    = 1e-10                    # Absolute residual tolerance
    max_iter = 1e5                      # Maximum no. of iteration
    ncheck   = 10000                     # Error check frequency
    CFL      = 0.9 / sqrt(2.0)          # Stability criterion
    dx, dy   = Lx / nx, Ly / ny         # Grid spacing
    Ṽ        = min(dx,dy) * CFL         
    Re       = 3.0 * sqrt(10.0) / 2.0 * π # Reynolds number
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
    Gdτ      = zeros(Float64, nx    , ny    )
    dτ_ρ     = zeros(Float64, nx    , ny    )
    ∇v       = zeros(Float64, nx    , ny    )
    ηs       = zeros(Float64, nx    , ny    ) .+ η_bg
    K        = zeros(Float64, nx    , ny    ) .+ K_r
    ρ        = zeros(Float64, nx    , ny    ) .+ ρ_r
    ρ_old    = zeros(Float64, nx    , ny    ) .+ ρ_r
    p        = zeros(Float64, nx, ny) #-cumsum(ρ .* gy, dims = 2) .* dy
    p_old    = deepcopy(p)
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
        vy[idx, :] .=  ε̇_bg .* LinRange(-0.5*Ly, 0.5*Ly, ny+1)
    end

    # Set circular inclusion
    for idx_y in eachindex(y)
        for idx_x in eachindex(x)
            if ( (x2D[idx_x, idx_y] - x_inc)^2 + (y2D[idx_x, idx_y] - y_inc)^2 < r^2 )
                ηs[idx_x, idx_y] = ηs_inc
                K[idx_x, idx_y] = K_inc
            end
        end
    end

    smooth_2DArray_diffusion!(ηs, 10, dx, dy)

    ηs_ini = deepcopy(ηs)

    # Visualize
    fg1   = Figure(size = (1600, 1600))
    ax1   = Axis(fg1[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Pt")
    ax2   = Axis(fg1[2, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vx")
    ax3   = Axis(fg1[2, 2], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vy")
    ax4   = Axis(fg1[3, 1][1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "ρ")
    hm1   = heatmap!(ax1, x, y, p)
    hm2   = heatmap!(ax2, xv, yv, vx)
    hm3   = heatmap!(ax3, xv, yv, vy)
    hm4   = heatmap!(ax4, x, y, ρ)
    # display(fg1)

    # Time loop
    time = 0.0
    for idx_Time in 1:nt
        
        # Update time
        dt   = min(η_bg / K_r / 10.0, min(dx, dy) / maximum(max.(abs.(@av_xa(vx)), abs.(@av_ya(vy)))) / 4.1)
        time = time + dt

        # Calculate pseudo-transient params
        dτ_ρ .= Ṽ .* max_LxLy ./ Re ./ (1.0 ./ (1.0 ./ (K_r .* dt) .+ 1.0 ./ (ηs)))
        # Gdτ  .= Ṽ^2 ./ dτ_ρ ./ (K_r)
        Gdτ  .= Ṽ^2 ./ dτ_ρ ./ (r̃ + 2.0)

        # Update old
        p_old .= p
        ρ_old .= ρ

        # Pseudo - transient solver loop
        max_Res_all = []; err = 2ϵ_tol
        for iter in 1:max_iter

            # Calculate Density
            ρ .= ρ_r .* (1.0 .+ 1.0 ./ K .* (p .- p0))

            # Conservation of mass - pressure update
            ∇v    .= diff(vx, dims = 1) ./ dx .+ diff(vy, dims = 2) ./ dy
            Res_p .= - r̃ .* Gdτ .* (∇v .+ 1.0 ./ K .* (p .- p_old) ./ dt)
            p    .+= Res_p

            # Constitutive equation - stress update
            τxx   .= ( τxx .+ 2.0 .* Gdτ .* (diff(vx, dims = 1) ./ dx .- 0.3 .* (diff(vx, dims = 1) ./ dx + diff(vy, dims = 2) ./ dy)) ) ./ (Gdτ ./ ηs .+ 1.0)
            τyy   .= ( τyy .+ 2.0 .* Gdτ .* (diff(vy, dims = 2) ./ dy .- 0.3 .* (diff(vx, dims = 1) ./ dx + diff(vy, dims = 2) ./ dy)) ) ./ (Gdτ ./ ηs .+ 1.0)
            τxy   .= ( τxy .+ 2.0 .* @av(Gdτ) .* (0.5 .* (diff(vx[2:end-1,:], dims = 2) ./ dy + diff(vy[:, 2:end-1], dims = 1) ./ dx)) ) ./ (@av(Gdτ) ./ @av(ηs) .+ 1.0)

            # Conservation of linear momentum - velocity update
            Res_vx .= @av_xi(dτ_ρ) .* ( diff(τxx[:, 2:end-1], dims = 1) ./ dx .+ diff(τxy, dims = 2) ./ dy .- diff(p[:, 2:end-1], dims = 1) ./ dx ) 
            Res_vy .= @av_yi(dτ_ρ) .* ( diff(τyy[2:end-1, :], dims = 2) ./ dy .+ diff(τxy, dims = 1) ./ dx .- diff(p[2:end-1, :], dims = 2) ./ dy )
            vx[2:end-1, 2:end-1]    .+= Res_vx
            vy[2:end-1, 2:end-1]    .+= Res_vy

            # Boundary conditions
            vx[:,   1] .= 0.0 # vx[:,       2]
            vx[:, end] .= 0.0 # vx[:, end - 1]
            vy[1,   :] .= 0.0 # vy[2,       :]
            vy[end, :] .= 0.0 # vy[end - 1, :]

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
        empty!(ax4)
        hm1 = heatmap!(ax1, x, y, p, colormap = :viridis)
        hm2 = heatmap!(ax2, xv, yv, vx, colormap = :roma)
        hm3 = heatmap!(ax3, xv, yv, vy, colormap = :roma)
        hm4 = contourf!(ax4, x, y, ρ, colormap = :roma)
        # ct1 = contour!(ax4, x, y, ηs_ini)

        Colorbar(fg1[1,2][1,1], hm1, label = "P [Pa]")
        Colorbar(fg1[1,2][1,2], hm2, label = "Vx [m.s⁻¹]")
        Colorbar(fg1[3,1][1,2], hm4, label = "ρ [kg.m-3]")

        display(fg1)
    end

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
compressible_viscous_stokes2D()

