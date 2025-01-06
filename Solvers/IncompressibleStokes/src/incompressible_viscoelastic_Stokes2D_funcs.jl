# Solving 2D stokes flow of an incompressible viscoelastic fluid putting the building blocks into functions for speeding things up
using CairoMakie, Printf, ColorSchemes

macro av(A)    esc(:( ($A[1:end-1, 1:end-1] .+ $A[2:end, 1:end-1] .+ $A[1:end-1, 2:end] .+ $A[2:end, 2:end]) .* 0.25 )) end
macro av_xi(A) esc(:( ($A[1:end-1, 2:end-1] .+ $A[2:end, 2:end-1]) .* 0.5 )) end
macro av_yi(A) esc(:( ($A[2:end-1, 1:end-1] .+ $A[2:end-1, 2:end]) .* 0.5 )) end
macro av_xa(A) esc(:( ($A[1:end-1, :] .+ $A[2:end, :]) .* 0.5 )) end
macro av_ya(A) esc(:( ($A[:, 1:end-1] .+ $A[:, 2:end]) .* 0.5 )) end

@views function incompressible_viscoelastic_stokes2D()
    # Real values
    ρ_r      = 2800.0
    gy       = 9.81
    gx       = 0.0
    sec_year = 3600.0 * 24.0 * 365.25
    ηs_inc   = 1e20
    G_r      = 1e10

    # Indepentent physics
    ρg       = ρ_r * gy
    ε̇_bg     = 1e-15
    η_bg     = 1e23

    # Nondimensional number
    Lx_Ly    = 1.0
    Lx_Lc    = 10.0

    # Dependent scales
    Psc      = η_bg * ε̇_bg
    Lsc      = 1.0 / ( ρg / (η_bg * ε̇_bg) )
    tsc      = 1.0 / ε̇_bg
    vsc      = Lsc / tsc

    # Dependent physics
    Lx       = Lsc * Lx_Lc
    Ly       = Lx / Lx_Ly
    Lx_r     = 10.0
    ρgx      = ρ_r * gx
    ρgy      = ρ_r * gy
    max_LxLy = max(Lx, Ly)

    # Numerics
    nx, ny   = 255, 255
    nt       = 15
    ϵ_tol    = 1e-6
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
    τII      = zeros(Float64, nx - 1, ny - 1)
    Jxx      = zeros(Float64, nx    , ny    )
    Jyy      = zeros(Float64, nx    , ny    )
    Jxy      = zeros(Float64, nx - 1, ny - 1)
    Jyx      = zeros(Float64, nx - 1, ny - 1)
    ωxx      = zeros(Float64, nx    , ny    )
    ωyy      = zeros(Float64, nx    , ny    )
    ωxy      = zeros(Float64, nx - 1, ny - 1)
    ωyx      = zeros(Float64, nx - 1, ny - 1)
    Res_vx   = zeros(Float64, nx - 1, ny - 2)
    Res_vy   = zeros(Float64, nx - 2, ny - 1)
    Res_p    = zeros(Float64, nx    , ny    )
    G̃dτ      = zeros(Float64, nx    , ny    )
    dτ_ρ     = zeros(Float64, nx    , ny    )
    p        = zeros(Float64, nx    , ny    )
    ηs       = zeros(Float64, nx    , ny    ) .+ η_bg
    ηve      = zeros(Float64, nx    , ny    )
    G        = zeros(Float64, nx    , ny    ) .+ G_r
    xc       = LinRange(0 + dx / 2.0, Lx - dx / 2.0, nx)
    yc       = LinRange(0 + dy / 2.0, Ly - dy / 2.0, ny)
    xv       = LinRange(0, Lx, nx + 1)
    yv       = LinRange(0, Ly, ny + 1)
    x_inc    = Lx / 2.0
    y_inc    = Ly / 2.0
    r        = Lx / Lx_r
    x2D      = reshape(repeat(xc , ny), nx, ny)
    y2D      = reshape(repeat(yc', nx), nx, ny)
    # Old fields
    τxx_old  = deepcopy(τxx)
    τyy_old  = deepcopy(τyy)
    τxy_old  = deepcopy(τxy)

    # Boundary conditions
    for idy in axes(vx, 2)
        vx[:, idy] .= -ε̇_bg .* LinRange(-0.5*Lx, 0.5*Lx, nx+1)
    end
    for idx in axes(vy, 1)
        vy[idx, :] .= ε̇_bg .* LinRange(-0.5*Ly, 0.5*Ly, ny+1)
    end

    # Set circular inclusion
    for idx_y in eachindex(yc)
        for idx_x in eachindex(xc)
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
    ax4   = Axis(fg1[1, 2], xlabel = "nt", ylabel = "τII", title = "Stress build up")
    hm1   = heatmap!(ax1, xc, yc, ηs)
    hm2   = heatmap!(ax2, xv, yv, vx)
    hm3   = heatmap!(ax3, xv, yv, vy)
    display(fg1)

    # TIME LOOP
    time = 0.0; τII_max_time = []; ntimes = []
    for idx_Time in 1:nt

        # Update time
        dt   = min(η_bg / G_r / 10.0, min(dx, dy) / maximum(max.(abs.(@av_xa(vx)), abs.(@av_ya(vy)))) / 4.1)
        time = time + dt

        @printf("\ndt = %.5e\n", dt)

        # Update old physical fields
        τxx_old .= τxx
        τyy_old .= τyy
        τxy_old .= τxy

        # Pseudo - transient solver loop
        max_Res_all = []; err = 2ϵ_tol
        for iter in 1:max_iter

            # Evalutate viscoelastic viscosity
            compute_viscoelastic_rheology!(ηve, ηs, G, dt)

            # Calculate pseudo-transient params
            compute_pseudotransient_parameters!(ηve, dτ_ρ, G̃dτ, Ṽ, max_LxLy, Re, r̃)

            # Conservation of mass - pressure update
            compute_pressure_residual!(vx, vy, Res_p, G̃dτ, dx, dy, r̃)
            update_pressure!(p, Res_p)
    
            # Compute stresses
            compute_vorticity!(ωxx, ωyy, ωxy, ωyx, vx, vy, dx, dy)
            compute_Jaumann_derivative!(Jxx, Jyy, Jxy, ωxx, ωyy, ωxy, ωyx, τxx, τyy, τxy, τxx_old, τyy_old, τxy_old, vx, vy, dt, dx, dy)
            update_stresses!(τxx, τyy, τxy, Jxx, Jyy, Jxy, vx, vy, ηve, G,  G̃dτ, dx, dy)

            # Second invariant
            τII .= sqrt.( 0.5 .* (@av(τxx).^2 .+ @av(τyy).^2 .+ 2.0 .* τxy.^2) )

            # Conservation of linear momentum - velocity update
            compute_velocity_residual!(τxx, τyy, τxy, p, Res_vx, Res_vy, dτ_ρ, dx, dy)
            update_velocity!(vx, vy, Res_vx, Res_vy)
    
            # Boundary conditions
            set_BC_velocity!(vx, vy, Lx, Ly, ε̇_bg, "ε̇_bg const.")
    
            # Monitor residuals 
            if iter % ncheck == 0
                @printf("\n")
                @printf("Iteration = %d \n", iter)
                @printf("Residual p  = %.2e \n", maximum(Res_p ) / Psc)
                @printf("Residual vx = %.2e \n", maximum(Res_vx) / vsc)
                @printf("Residual vy = %.2e \n", maximum(Res_vy) / vsc)
                @printf("Physics:")
                @printf("min(vx) = %.2e; max(vx) = %.2e; min(vy) = %.2e; max(vy) = %.2e\n", minimum(vx), maximum(vx), minimum(vy), maximum(vy))
                err = maximum([maximum(Res_p ) / Psc maximum(Res_vx) / vsc maximum(Res_vy) / vsc])
            end
    
            if err < ϵ_tol
                break
            end
        end
    
        push!(τII_max_time, maximum(τII))
        push!(ntimes, idx_Time)

        println(τII_max_time)
        # Visualize
        empty!(ax1)
        empty!(ax2)
        empty!(ax3)
        empty!(ax4)
        hm1 = heatmap!(ax1, xc, yc, p, colormap = :viridis)
        hm2 = heatmap!(ax2, xv, yv, vx, colormap = :roma)
        hm3 = heatmap!(ax3, xv, yv, vy, colormap = :roma)
        ln  = scatter!(ax4, ntimes, τII_max_time, color = :black)
        # Colorbar(fg1[1,2][1,1], hm1, label = "P [Pa]")
        # Colorbar(fg1[1,2][1,2], hm2, label = "Vx [m.s⁻¹]")
    
        display(fg1)
    end

    # Return
    return nothing
end

# --------------------- #
#|   Helper functions  |#
# --------------------- #

# Diffuse array
@views function smooth_2DArray_diffusion!(A :: Matrix{Float64}, nsteps :: Int64, dx :: Float64, dy :: Float64)
    for _ in 1:nsteps
        A[2:end-1, 2:end-1] .+= (min.(dx.^2, dy.^2) ./ 1.0 ./ 4.1) .* ( diff( 1.0 .* diff(A[:, 2:end-1], dims = 1), dims = 1) ./ dx^2 + diff(1.0 .* diff(A[2:end-1, :], dims = 2), dims = 2) ./ dy^2)
    end

    # Return
    return nothing

end

# Viscoelastic rheology
@views function compute_viscoelastic_rheology!(
    ηve :: Matrix{Float64},
    ηs  :: Matrix{Float64},
    G   :: Matrix{Float64},
    dt  :: Float64
)

    # Viscoelastic rheology
    ηve .= 1.0 ./ (1.0 ./ (G .* dt) .+ 1.0 ./ (ηs))

    # Return
    return nothing

end

# Pseudo-transient parameters
@views function compute_pseudotransient_parameters!(
    ηve      :: Matrix{Float64},
    dτ_ρ     :: Matrix{Float64},
    G̃dτ      :: Matrix{Float64},
    Ṽ        :: Float64,
    max_LxLy :: Float64,
    Re       :: Float64,
    r̃        :: Float64
)
    # Compute PT params
    dτ_ρ .= Ṽ .* max_LxLy ./ Re ./ ηve
    G̃dτ  .= Ṽ^2 ./ dτ_ρ ./ (r̃ + 2.0)

    # Return
    return nothing

end

# Conservation of mass
@views function compute_pressure_residual!(
    vx    :: Matrix{Float64},
    vy    :: Matrix{Float64},
    Res_p :: Matrix{Float64},
    G̃dτ   :: Matrix{Float64},
    dx    :: Float64,
    dy    :: Float64,
    r̃     :: Float64
)
    # Pressure residual
    Res_p .= - r̃ .* G̃dτ .* (diff(vx, dims = 1) ./ dx .+ diff(vy, dims = 2) ./ dy)

    # Return
    return nothing

end
@views function update_pressure!(
    p     :: Matrix{Float64},
    Res_p :: Matrix{Float64}
)
    # Update pressure
    p .+= Res_p

    # Return
    return nothing
end

# Vorticity
@views function compute_vorticity!(
    ωxx :: Matrix{Float64},
    ωyy :: Matrix{Float64},
    ωxy :: Matrix{Float64},
    ωyx :: Matrix{Float64},
    vx  :: Matrix{Float64},
    vy  :: Matrix{Float64},
    dx  :: Float64,
    dy  :: Float64
)
    # Compute vorticity
    ωxx   .= 0.5 .* (diff(vx, dims = 1) ./ dx - diff(vx, dims = 1) ./dx)
    ωyy   .= 0.5 .* (diff(vy, dims = 2) ./ dy - diff(vy, dims = 2) ./dy)
    ωxy   .= 0.5 .* (diff(vy[:, 2:end-1], dims = 1) ./ dx - diff(vx[2:end-1, :], dims = 2) ./dy)
    ωyx   .= 0.5 .* (diff(vx[2:end-1, :], dims = 2) ./ dy - diff(vy[:, 2:end-1], dims = 1) ./dx)

    # Return
    return nothing

end

# Jaumann derivative
@views function compute_Jaumann_derivative!(
    Jxx     :: Matrix{Float64},
    Jyy     :: Matrix{Float64},
    Jxy     :: Matrix{Float64},
    ωxx     :: Matrix{Float64},
    ωyy     :: Matrix{Float64},
    ωxy     :: Matrix{Float64},
    ωyx     :: Matrix{Float64},
    τxx     :: Matrix{Float64},
    τyy     :: Matrix{Float64},
    τxy     :: Matrix{Float64},
    τxx_old :: Matrix{Float64},
    τyy_old :: Matrix{Float64},
    τxy_old :: Matrix{Float64},
    vx      :: Matrix{Float64},
    vy      :: Matrix{Float64},
    dt      :: Float64,
    dx      :: Float64,
    dy      :: Float64
)
    # Normal component 1
    Jxx              .= -τxx_old ./ dt                                                  # Previous time step
    Jxx[2:end,   :] .+= max.(0.0, vx[2:end-1,      :]) .* diff(τxx, dims = 1) ./ dx     # Advection
    Jxx[1:end-1, :] .+= min.(vx[2:end-1,      :], 0.0) .* diff(τxx, dims = 1) ./ dx
    Jxx[:,   2:end] .+= max.(0.0, vy[:,      2:end-1]) .* diff(τxx, dims = 2) ./ dy
    Jxx[:, 1:end-1] .+= min.(vy[:,      2:end-1], 0.0) .* diff(τxx, dims = 2) ./ dy
    Jxx             .-= 2.0 .* ωxx .* τxx                                               # Rotation
    Jxx[2:end-1, 2:end-1] .-= @av(ωxy .* τxy .+ ωxy .* τxy)

    # Normal component 2
    Jyy              .= -τyy_old ./ dt
    Jyy[2:end,   :] .+= max.(0.0, vx[2:end-1,      :]) .* diff(τyy, dims = 1) ./ dx
    Jyy[1:end-1, :] .+= min.(vx[2:end-1,      :], 0.0) .* diff(τyy, dims = 1) ./ dx
    Jyy[:,   2:end] .+= max.(0.0, vy[:,      2:end-1]) .* diff(τyy, dims = 2) ./ dy
    Jyy[:, 1:end-1] .+= min.(vy[:,      2:end-1], 0.0) .* diff(τyy, dims = 2) ./ dy
    Jyy             .-= 2.0 .* ωyy .* τyy
    Jyy[2:end-1, 2:end-1] .-= @av(ωyx .* τxy .+ ωyx .* τxy)

    # Shear component 1 (and 2 for viscous isotropic materials)
    Jxy              .= -τxy_old ./ dt
    Jxy[2:end,   :] .+= max.(0.0, @av(vx[2:end-1, :])) .* diff(τxy, dims = 1) ./ dx
    Jxy[1:end-1, :] .+= min.(@av(vx[2:end-1, :]), 0.0) .* diff(τxy, dims = 1) ./ dx
    Jxy[:,   2:end] .+= max.(0.0, @av(vy[:, 2:end-1])) .* diff(τxy, dims = 2) ./ dy
    Jxy[:, 1:end-1] .+= min.(@av(vy[:, 2:end-1]), 0.0) .* diff(τxy, dims = 2) ./ dy
    Jxy             .-= @av(ωxx) .* τxy .+ ωyx .* @av(τxx)
    Jxy             .-= ωxy .* @av(τyy) .+ @av(ωyy) .* τxy

    # Return
    return nothing

end

# Update stresses
@views function update_stresses!(
    τxx :: Matrix{Float64},
    τyy :: Matrix{Float64},
    τxy :: Matrix{Float64},
    Jxx :: Matrix{Float64},
    Jyy :: Matrix{Float64},
    Jxy :: Matrix{Float64},
    vx  :: Matrix{Float64},
    vy  :: Matrix{Float64},
    ηve :: Matrix{Float64},
    G   :: Matrix{Float64},
    G̃dτ :: Matrix{Float64},
    dx  :: Float64,
    dy  :: Float64
)
    # Update stresses
    τxx   .= ( τxx .+ 2.0 .* G̃dτ      .* ( -1.0 ./ (2.0 .* G)      .* Jxx .+ diff(vx, dims = 1) ./ dx) ) ./ (G̃dτ ./ ηve .+ 1.0)
    τyy   .= ( τyy .+ 2.0 .* G̃dτ      .* ( -1.0 ./ (2.0 .* G)      .* Jyy .+ diff(vy, dims = 2) ./ dy) ) ./ (G̃dτ ./ ηve .+ 1.0)
    τxy   .= ( τxy .+ 2.0 .* @av(G̃dτ) .* ( -1.0 ./ (2.0 .* @av(G)) .* Jxy .+ 0.5 .* (diff(vx[2:end-1,:], dims = 2) ./ dy .+ diff(vy[:, 2:end-1], dims = 1) ./ dx)) ) ./ (@av(G̃dτ) ./ @av(ηve) .+ 1.0)

    # Return
    return nothing

end

# Conservation of linear momentum
@views function compute_velocity_residual!(
    τxx    :: Matrix{Float64},
    τyy    :: Matrix{Float64},
    τxy    :: Matrix{Float64},
    p      :: Matrix{Float64},
    Res_vx :: Matrix{Float64},
    Res_vy :: Matrix{Float64},
    dτ_ρ   :: Matrix{Float64},
    dx     :: Float64,
    dy     :: Float64
)
    # Velocity residual
    Res_vx .= @av_xi(dτ_ρ) .* ( diff(τxx[:, 2:end-1], dims = 1) ./ dx .+ diff(τxy, dims = 2) ./ dy .- diff(p[:, 2:end-1], dims = 1) ./ dx )
    Res_vy .= @av_yi(dτ_ρ) .* ( diff(τyy[2:end-1, :], dims = 2) ./ dy .+ diff(τxy, dims = 1) ./ dx .- diff(p[2:end-1, :], dims = 2) ./ dy )

    # Return
    return nothing

end
@views function update_velocity!(
    vx     :: Matrix{Float64},
    vy     :: Matrix{Float64},
    Res_vx :: Matrix{Float64},
    Res_vy :: Matrix{Float64}
)
    # Update velocity
    vx[2:end-1, 2:end-1] .+= Res_vx
    vy[2:end-1, 2:end-1] .+= Res_vy

    # Return
    return nothing

end

# Boundary conditions
@views function set_BC_velocity!(
    vx     :: Matrix,
    vy     :: Matrix,
    Lx     :: Float64,
    Ly     :: Float64,
    ε̇_bg   :: Float64,
    BCtype :: String
)
    # No stress
    if BCtype == "No stress"
        vx[:,   1] .= vx[:,       2]
        vx[:, end] .= vx[:, end - 1]
        vy[1,   :] .= vy[2,       :]
        vy[end, :] .= vy[end - 1, :]
    end

    # Constant background strain rate
    if BCtype == "ε̇_bg const."
        vx[1,   :] .=  0.5 * Lx * ε̇_bg
        vx[end, :] .= -0.5 * Lx * ε̇_bg
        vy[:,   1] .= -0.5 * Ly * ε̇_bg
        vy[:, end] .=  0.5 * Ly * ε̇_bg
    end

    # Return
    return nothing

end

# --------------------- #
#|      Run main       |#
# --------------------- #
incompressible_viscoelastic_stokes2D()

