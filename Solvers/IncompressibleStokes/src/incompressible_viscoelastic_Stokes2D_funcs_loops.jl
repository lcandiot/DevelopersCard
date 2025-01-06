# Solving 2D stokes flow of an incompressible viscoelastic fluid putting the building blocks into functions for speeding things up
using CairoMakie, Printf, ColorSchemes

macro av_arr(A)    esc(:( ($A[1:end-1, 1:end-1] .+ $A[2:end, 1:end-1] .+ $A[1:end-1, 2:end] .+ $A[2:end, 2:end]) .* 0.25 )) end
macro av_xi_arr(A) esc(:( ($A[1:end-1, 2:end-1] .+ $A[2:end, 2:end-1]) .* 0.5 )) end
macro av_yi_arr(A) esc(:( ($A[2:end-1, 1:end-1] .+ $A[2:end-1, 2:end]) .* 0.5 )) end
macro av_xa_arr(A) esc(:( ($A[1:end-1, :] .+ $A[2:end, :]) .* 0.5 )) end
macro av_ya_arr(A) esc(:( ($A[:, 1:end-1] .+ $A[:, 2:end]) .* 0.5 )) end

macro d_dx(A)  esc(:( ($A[idx + 1, idy    ] - $A[idx, idy     ]))) end
macro d_dy(A)  esc(:( ($A[idx    , idy + 1] - $A[idx, idy     ]))) end
macro d_dxi(A)  esc(:( ($A[idx + 1, idy + 1] - $A[idx, idy + 1]))) end
macro d_dyi(A)  esc(:( ($A[idx + 1, idy + 1] - $A[idx + 1, idy]))) end
macro d2_dx(A) esc(:( ($A[idx + 1, idy    ] - 2.0 * $A[idx, idy] + $A[idx - 1, idy    ]) )) end
macro d2_dy(A) esc(:( ($A[idx,     idy + 1] - 2.0 * $A[idx, idy] + $A[idx,     idy - 1]) )) end
macro av(A)    esc(:( ($A[idx, idy] + $A[idx + 1, idy] + $A[idx, idy + 1] + $A[idx + 1, idy + 1]) * 0.25 )) end
macro av_i(A)    esc(:( ($A[idx, idy] + $A[idx - 1, idy] + $A[idx, idy - 1] + $A[idx - 1, idy - 1]) * 0.25 )) end
macro av_xa(A) esc(:( ($A[idx, idy] + $A[idx + 1, idy]) * 0.5 )) end
macro av_ya(A) esc(:( ($A[idx, idy] + $A[idx, idy + 1]) * 0.5 )) end
macro av_xi(A) esc(:( ($A[idx, idy + 1] + $A[idx + 1, idy + 1]) * 0.5 )) end
macro av_yi(A) esc(:( ($A[idx + 1, idy] + $A[idx + 1, idy + 1]) * 0.5 )) end

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
    # hm1   = heatmap!(ax1, xc, yc, ηs)
    # hm2   = heatmap!(ax2, xv, yv, vx)
    # hm3   = heatmap!(ax3, xv, yv, vy)
    # display(fg1)

    # TIME LOOP
    time = 0.0; τII_max_time = []; ntimes = []
    for idx_Time in 1:nt

        # Update time
        # dt   = min(η_bg / G_r / 10.0, min(dx, dy) / maximum(max.(abs.(@av_xa(vx)), abs.(@av_ya(vy)))) / 4.1)
        dt   = η_bg / G_r / 10.0
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
            τII .= sqrt.( 0.5 .* (@av_arr(τxx).^2 .+ @av_arr(τyy).^2 .+ 2.0 .* τxy.^2) )

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
    
        # Store maximum shear stress to monitor build up
        push!(τII_max_time, maximum(τII))
        push!(ntimes, idx_Time)

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
#|  Compute functions  |#
# --------------------- #

# Diffuse array
function smooth_2DArray_diffusion!(A :: Matrix{Float64}, nsteps :: Int64, dx :: Float64, dy :: Float64)
    nx, ny = size(A)
    for _ in 1:nsteps
        for idy in 1:ny
            for idx in 1:nx
                if idx > 1 && idx < nx && idy > 1 && idy < ny
                    A[idx, idy] += (min(dx^2, dy^2) / 1.0 / 4.1) * ( @d2_dx(A) / dx^2 + @d2_dy(A) / dy^2 )
                end
            end
        end
    end

    # Return
    return nothing

end

# Viscoelastic rheology
function compute_viscoelastic_rheology!(
    ηve :: Matrix{Float64},
    ηs  :: Matrix{Float64},
    G   :: Matrix{Float64},
    dt  :: Float64
)
    # Get array sizes
    nx, ny = size(ηve)

    # Viscoelastic rheology
    for idy in 1:ny
        for idx in 1:nx
            ηve[idx, idy] = 1.0 / (1.0 / (G[idx, idy] * dt) + 1.0 / (ηs[idx, idy]))
        end
    end

    # Return
    return nothing

end

# Pseudo-transient parameters
function compute_pseudotransient_parameters!(
    ηve      :: Matrix{Float64},
    dτ_ρ     :: Matrix{Float64},
    G̃dτ      :: Matrix{Float64},
    Ṽ        :: Float64,
    max_LxLy :: Float64,
    Re       :: Float64,
    r̃        :: Float64
)
    # Get array sizes
    nx, ny = size(ηve)

    # Compute PT params
    for idy in 1:ny
        for idx in 1:nx
            dτ_ρ[idx, idy] = Ṽ * max_LxLy / Re / ηve[idx, idy]
            G̃dτ[idx, idy]  = Ṽ^2 / dτ_ρ[idx, idy] / (r̃ + 2.0)
        end
    end

    # Return
    return nothing

end

# Conservation of mass
function compute_pressure_residual!(
    vx    :: Matrix{Float64},
    vy    :: Matrix{Float64},
    Res_p :: Matrix{Float64},
    G̃dτ   :: Matrix{Float64},
    dx    :: Float64,
    dy    :: Float64,
    r̃     :: Float64
)
    # Get sizes
    nx, ny = size(Res_p)

    # Pressure residual
    for idy in 1:ny
        for idx in 1:nx
            Res_p[idx, idy] = - r̃ * G̃dτ[idx, idy] * (@d_dx(vx) / dx + @d_dy(vy) / dy)
        end
    end

    # Return
    return nothing

end
function update_pressure!(
    p     :: Matrix{Float64},
    Res_p :: Matrix{Float64}
)
    # Get sizes
    nx, ny = size(p)

    # Update pressure
    for idy in 1:ny
        for idx in 1:nx
            p[idx, idy] += Res_p[idx, idy]
        end
    end

    # Return
    return nothing
end

# Vorticity
function compute_vorticity!(
    ωxx :: Matrix{Float64},
    ωyy :: Matrix{Float64},
    ωxy :: Matrix{Float64},
    ωyx :: Matrix{Float64},
    vx  :: Matrix{Float64},
    vy  :: Matrix{Float64},
    dx  :: Float64,
    dy  :: Float64
)
    # Get sizes
    nx, ny = size(ωxx)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx
            ωxx[idx, idy] = 0.5 * (@d_dx(vx) / dx - @d_dx(vx) / dx)
            ωyy[idx, idy] = 0.5 * (@d_dy(vy) / dy - @d_dy(vy) / dy)
        end
    end

    # Shear components
    for idy in 1:ny - 1
        for idx in 1:nx - 1
            ωxy[idx, idy] = 0.0
            ωyx[idx, idy] = 0.0
            if idy > 1 && idy < ny
                ωxy[idx, idy] += 0.5 * (@d_dx(vy) / dx)
                ωyx[idx, idy] -= 0.5 * (@d_dx(vy) / dx)
            end
            if idx > 1 && idx < nx
                ωxy[idx, idy] -= 0.5 * (@d_dy(vx) / dy)
                ωyx[idx, idy] += 0.5 * (@d_dy(vx) / dy)
            end
        end
    end

    # Return
    return nothing

end

# Jaumann derivative
function compute_Jaumann_derivative!(
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
    # Get sizes
    nx, ny = size(Jxx)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx

            # Old values
            Jxx[idx, idy] = -τxx_old[idx, idy] / dt
            Jyy[idx, idy] = -τyy_old[idx, idy] / dt

            # Advection
            if idx > 1 && idx < nx && idy > 1 && idy < ny
                Jxx[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τxx) / dx
                Jyy[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τyy) / dx
                Jxx[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τxx) / dy
                Jyy[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τyy) / dy
            end
            if idx < nx && idy < ny
                Jxx[idx, idy] += min(vx[idx, idy], 0.0) * @d_dx(τxx) / dx
                Jyy[idx, idy] += min(vx[idx, idy], 0.0) * @d_dx(τyy) / dx
                Jxx[idx, idy] += min(vy[idx, idy], 0.0) * @d_dy(τxx) / dy
                Jyy[idx, idy] += min(vy[idx, idy], 0.0) * @d_dy(τyy) / dy
            end

            # Rotation
            Jxx[idx, idy] -= 2.0 * ωxx[idx, idy] * τxx[idx, idy]
            Jyy[idx, idy] -= 2.0 * ωyy[idx, idy] * τyy[idx, idy]
            if idx > 1 && idx < nx && idy > 1 && idy < ny
                Jxx[idx, idy] -= @av_i(ωxy) * @av_i(τxy) + @av_i(ωxy) * @av_i(τxy)
                Jyy[idx, idy] -= @av_i(ωyx) * @av_i(τxy) + @av_i(ωyx) * @av_i(τxy)
            end
        end
    end

    # Shear component
    for idy in 1:ny - 1
        for idx in 1:nx - 1

            # Old values
            Jxy[idx, idy] = -τxy_old[idx, idy] / dt

            # Advection
            if idx > 1 && idx < nx - 1 && idy > 1 && idy < ny - 1
                Jxy[idx, idy] += max(0.0, @av(vx)) * @d_dx(τxy) / dx
                Jxy[idx, idy] += max(0.0, @av(vy)) * @d_dy(τxy) / dy
            end
            if idx < nx - 1 && idy < ny - 1
                Jxy[idx, idy] += min(@av(vx), 0.0) * @d_dx(τxy) / dx
                Jxy[idx, idy] += min(@av(vy), 0.0) * @d_dy(τxy) / dy
            end

            # Rotation
            Jxy[idx, idy] -= @av(ωxx) * τxy[idx, idy] + ωyx[idx, idy] .* @av(τxx)
            Jxy[idx, idy] -= ωxy[idx, idy] * @av(τyy) + @av(ωyy) * τxy[idx, idy]
        end
    end

    # Return
    return nothing

end

# Update stresses
function update_stresses!(
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
    # Get sizes
    nx, ny = size(G)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx
            τxx[idx, idy] = ( τxx[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jxx[idx, idy] + @d_dx(vx) / dx) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            τyy[idx, idy] = ( τyy[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jyy[idx, idy] + @d_dy(vy) / dy) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
        end
    end

    # Shear component
    for idy in 1:ny - 1
        for idx in 1:nx - 1
            τxy[idx, idy] = ( τxy[idx, idy] + 2.0 * @av(G̃dτ) * ( -1.0 / (2.0 * @av(G)) * Jxy[idx, idy] + 0.5 * (@d_dyi(vx) / dy + @d_dxi(vy) / dx)) ) / (@av(G̃dτ) / @av(ηve) + 1.0)
        end
    end

    # Return
    return nothing

end

# Conservation of linear momentum
function compute_velocity_residual!(
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
    # Get sizes
    nx , ny = size(p)

    # Velocity residual
    for idy in 1:ny - 1
        for idx in 1:nx - 1
            if idy < ny - 1
                Res_vx[idx, idy] = @av_xi(dτ_ρ) * ( @d_dxi(τxx) / dx + @d_dy(τxy) / dy - @d_dxi(p) / dx )
            end
            if idx < nx - 1
                Res_vy[idx, idy] = @av_yi(dτ_ρ) * ( @d_dyi(τyy) / dy + @d_dx(τxy) / dx - @d_dyi(p) / dy )
            end
        end
    end

    # Return
    return nothing

end
function update_velocity!(
    vx     :: Matrix{Float64},
    vy     :: Matrix{Float64},
    Res_vx :: Matrix{Float64},
    Res_vy :: Matrix{Float64}
)
    # Get sizes
    nx, ny = size(vx)[1] - 1, size(vx)[2]

    # Update velocity
    for idy in 1:ny-1
        for idx in 1:nx-1
            if idy < ny - 1
                vx[idx + 1, idy + 1] += Res_vx[idx, idy]
            end
            if idx < nx - 1
                vy[idx + 1, idy + 1] += Res_vy[idx, idy]
            end
        end
    end

    # Return
    return nothing

end

# Boundary conditions
function set_BC_velocity!(
    vx     :: Matrix,
    vy     :: Matrix,
    Lx     :: Float64,
    Ly     :: Float64,
    ε̇_bg   :: Float64,
    BCtype :: String
)
    # Get sizes
    nx, ny = size(vx)[1] - 1, size(vx)[2]

    # Velocity in x
    for idy in 1:ny
        for _ in 1:nx + 1
             # Constant background strain rate
            if BCtype == "ε̇_bg const."
                vx[1,      idy] =  0.5 * Lx * ε̇_bg
                vx[nx + 1, idy] = -0.5 * Lx * ε̇_bg
            end
        end
    end

    # Velocity in y
    for _ in 1:ny + 1
        for idx in 1:nx
            # Constant background strain rate
            if BCtype == "ε̇_bg const."
                vy[idx,      1] = -0.5 * Ly * ε̇_bg
                vy[idx, ny + 1] =  0.5 * Ly * ε̇_bg
            end
        end
    end

    # Return
    return nothing

end

# --------------------- #
#|      Run main       |#
# --------------------- #
incompressible_viscoelastic_stokes2D()