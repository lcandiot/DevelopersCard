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
macro d2_dx2(A) esc(:( ($A[idx + 1, idy    ] - 2.0 * $A[idx, idy] + $A[idx - 1, idy    ]) )) end
macro d2_dy2(A) esc(:( ($A[idx,     idy + 1] - 2.0 * $A[idx, idy] + $A[idx,     idy - 1]) )) end
macro av(A)    esc(:( ($A[idx, idy] + $A[idx + 1, idy] + $A[idx, idy + 1] + $A[idx + 1, idy + 1]) * 0.25 )) end
macro av_i(A)    esc(:( ($A[idx, idy] + $A[idx - 1, idy] + $A[idx, idy - 1] + $A[idx - 1, idy - 1]) * 0.25 )) end
macro av_xa(A) esc(:( ($A[idx, idy] + $A[idx + 1, idy]) * 0.5 )) end
macro av_ya(A) esc(:( ($A[idx, idy] + $A[idx, idy + 1]) * 0.5 )) end
macro av_xi(A) esc(:( ($A[idx, idy + 1] + $A[idx + 1, idy + 1]) * 0.5 )) end
macro av_yi(A) esc(:( ($A[idx + 1, idy] + $A[idx + 1, idy + 1]) * 0.5 )) end

@views function compressible_vevp_stokes2D()
    # Real values
    ρ_r      = 2800.0
    gy       = 9.81
    gx       = 0.0
    sec_year = 3600.0 * 24.0 * 365.25
    ηs_inc   = 1e20
    ηvp      = 1e18
    G_r      = 1e10
    K_r      = 2e11
    H_r      = -0e8      # Cohesion weakening factor
    C_r      = 1e70
    p0       = 101300.0
    pconf    = 0e6
    ρ0       = 2800.0
    φ_r      = 30.0 * π / 180.0
    ψ_r      =  0.0 * π / 180.0

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
    nx, ny   = 102, 102
    nt       = 15
    ϵ_tol    = 1e-8
    max_iter = 1e5
    ncheck   = 1000
    CFL      = 0.9 / sqrt(2.0) / 10.0
    dx, dy   = Lx / nx, Ly / ny
    _dx, _dy = 1.0 / dx, 1.0 / dy
    Ṽ        = min(dx,dy) * CFL
    Re       = 3.0 * sqrt(10.0) / 2.0 * π
    r̃        = 0.5
    rel      = 0.5

    # Mechanics switches
    is_buoyant = 0

    # Initialization
    AvNW     = zeros(Float64, nx + 1, ny + 1)
    AvNE     = zeros(Float64, nx + 1, ny + 1)
    AvSW     = zeros(Float64, nx + 1, ny + 1)
    AvSE     = zeros(Float64, nx + 1, ny + 1)
    wNW      =  ones(Float64, nx + 1, ny + 1)
    wNE      =  ones(Float64, nx + 1, ny + 1)
    wSW      =  ones(Float64, nx + 1, ny + 1)
    wSE      =  ones(Float64, nx + 1, ny + 1)
    vx       = zeros(Float64, nx + 1, ny + 2)
    vy       = zeros(Float64, nx + 2, ny + 1)
    τxx      = zeros(Float64, nx    , ny    )
    τyy      = zeros(Float64, nx    , ny    )
    τzz      = zeros(Float64, nx    , ny    )
    τxy      = zeros(Float64, nx    , ny    )
    τxyv     = zeros(Float64, nx + 1, ny + 1)
    τII      = zeros(Float64, nx    , ny    )
    ε̇xx      = zeros(Float64, nx    , ny    )
    ε̇yy      = zeros(Float64, nx    , ny    )
    ε̇zz      = zeros(Float64, nx    , ny    )
    ε̇xy      = zeros(Float64, nx    , ny    )
    ε̇xyv     = zeros(Float64, nx + 1, ny + 1)
    ε̇II      = zeros(Float64, nx    , ny    )
    Jxx      = zeros(Float64, nx    , ny    )
    Jyy      = zeros(Float64, nx    , ny    )
    Jzz      = zeros(Float64, nx    , ny    )
    Jxy      = zeros(Float64, nx    , ny    )
    Jxyv     = zeros(Float64, nx + 1, ny + 1)
    Jyx      = zeros(Float64, nx    , ny    )
    Jyxv     = zeros(Float64, nx + 1, ny + 1)
    ωxx      = zeros(Float64, nx    , ny    )
    ωyy      = zeros(Float64, nx    , ny    )
    ωxy      = zeros(Float64, nx    , ny    )
    ωxyv     = zeros(Float64, nx + 1, ny + 1)
    ωyx      = zeros(Float64, nx    , ny    )
    ωyxv     = zeros(Float64, nx + 1, ny + 1)
    dvxdτ    = zeros(Float64, nx - 1, ny    )
    dvydτ    = zeros(Float64, nx    , ny - 1)
    Res_vx   = zeros(Float64, nx - 1, ny    )
    Res_vy   = zeros(Float64, nx    , ny - 1)
    Res_p    = zeros(Float64, nx    , ny    )
    dpdτ     = zeros(Float64, nx    , ny    )
    G̃dτ      = zeros(Float64, nx    , ny    )
    G̃dτv     = zeros(Float64, nx + 1, ny + 1)
    dτ_ρ     = zeros(Float64, nx    , ny    )
    p        = zeros(Float64, nx    , ny    ) .+ pconf
    pcorr    = zeros(Float64, nx    , ny    )
    ∇v       = zeros(Float64, nx    , ny    )
    ηs       = zeros(Float64, nx    , ny    ) .+ η_bg
    ηve      = zeros(Float64, nx    , ny    )
    ηvev     = zeros(Float64, nx + 1, ny + 1)
    ηvep     = zeros(Float64, nx    , ny    )
    ρ        = zeros(Float64, nx    , ny    ) .+ ρ_r
    G        = zeros(Float64, nx    , ny    ) .+ G_r
    Gv       = zeros(Float64, nx + 1, ny + 1) .+ G_r
    K        = zeros(Float64, nx    , ny    ) .+ K_r
    H        = zeros(Float64, nx    , ny    ) .+ H_r
    C        = zeros(Float64, nx    , ny    ) .+ C_r
    C0       = zeros(Float64, nx    , ny    ) .+ C_r
    Pla      = zeros(Float64, nx    , ny    )
    F        = zeros(Float64, nx    , ny    )
    λ̇        = zeros(Float64, nx    , ny    )
    λ̇rel     = zeros(Float64, nx    , ny    )
    λ̇_old    = zeros(Float64, nx    , ny    )
    ∂Q∂τxx   = zeros(Float64, nx    , ny    )
    ∂Q∂τyy   = zeros(Float64, nx    , ny    )
    ∂Q∂τzz   = zeros(Float64, nx    , ny    )
    ∂Q∂τxy   = zeros(Float64, nx    , ny    )
    sinφ     = sin(φ_r) .* ones(Float64, nx, ny)
    sinψ     = sin(ψ_r) .* ones(Float64, nx, ny)
    cosφ     = cos(φ_r) .* ones(Float64, nx, ny)
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
    τxyv_old = deepcopy(τxyv)
    p_old    = deepcopy(p  )

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

    # Smooth out the inclusion
    smooth_2DArray_diffusion!(ηs, 10, dx, dy)

    # Set weights
    set_weigths!(wNW, wNE, wSW, wSE)

    # Visualize
    fg1   = Figure(size = (1600, 1600))
    ax1   = Axis(fg1[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Pt")
    ax2   = Axis(fg1[2, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vx")
    ax3   = Axis(fg1[2, 2], xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Vy")
    ax4   = Axis(fg1[1, 2], xlabel = "x", ylabel = "y", title = "ρ")
    ax5   = Axis(fg1[3, 1], xlabel = "nt", ylabel = "τII", title = "Stress build up")
    ax6   = Axis(fg1[3, 2], xlabel = "nt", ylabel = "τII", title = "∇v")
    hm1   = heatmap!(ax1, xc, yc, ηs)
    hm2   = heatmap!(ax2, xv, yv, vx)
    hm3   = heatmap!(ax3, xv, yv, vy)
    # display(fg1)

    # TIME LOOP
    time = 0.0; τII_max_time = []; ntimes = []
    for idx_Time in 1:nt

        # Update time
        # dt   = min(η_bg / G_r / 10.0, min(dx, dy) / maximum(max.(abs.(@av_arr(vx)), abs.(@av_arr(vy)))) / 4.1)
        dt   = η_bg / G_r / 10.0
        time = time + dt

        @printf("\ndt = %.5e\n", dt)

        # Update old physical fields
        τxx_old  .= τxx
        τyy_old  .= τyy
        τxyv_old .= τxyv
        τxy_old  .= @av_arr(τxyv)
        p_old    .= p

        # Pseudo - transient solver loop
        max_Res_all = []; err = 2ϵ_tol
        # compute_viscoelastic_rheology!(ηve, ηs, G, dt)
        for iter in 1:max_iter
            
            # λ̇_old .= λ̇rel
            compute_viscoelastic_rheology!(ηve, ηs, G, dt)
            compute_pseudotransient_parameters!(ηve, dτ_ρ, G̃dτ, Ṽ, max_LxLy, Re, r̃)
            c2v(G̃dτ, G̃dτv, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            c2v(ηve, ηvev, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            compute_divergence_vs!(∇v, vx, vy, _dx, _dy)
            
            compute_pressure_change!(∇v, dpdτ, p, p_old, K, G̃dτ, dt, r̃)
            update_pressure!(p, dpdτ)
            
            compute_vorticity!(ωxx, ωyy, ωxyv, ωyxv, vx, vy, _dx, _dy)
            # v2c(τxyv, τxy)
            # v2c(ωxyv, ωxy)
            # v2c(ωyxv, ωyx)
            compute_Jaumann_derivative!(Jxx, Jyy, Jxy, ωxx, ωyy, ωxy, ωyx, τxx, τyy, τxy, τxx_old, τyy_old, τxy_old, vx, vy, dt, _dx, _dy)
            c2v(Jxy, Jxyv, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            c2v(G̃dτ, G̃dτv, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            c2v(ηve, ηvev, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            # compute_2ndInv(τII, τxx, τyy, τzz, τxy)
            # @printf "min(τII ) = %.2e \t max(τII ) = %.2e\n" minimum(τII) maximum(τII)
            compute_strainrates(ε̇xx, ε̇yy, ε̇zz, ε̇xyv, vx, vy, ∇v, _dx, _dy)
            v2c(ε̇xyv, ε̇xy)
            # @printf "min(ε̇xyv) = %.2e \t max(ε̇xyv) = %.2e\n" minimum(ε̇xyv) maximum(ε̇xyv)
            # @printf "min(ε̇xy ) = %.2e \t max(ε̇xy ) = %.2e\n" minimum(ε̇xy) maximum(ε̇xy)
            update_stresses!(τxx, τyy, τzz, τxyv, ε̇xx, ε̇yy, ε̇zz, ε̇xyv, Jxx, Jyy, Jzz, Jxyv, ∇v, vx, vy, ηve, ηvev, G, Gv, G̃dτ, G̃dτv, _dx, _dy)
            # compute_shearstress_vertices!(τxyv, ε̇xyv, ηvev)
            # v2c(τxyv, τxy)
            # c2v(τxy, τxyv, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            compute_2ndInv(τII, τxx, τyy, τzz, τxyv)
            compute_velocity_change!(τxx, τyy, τxyv, p, dvxdτ, dvydτ, dτ_ρ, _dx, _dy)
            update_velocity!(vx, vy, dvxdτ, dvydτ)

            # @printf "min(τxx ) = %.2e \t max(τxx ) = %.2e\n" minimum(τxx) maximum(τxx)
            # @printf "min(τII ) = %.2e \t max(τII ) = %.2e\n" minimum(τII) maximum(τII)

            compute_2ndInv(ε̇II, ε̇xx, ε̇yy, ε̇zz, ε̇xyv)
            # check_yield_1(τII, p, C, F, sinφ, cosφ)
            # compute_plastic_flow_potential(Pla, F, H, K, ηve, sinφ, sinψ, cosφ, τII, τxx, τyy, τzz, τxy, λ̇, λ̇rel, λ̇_old, ∂Q∂τxx, ∂Q∂τyy,∂Q∂τzz, ∂Q∂τxy, ηvp, dt, rel)
            # correct_stress_pressure(K, τxx, τyy, τzz, τxy, ηve, ε̇xx, ε̇yy, ε̇zz, ε̇xy, λ̇rel, ∂Q∂τxx, ∂Q∂τyy, ∂Q∂τzz, ∂Q∂τxy, sinψ, pcorr, p, H, C, C0, dt)
            # @printf "min(λ̇   ) = %.2e \t max(λ̇   ) = %.2e\n" minimum(λ̇rel) maximum(λ̇rel)
            # @printf "min(Qxx ) = %.2e \t max(Qxx ) = %.2e\n" minimum(∂Q∂τxx) maximum(∂Q∂τxx)
            # @printf "min(τxx ) = %.2e \t max(τxx ) = %.2e\n" minimum(τxx) maximum(τxx)
            # @printf "min(τyy ) = %.2e \t max(τyy ) = %.2e\n" minimum(τyy) maximum(τyy)
            # @printf "min(τzz ) = %.2e \t max(τzz ) = %.2e\n" minimum(τzz) maximum(τzz)
            # @printf "min(τxy ) = %.2e \t max(τxy ) = %.2e\n" minimum(τxy) maximum(τxy)

            # compute_2ndInv(τII, τxx, τyy, τzz, τxy)
            # check_yield_2(τII, F, pcorr, sinφ, cosφ, C, λ̇rel, ηvp)
            # compute_viscoelastoplastic_rheology(ηvep, τII, ε̇II)
            # compute_pseudotransient_parameters!(ηve, dτ_ρ, G̃dτ, Ṽ, max_LxLy, Re, r̃)

            # Compute density
            compute_density!(ρ, p, K, ρ0, p0)

            # Boundary conditions
            set_BC_velocity!(vx, vy, Lx, Ly, ε̇_bg, "ε̇_bg const.")

            # Residuals
            compute_pressure_residual!(p, p_old, K, ∇v, Res_p, dt)
            compute_velocity_residual!(τxx, τyy, τxyv, p, Res_vx, Res_vy, _dx, _dy)

            # Monitor residuals 
            if iter % ncheck == 0
                @printf("\n")
                @printf("Iteration = %d \n", iter)
                @printf("Residual p  = %.2e \n", maximum(Res_p ) * dt)
                @printf("Residual vx = %.2e \n", maximum(Res_vx) / Psc * Lsc)
                @printf("Residual vy = %.2e \n", maximum(Res_vy) / Psc * Lsc)
                @printf("Physics:\n")
                @printf("min(vx) = %.2e; max(vx) = %.2e; min(vy) = %.2e; max(vy) = %.2e\n", minimum(vx), maximum(vx), minimum(vy), maximum(vy))
                @printf "min(ηvep) = %.2e \t max(ηvep) = %.2e\n" minimum(ηvp) maximum(ηvp)
                @printf "min(λ̇   ) = %.2e \t max(λ̇   ) = %.2e\n" minimum(λ̇) maximum(λ̇)
                @printf "min(C   ) = %.2e \t max(C   ) = %.2e\n" minimum(C) maximum(C)
                @printf "min(τxx ) = %.2e \t max(τxx ) = %.2e\n" minimum(τxx) maximum(τxx)
                @printf "min(τyy ) = %.2e \t max(τyy ) = %.2e\n" minimum(τyy) maximum(τyy)
                @printf "min(τzz ) = %.2e \t max(τzz ) = %.2e\n" minimum(τzz) maximum(τzz)
                @printf "min(τxy ) = %.2e \t max(τxy ) = %.2e\n" minimum(τxy) maximum(τxy)
                err = maximum([maximum(Res_p ) * dt maximum(Res_vx) / Psc * Lsc maximum(Res_vy) / Psc * Lsc])
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
        empty!(ax5)
        empty!(ax6)
        hm1 = heatmap!(ax1, xc, yc, p, colormap = :viridis)
        hm2 = heatmap!(ax2, xv, yv, vx, colormap = :roma)
        hm3 = heatmap!(ax3, xv, yv, vy, colormap = :roma)
        hm4 = heatmap!(ax4, xv, yv, ρ, colormap = :roma)
        ln  = scatter!(ax5, ntimes, τII_max_time, color = :black)
        hm4 = heatmap!(ax6, xc, yc, ∇v, colormap = :roma)
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

function set_weigths!(
    wNW :: Matrix{Float64},
    wNE :: Matrix{Float64},
    wSW :: Matrix{Float64},
    wSE :: Matrix{Float64},
)
    # Get sizes
    nx, ny = size(wNW)

    # Set weights
    for idy in 1:ny
        for idx in 1:nx

            # North West
            (idx == 1 && idy == 1                             ) ? wNW[idx, idy] = 4.0 : nothing # Left corner
            (idx >= 2 && idx <=nx-1 && idy == 1               ) ? wNW[idx, idy] = 2.0 : nothing # Left edge
            (idx == 1               && idy >= 2 && idy <= ny-1) ? wNW[idx, idy] = 2.0 : nothing # Top edge
            (idx == nx              && idy >= 1               ) ? wNW[idx, idy] = 0.0 : nothing # Bottom edge
            (idx >= 1               && idy == ny              ) ? wNW[idx, idy] = 0.0 : nothing # Right edge

            # North East
            (idx == 1 && idy == ny                            ) ? wNE[idx, idy] = 4.0 : nothing # Right corner
            (idx >= 2 && idx <=nx-1 && idy == ny              ) ? wNE[idx, idy] = 2.0 : nothing # Right edge
            (idx == 1               && idy >= 2 && idy <= ny-1) ? wNE[idx, idy] = 2.0 : nothing # Top edge
            (idx == nx              && idy >= 1               ) ? wNE[idx, idy] = 0.0 : nothing # Bottom edge
            (idx >= 1               && idy == 1               ) ? wNE[idx, idy] = 0.0 : nothing # Left edge

            # South West
            (idx == nx && idy == 1                             ) ? wSW[idx, idy] = 4.0 : nothing # Left corner
            (idx >= 2  && idx <=nx-1 && idy == 1               ) ? wSW[idx, idy] = 2.0 : nothing # Left edge
            (idx == nx               && idy >= 2 && idy <= ny-1) ? wSW[idx, idy] = 2.0 : nothing # Bottom edge
            (idx == 1                && idy >= 1               ) ? wSW[idx, idy] = 0.0 : nothing # Top edge
            (idx >= 1                && idy == ny              ) ? wSW[idx, idy] = 0.0 : nothing # Right edge

            # South East
            (idx == nx && idy == ny                           ) ? wSE[idx, idy] = 4.0 : nothing # Right corner
            (idx >= 2 && idx <=nx-1 && idy == ny              ) ? wSE[idx, idy] = 2.0 : nothing # Right edge
            (idx == nx              && idy >= 2 && idy <= ny-1) ? wSE[idx, idy] = 2.0 : nothing # Bottom edge
            (idx == 1               && idy >= 1               ) ? wSE[idx, idy] = 0.0 : nothing # Top edge
            (idx >= 1               && idy == 1               ) ? wSE[idx, idy] = 0.0 : nothing # left edge
        end
    end

    # Return
    return nothing
end

# Interpolate from centers to vertices
function c2v(
    Ac   :: Matrix{Float64},
    Av   :: Matrix{Float64},
    AvNW :: Matrix{Float64},
    AvNE :: Matrix{Float64},
    AvSW :: Matrix{Float64},
    AvSE :: Matrix{Float64},
    wNW  :: Matrix{Float64},
    wNE  :: Matrix{Float64},
    wSW  :: Matrix{Float64},
    wSE  :: Matrix{Float64}
)
    # Get sizes
    nvx, nvy = size(Av)

    # Interpolate
    for idy in 1:nvy
        for idx in 1:nvx
            (idx >= 1 && idx <= nvx-1 && idy >= 1 && idy <= nvy-1) ? AvNW[idx, idy] = Ac[idx,   idy  ] : nothing # NW
            (idx >= 1 && idx <= nvx-1 && idy >= 2 && idy <= nvy  ) ? AvNE[idx, idy] = Ac[idx,   idy-1] : nothing # NE
            (idx >= 2 && idx <= nvx   && idy >= 1 && idy <= nvy-1) ? AvSW[idx, idy] = Ac[idx-1, idy  ] : nothing # SW
            (idx >= 2 && idx <= nvx   && idy >= 2 && idy <= nvy  ) ? AvSE[idx, idy] = Ac[idx-1, idy-1] : nothing # SE
            Av[idx, idy] = 0.25 * (wNW[idx, idy] * AvNW[idx, idy] + wNE[idx, idy] * AvNE[idx, idy] + wSW[idx, idy] * AvSW[idx, idy] + wSE[idx, idy] * AvSE[idx, idy])
        end
    end

    # Return
    return nothing
end

# Interpolate from vertices to centers
function v2c(
    Av :: Matrix{Float64},
    Ac :: Matrix{Float64}
)
    # Get sizes
    ncx, ncy = size(Ac)
    
    # Interpolate
    for idy in 1:ncy
        for idx in 1:ncx
            Ac[idx, idy] = 0.25 * (Av[idx, idy] + Av[idx+1, idy] + Av[idx, idy+1] + Av[idx+1, idy+1])
        end
    end

    # Return
    return nothing
end

# Diffuse array
function smooth_2DArray_diffusion!(A :: Matrix{Float64}, nsteps :: Int64, dx :: Float64, dy :: Float64)
    nx, ny = size(A)
    for _ in 1:nsteps
        for idy in 1:ny
            for idx in 1:nx
                if idx > 1 && idx < nx && idy > 1 && idy < ny
                    A[idx, idy] += (min(dx^2, dy^2) / 1.0 / 4.1) * ( @d2_dx2(A) / dx^2 + @d2_dy2(A) / dy^2 )
                end
            end
        end
    end

    # Return
    return nothing

end

# Compute sqrt of 2nd tensor invariant
@views function compute_2ndInv(
    TII  :: Matrix{Float64},
    Txx  :: Matrix{Float64},
    Tyy  :: Matrix{Float64},
    Tzz  :: Matrix{Float64},
    Txyv :: Matrix{Float64},
)
    # Get size
    nx, ny = size(TII)

    # Compute sqrt of 2nd invariant
    for idy in 1:ny
        for idx in 1:nx
            TII[idx, idy] = sqrt( 0.5 * (Txx[idx, idy]^2 + Tyy[idx, idy]^2 + Tzz[idx, idy]^2 + 2.0 * @av(Txyv)^2))
        end
    end
end

# Divergence
function compute_divergence_vs!(
    ∇v :: Matrix{Float64},
    vx :: Matrix{Float64},
    vy :: Matrix{Float64},
    _dx :: Float64,
    _dy :: Float64
)
    # Get array sizes
    nx, ny = size(∇v)

    # Calculate divergence of solid velocity
    for idy in 1:ny
        for idx in 1:nx
            ∇v[idx, idy] = @d_dxi(vx) * _dx + @d_dyi(vy) * _dy
        end
    end

    # Return
    return nothing

end

# Density EOS
function compute_density!(
    ρ  :: Matrix{Float64},
    p  :: Matrix{Float64},
    K  :: Matrix{Float64},
    ρ0 :: Float64,
    p0 :: Float64
)
    # Get array sizes
    nx, ny = size(ρ)

    # Compute density with linearised EOS
    for idy in 1:ny
        for idx in 1:nx
            ρ[idx, idy] = ρ0 * (1.0 + 1.0 / K[idx, idy] * (p[idx, idy] - p0))
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

# Check yielding - part I
function check_yield_1(
    τII  :: Matrix{Float64},
    p    :: Matrix{Float64},
    C    :: Matrix{Float64},
    F    :: Matrix{Float64},
    sinφ :: Matrix{Float64},
    cosφ :: Matrix{Float64}
)
    # Get sizes
    nx, ny = size(p)

    # Evaluate yield
    for idy in 1:ny
        for idx in 1:nx
            F[idx, idy] = τII[idx, idy] - p[idx, idy] * sinφ[idx, idy] - C[idx, idy] * cosφ[idx, idy]
        end
    end
end

# Plastic multiplier and flow potential
function compute_plastic_flow_potential(
    Pla    :: Matrix{Float64},
    F      :: Matrix{Float64},
    H      :: Matrix{Float64},
    K      :: Matrix{Float64},
    ηve    :: Matrix{Float64},
    sinφ   :: Matrix{Float64},
    sinψ   :: Matrix{Float64},
    cosφ   :: Matrix{Float64},
    τII    :: Matrix{Float64},
    τxx    :: Matrix{Float64},
    τyy    :: Matrix{Float64},
    τzz    :: Matrix{Float64},
    τxy    :: Matrix{Float64},
    λ̇      :: Matrix{Float64},
    λ̇rel   :: Matrix{Float64},
    λ̇_old  :: Matrix{Float64},
    ∂Q∂τxx :: Matrix{Float64},
    ∂Q∂τyy :: Matrix{Float64},
    ∂Q∂τzz :: Matrix{Float64},
    ∂Q∂τxy :: Matrix{Float64},
    ηvp    :: Float64, 
    dt     :: Float64,
    rel    :: Float64
)
    # Get sizes
    nx, ny = size(Pla)

    # Flag plastic nodes
    Pla .= 0.0
    for idy in 1:ny
        for idx in 1:nx
            F[idx, idy] > 0.0 ? Pla[idx, idy] = 1.0 : nothing
        end
    end

    # Plastic multiplier
    for idy in 1:ny
        for idx in 1:nx
            λ̇[idx, idy] = Pla[idx, idy] * (F[idx, idy] / (ηve[idx, idy] + ηvp + K[idx, idy]*dt*sinφ[idx, idy]*sinψ[idx, idy] + H[idx, idy]*cosφ[idx, idy]*dt))
            λ̇rel[idx, idy] = rel * Pla[idx, idy] * (F[idx, idy] / (ηve[idx, idy] + ηvp + K[idx, idy]*dt*sinφ[idx, idy]*sinψ[idx, idy] + H[idx, idy]*cosφ[idx, idy]*dt)) + (1.0 - rel) * λ̇_old[idx, idy]
        end
    end

    # Derivative of flow potential
    for idy in 1:ny
        for idx in 1:nx
            ∂Q∂τxx[idx, idy] = 0.5 / τII[idx, idy] * τxx[idx, idy]
            ∂Q∂τyy[idx, idy] = 0.5 / τII[idx, idy] * τyy[idx, idy]
            ∂Q∂τzz[idx, idy] = 0.5 / τII[idx, idy] * τzz[idx, idy]
            if idx < nx && idy < ny
                ∂Q∂τxy[idx, idy] = 1.0 / @av(τII) * τxy[idx, idy]
            end
        end
    end

end

# Plastic corrections
function correct_stress_pressure(
    K   :: Matrix{Float64},
    τxx :: Matrix{Float64},
    τyy :: Matrix{Float64},
    τzz :: Matrix{Float64},
    τxy :: Matrix{Float64},
    ηve :: Matrix{Float64},
    ε̇xx :: Matrix{Float64},
    ε̇yy :: Matrix{Float64},
    ε̇zz :: Matrix{Float64},
    ε̇xy :: Matrix{Float64},
    λ̇   :: Matrix{Float64},
    ∂Q∂τxx :: Matrix{Float64},
    ∂Q∂τyy :: Matrix{Float64},
    ∂Q∂τzz :: Matrix{Float64},
    ∂Q∂τxy :: Matrix{Float64},
    sinψ   :: Matrix{Float64},
    pcorr  :: Matrix{Float64},
    p      :: Matrix{Float64},
    H      :: Matrix{Float64},
    C      :: Matrix{Float64},
    C0     :: Matrix{Float64},
    dt     :: Float64
)
    # Get sizes
    nx, ny = size(τxx)

    # Corrections
    for idy in 1:ny
        for idx in 1:nx
            τxx[idx, idy] = 2.0 * ηve[idx, idy] * (ε̇xx[idx, idy] - λ̇[idx, idy] * ∂Q∂τxx[idx, idy])
            τyy[idx, idy] = 2.0 * ηve[idx, idy] * (ε̇yy[idx, idy] - λ̇[idx, idy] * ∂Q∂τyy[idx, idy])
            τzz[idx, idy] = 2.0 * ηve[idx, idy] * (ε̇zz[idx, idy] - λ̇[idx, idy] * ∂Q∂τzz[idx, idy])
            pcorr[idx, idy] = p[idx, idy] + K[idx, idy] * dt * λ̇[idx, idy] * sinψ[idx, idy]
            C[idx, idy] = C0[idx, idy] + dt * sqrt(2/3) * λ̇[idx, idy] * H[idx, idy]
            if idx < nx && idy < ny
                τxy[idx, idy] = 2.0 * @av(ηve) * (ε̇xy[idx, idy] - 0.5 * @av(λ̇) * ∂Q∂τxy[idx, idy])
            end
        end
    end
end

# Check yield - part II
function check_yield_2(
    τII :: Matrix{Float64},
    F   :: Matrix{Float64},
    pcorr :: Matrix{Float64},
    sinφ  :: Matrix{Float64},
    cosφ  :: Matrix{Float64},
    C     :: Matrix{Float64},
    λ̇     :: Matrix{Float64},
    ηvp   :: Float64
)
    # Get sizes
    nx, ny = size(τII)

    # Check yield
    for idy in 1:ny
        for idx in 1:nx
            F[idx, idy] = τII[idx, idy] - pcorr[idx, idy] * sinφ[idx, idy] - C[idx, idy] * cosφ[idx, idy] - ηvp * λ̇[idx, idy]
        end
    end
end

# Visco-elasto-plastic rheology
function compute_viscoelastoplastic_rheology(
    ηvep :: Matrix{Float64},
    τII  :: Matrix{Float64},
    ε̇II  :: Matrix{Float64},
)
    ηvep .= τII ./ 2.0 ./ ε̇II
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
function compute_pressure_change!(
    ∇v    :: Matrix{Float64},
    dpdτ :: Matrix{Float64},
    p     :: Matrix{Float64},
    p_old :: Matrix{Float64},
    K     :: Matrix{Float64},
    G̃dτ   :: Matrix{Float64},
    dt    :: Float64,
    r̃     :: Float64
)
    # Get sizes
    nx, ny = size(dpdτ)

    # Pressure residual
    for idy in 1:ny
        for idx in 1:nx
            dpdτ[idx, idy] = - r̃ * G̃dτ[idx, idy] * (∇v[idx, idy] + 1.0 / K[idx, idy] * (p[idx, idy] - p_old[idx, idy]) / dt)
        end
    end

    # Return
    return nothing

end

# Pressure residual
function compute_pressure_residual!(
    p     :: Matrix{Float64},
    p_old :: Matrix{Float64},
    K     :: Matrix{Float64},
    ∇v    :: Matrix{Float64},
    Res_p :: Matrix{Float64},
    dt    :: Float64
)
    # Get size
    nx, ny = size(p)

    # Compute residual
    for idy in 1:ny
        for idx in 1:nx
            Res_p[idx, idy] = -∇v[idx, idy] - 1.0 / (K[idx, idy] * dt) * (p[idx, idy] - p_old[idx, idy])
        end
    end
end

# Pressure
function update_pressure!(
    p    :: Matrix{Float64},
    dpdτ :: Matrix{Float64}
)
    # Get sizes
    nx, ny = size(p)

    # Update pressure
    for idy in 1:ny
        for idx in 1:nx
            p[idx, idy] += dpdτ[idx, idy]
        end
    end

    # Return
    return nothing
end

# Vorticity
function compute_vorticity!(
    ωxx  :: Matrix{Float64},
    ωyy  :: Matrix{Float64},
    ωxyv :: Matrix{Float64},
    ωyxv :: Matrix{Float64},
    vx   :: Matrix{Float64},
    vy   :: Matrix{Float64},
    _dx  :: Float64,
    _dy  :: Float64
)
    # Get sizes
    nx, ny = size(ωxx)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx
            ωxx[idx, idy] = 0.5 * (@d_dxi(vx) * _dx - @d_dxi(vx) * _dx)
            ωyy[idx, idy] = 0.5 * (@d_dyi(vy) * _dy - @d_dyi(vy) * _dy)
        end
    end

    # Shear components
    for idy in 1:ny + 1
        for idx in 1:nx + 1
            ωxyv[idx, idy] = 0.5 * (@d_dx(vy) * _dx - @d_dy(vx) * _dy)
            ωyxv[idx, idy] = 0.5 * (@d_dy(vx) * _dy - @d_dx(vy) * _dx)
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
    _dx      :: Float64,
    _dy      :: Float64
)
    # Get sizes
    nx, ny = size(Jxx)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx

            # Old values
            Jxx[idx, idy] = -τxx_old[idx, idy] / dt
            Jyy[idx, idy] = -τyy_old[idx, idy] / dt
            Jxy[idx, idy] = -τxy_old[idx, idy] / dt

            # Advection
            if idx > 1 && idx < nx && idy > 1 && idy < ny
                Jxx[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τxx) * _dx
                Jyy[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τyy) * _dx
                Jxx[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τxx) * _dy
                Jyy[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τyy) * _dy
                Jxy[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τxy) * _dx
                Jxy[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τxy) * _dy
            end
            if idx < nx && idy < ny
                Jxx[idx, idy] += min(vx[idx, idy], 0.0) * @d_dx(τxx) * _dx
                Jyy[idx, idy] += min(vx[idx, idy], 0.0) * @d_dx(τyy) * _dx
                Jxx[idx, idy] += min(vy[idx, idy], 0.0) * @d_dy(τxx) * _dy
                Jyy[idx, idy] += min(vy[idx, idy], 0.0) * @d_dy(τyy) * _dy
                Jxy[idx, idy] += min(vx[idx, idy], 0.0) * @d_dx(τxy) * _dx
                Jxy[idx, idy] += min(vy[idx, idy], 0.0) * @d_dy(τxy) * _dy
            end

            # Rotation
            Jxx[idx, idy] -= 2.0 * ωxx[idx, idy] * τxx[idx, idy]
            Jyy[idx, idy] -= 2.0 * ωyy[idx, idy] * τyy[idx, idy]
            Jxx[idx, idy] -= ωxy[idx, idy] * τxy[idx, idy] + ωxy[idx, idy] * τxy[idx, idy]
            Jyy[idx, idy] -= ωyx[idx, idy] * τxy[idx, idy] + ωyx[idx, idy] * τxy[idx, idy]
            Jxy[idx, idy] -= ωxx[idx, idy] * τxy[idx, idy] + ωyx[idx, idy] * τxx[idx, idy]
            Jxy[idx, idy] -= ωxy[idx, idy] * τyy[idx, idy] + ωyy[idx, idy] * τxy[idx, idy]
        end
    end

    # # Shear component
    # for idy in 1:ny
    #     for idx in 1:nx

    #         # Old values
    #         Jxy[idx, idy] = -τxy_old[idx, idy] / dt

    #         # Advection
    #         if idx > 1 && idx < nx - 1 && idy > 1 && idy < ny - 1
    #             Jxy[idx, idy] += max(0.0, @av(vx)) * @d_dx(τxy) * _dx
    #             Jxy[idx, idy] += max(0.0, @av(vy)) * @d_dy(τxy) * _dy
    #         end
    #         if idx < nx - 1 && idy < ny - 1
    #             Jxy[idx, idy] += min(@av(vx), 0.0) * @d_dx(τxy) * _dx
    #             Jxy[idx, idy] += min(@av(vy), 0.0) * @d_dy(τxy) * _dy
    #         end

    #         # Rotation
    #         Jxy[idx, idy] -= @av(ωxx) * τxy[idx, idy] + ωyx[idx, idy] .* @av(τxx)
    #         Jxy[idx, idy] -= ωxy[idx, idy] * @av(τyy) + @av(ωyy) * τxy[idx, idy]
    #     end
    # end

    # Return
    return nothing

end

# Strain rates
function compute_strainrates(
    ε̇xx  :: Matrix{Float64},
    ε̇yy  :: Matrix{Float64},
    ε̇zz  :: Matrix{Float64},
    ε̇xyv :: Matrix{Float64},
    vx   :: Matrix{Float64},
    vy   :: Matrix{Float64},
    ∇v   :: Matrix{Float64},
    _dx  :: Float64,
    _dy  :: Float64
)
    # Get sizes
    nx, ny = size(ε̇xx)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx
            ε̇xx[idx, idy] = @d_dxi(vx) * _dx - ∇v[idx, idy] / 3.0
            ε̇yy[idx, idy] = @d_dyi(vy) * _dy - ∇v[idx, idy] / 3.0
            ε̇zz[idx, idy] =                  - ∇v[idx, idy] / 3.0
        end
    end

    # Shear component
    for idy in 1:ny+1
        for idx in 1:nx+1
            ε̇xyv[idx, idy] = 0.5 * (@d_dx(vy) * _dx + @d_dy(vx) * _dy)
        end
    end
end

# Update stresses
function update_stresses!(
    τxx :: Matrix{Float64},
    τyy :: Matrix{Float64},
    τzz :: Matrix{Float64},
    τxyv :: Matrix{Float64},
    ε̇xx :: Matrix{Float64},
    ε̇yy :: Matrix{Float64},
    ε̇zz :: Matrix{Float64},
    ε̇xyv :: Matrix{Float64},
    Jxx :: Matrix{Float64},
    Jyy :: Matrix{Float64},
    Jzz :: Matrix{Float64},
    Jxyv :: Matrix{Float64},
    ∇v  :: Matrix{Float64},
    vx  :: Matrix{Float64},
    vy  :: Matrix{Float64},
    ηve :: Matrix{Float64},
    ηvev:: Matrix{Float64},
    G   :: Matrix{Float64},
    Gv  :: Matrix{Float64},
    G̃dτ :: Matrix{Float64},
    G̃dτv:: Matrix{Float64},
    _dx :: Float64,
    _dy :: Float64
)
    # Get sizes
    nx, ny = size(G)

    # Normal components
    for idy in 1:ny
        for idx in 1:nx
            τxx[idx, idy] = ( τxx[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jxx[idx, idy] + ε̇xx[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            τyy[idx, idy] = ( τyy[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jyy[idx, idy] + ε̇yy[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            # τzz[idx, idy] = ( τyy[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jzz[idx, idy] - 1.0/3.0 * ∇v[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
        end
    end

    # Shear component
    for idy in 1:ny+1
        for idx in 1:nx+1
            τxyv[idx, idy] = ( τxyv[idx, idy] + 2.0 * G̃dτv[idx, idy] * ( -1.0 / (2.0 * Gv[idx, idy]) * Jxyv[idx, idy] + ε̇xyv[idx, idy]) ) / (G̃dτv[idx, idy] / ηvev[idx, idy] + 1.0)
        end
    end

    # Return
    return nothing
end

# Shear stress on vertices
function compute_shearstress_vertices!(
    τxyv :: Matrix{Float64},
    ε̇xyv :: Matrix{Float64},
    ηvev :: Matrix{Float64}
)
    # Get size
    nxv, nyv = size(τxyv)

    # Compute
    for idy in 1:nyv
        for idx in 1:nxv
            τxyv[idx, idy] = 2.0 * ηvev[idx, idy] * ε̇xyv[idx, idy]
        end
    end

    # Return
    return
end

# Conservation of linear momentum
function compute_velocity_change!(
    τxx    :: Matrix{Float64},
    τyy    :: Matrix{Float64},
    τxyv   :: Matrix{Float64},
    p      :: Matrix{Float64},
    dvxdτ  :: Matrix{Float64},
    dvydτ  :: Matrix{Float64},
    dτ_ρ   :: Matrix{Float64},
    _dx    :: Float64,
    _dy    :: Float64
)
    # Get sizes
    nx , ny = size(p)

    # Velocity residual
    for idy in 1:ny
        for idx in 1:nx
            if idx < nx
                dvxdτ[idx, idy] = @av_xa(dτ_ρ) * ( @d_dx(τxx) * _dx + @d_dyi(τxyv) * _dy - @d_dx(p) * _dx )
            end
        end
    end
    for idy in 1:ny
        for idx in 1:nx
            if idy < ny
                dvydτ[idx, idy] = @av_ya(dτ_ρ) * ( @d_dy(τyy) * _dy + @d_dxi(τxyv) * _dx - @d_dy(p) * _dy )
            end
        end
    end
    # Return
    return nothing

end

# Velocity residual
function compute_velocity_residual!(
    τxx    :: Matrix{Float64},
    τyy    :: Matrix{Float64},
    τxyv   :: Matrix{Float64},
    p      :: Matrix{Float64},
    Res_vx :: Matrix{Float64},
    Res_vy :: Matrix{Float64},
    _dx     :: Float64,
    _dy     :: Float64
)
    # Get size
    nx, ny = size(p)

    # Compute residual
    for idy in 1:ny
        for idx in 1:nx
            if idx < nx
                Res_vx[idx, idy] = @d_dx(τxx) * _dx + @d_dyi(τxyv) * _dy - @d_dx(p) * _dx
            end
        end
    end
    for idy in 1:ny
        for idx in 1:nx
            if idy < ny
                Res_vy[idx, idy] = @d_dy(τyy) * _dy + @d_dxi(τxyv) * _dx - @d_dy(p) * _dy
            end
        end
    end
    # Return
    return nothing
end

function update_velocity!(
    vx     :: Matrix{Float64},
    vy     :: Matrix{Float64},
    dvxdτ  :: Matrix{Float64},
    dvydτ  :: Matrix{Float64}
)
    # Get sizes
    nx_vx, ny_vx = size(vx)[1], size(vx)[2]
    nx_vy, ny_vy = size(vy)[1], size(vy)[2]

    # @printf "nx_vx = %d \t ny_vx = %d \t nx_vy = %d \t ny_vy = %d\n" nx_vx nx_vy nx_vy ny_vy
    # error("Hello...")
    # Update velocity
    for idy in 1:ny_vx
        for idx in 1:nx_vx
            if idx < nx_vx-1 && idy < ny_vx-1
                vx[idx + 1, idy + 1] += dvxdτ[idx, idy]
            end
         end
    end
    for idy in 1:ny_vy
        for idx in 1:nx_vy
            if idy < ny_vy-1 && idx < nx_vy-1
                vy[idx + 1, idy + 1] += dvydτ[idx, idy]
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
    nx, ny = size(vx)[1]-1, size(vy)[2]-1

    # Velocity in x
    # vx[1   , 2:end-1] .=  0.5 * Lx * ε̇_bg
    # vx[nx+1, 2:end-1] .= -0.5 * Lx * ε̇_bg

    # # Velocity in y
    # vy[2:end-1,    1] .= -0.5 * Ly * ε̇_bg
    # vy[2:end-1, ny+1] .=  0.5 * Ly * ε̇_bg

    # vx[:, [1, end]] .= vx[:, [2, end-1]]
    # vy[[1, end], :] .= vy[[2, end-1], :]

    # Return
    return nothing
end

# --------------------- #
#|      Run main       |#
# --------------------- #
compressible_vevp_stokes2D();