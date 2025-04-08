# Solving 2D stokes flow of an incompressible viscoelastic fluid putting the building blocks into functions for speeding things up

# Options
const save_figure = false

# Load packages
using CairoMakie, Printf, ColorSchemes, Statistics

# Define macros
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
    ηs_inc   = 1e50
    ηvp      = 1e18
    G_r      = 1e10
    K_r      = 2e10
    H_r      = -2e8      # Cohesion weakening factor
    C_r      = 1.75e7
    p0       = 101300.0
    pconf    = 250e6
    ρ0       = 2800.0
    φ_r      = 30.0 * π / 180.0
    ψ_r      =  0.0 * π / 180.0

    # Indepentent physics
    ρg       = ρ_r * gy
    ε̇_bg     = 2e-13
    η_bg     = 1e50

    # Nondimensional number
    Lx_Ly    = 1.35
    Lx_Lc    = 1.0 #260.0

    # Dependent scales
    tsc      = 8e8
    Psc      = min(η_bg, G_r*tsc) * ε̇_bg
    Lsc      = 14.1e3#1.0 / ( ρg / (min(η_bg, G_r*tsc) * ε̇_bg) )
    vsc      = Lsc / tsc

    # Dependent physics
    Lx       = Lsc * Lx_Lc
    Ly       = Lx / Lx_Ly
    Lx_r     = 28.2
    ρgx      = ρ_r * gx
    ρgy      = ρ_r * gy
    max_LxLy = max(Lx, Ly)

    # Numerics
    nx, ny   = 4*16-2, 3*16-2
    nt       = 40
    ϵ_tol    = 5e-8
    max_iter = 5e5
    ncheck   = 500
    CFL      = 0.6 / sqrt(2.0) / 10.0
    dx, dy   = Lx / nx, Ly / ny
    _dx, _dy = 1.0 / dx, 1.0 / dy
    Ṽ        = min(dx,dy) * CFL
    Re       = 3.0 * sqrt(5.0) / 2.0 * π
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
    τxx      = zeros(Float64, nx    , ny    ) .- 0.5 .* pconf
    τyy      = zeros(Float64, nx    , ny    ) .+ 0.5 .* pconf
    τzz      = zeros(Float64, nx    , ny    )
    τxy      = zeros(Float64, nx    , ny    )
    τxyv     = zeros(Float64, nx + 1, ny + 1)
    τxx1     = zeros(Float64, nx    , ny    )
    τyy1     = zeros(Float64, nx    , ny    )
    τzz1     = zeros(Float64, nx    , ny    )
    τxy1     = zeros(Float64, nx    , ny    )
    τII      = zeros(Float64, nx    , ny    )
    ε̇xx      = zeros(Float64, nx    , ny    )
    ε̇yy      = zeros(Float64, nx    , ny    )
    ε̇zz      = zeros(Float64, nx    , ny    )
    ε̇xy      = zeros(Float64, nx    , ny    )
    ε̇xyv     = zeros(Float64, nx + 1, ny + 1)
    ε̇xx1     = zeros(Float64, nx    , ny    )
    ε̇yy1     = zeros(Float64, nx    , ny    )
    ε̇zz1     = zeros(Float64, nx    , ny    )
    ε̇xy1     = zeros(Float64, nx    , ny    )
    ε̇xyv1    = zeros(Float64, nx + 1, ny + 1)
    ε̇II      = zeros(Float64, nx    , ny    )
    ε̇II_pl   = zeros(Float64, nx    , ny    )
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
    pcorr    = zeros(Float64, nx    , ny    ) .+ pconf
    ∇v       = zeros(Float64, nx    , ny    )
    ηs       = zeros(Float64, nx    , ny    ) .+ η_bg
    ηsv      = zeros(Float64, nx+1  , ny+1  ) .+ η_bg
    ηve      = zeros(Float64, nx    , ny    )
    ηe       = zeros(Float64, nx    , ny    ) .+ G_r .* tsc
    ηvev     = zeros(Float64, nx + 1, ny + 1)
    ηvep     = zeros(Float64, nx    , ny    )
    ρ        = zeros(Float64, nx    , ny    ) .+ ρ_r
    G        = zeros(Float64, nx    , ny    ) .+ G_r
    Gv       = zeros(Float64, nx + 1, ny + 1) .+ G_r
    K        = zeros(Float64, nx    , ny    ) .+ K_r
    Hc       = zeros(Float64, nx    , ny    ) .+ H_r
    Cc       = zeros(Float64, nx    , ny    ) .+ C_r
    Pla      = zeros(Float64, nx    , ny    )
    Fc       = zeros(Float64, nx    , ny    )
    Fcrel    = zeros(Float64, nx    , ny    )
    λ̇c       = zeros(Float64, nx    , ny    )
    λ̇crel    = zeros(Float64, nx    , ny    )
    λ̇c_old   = zeros(Float64, nx    , ny    )
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
    x_inc    = 0.0
    y_inc    = 0.0
    r        = Lx / Lx_r
    x2D      = reshape(repeat(xc , ny), nx, ny)
    y2D      = reshape(repeat(yc', nx), nx, ny)

    # Old fields
    τxx_old    = deepcopy(τxx )
    τyy_old    = deepcopy(τyy )
    τzz_old    = deepcopy(τzz )
    τxy_old    = deepcopy(τxy )
    τxyv_old   = deepcopy(τxyv)
    p_old      = deepcopy(p   )
    Cc_old     = deepcopy(Cc  )
    ε̇II_pl_old = deepcopy(ε̇II_pl)

    # Boundary conditions
    for idy in axes(vx, 2)
        if idy > 1 && idy < ny + 2
            vx[:, idy] .= -ε̇_bg .* (xv .- 0.5 .* Lx)
        end
    end
    for idx in axes(vy, 1)
        if idx > 1 && idx < nx + 2
            vy[idx, :] .=  ε̇_bg .* (yv .- 1.0 .* Ly)
        end
    end
    set_BC_velocity!(vx, vy, collect(xv), collect(yv), Lx, Ly, ε̇_bg, "ε̇_bg const.")

    # Set circular inclusion
    for idy in eachindex(yc)
        for idx in eachindex(xc)
            if ( (xc[idx] - x_inc)^2 + (yc[idy] - y_inc)^2 < r^2 )
                # ηs[idx, idy] = ηs_inc
                G[idx, idy]  = G_r / 2.0
            end
        end
    end
    for idy in eachindex(yv)
        for idx in eachindex(xv)
            if ( (xv[idx] - x_inc)^2 + (yv[idy] - y_inc)^2 < r^2 )
                # ηs[idx, idy] = ηs_inc
                Gv[idx, idy]  = G_r / 2.0
            end
        end
    end

    ηs2   = deepcopy(ηs)
    G2    = deepcopy(G)
    sinφ2 = deepcopy(sinφ)
    cosφ2 = deepcopy(cosφ)
    sinψ2 = deepcopy(sinψ)

    # Smooth out the inclusion
    # smooth_2DArray_diffusion!(ηs, ηs2,     1)
    # smooth_2DArray_diffusion!(G, G2,     1)
    # smooth_2DArray_diffusion!(sinφ, sinφ2, 1)
    # smooth_2DArray_diffusion!(cosφ, cosφ2, 1)
    # smooth_2DArray_diffusion!(sinψ, sinψ2, 1)

    # Set weights
    set_weigths!(wNW, wNE, wSW, wSE)

    # Visualize
    fg1 = Figure(size = (1000, 600 / Lx_Ly), fontsize = 10, figure_padding = 10)
    gl1 = GridLayout(fg1[1,1])
    gl2 = GridLayout(fg1[1,2])
    gl3 = GridLayout(fg1[1,3])
    gl4 = GridLayout(fg1[2,1])
    gl5 = GridLayout(fg1[2,2])
    gl6 = GridLayout(fg1[2,3])
    ax1 = Axis(gl1[1, 1], xlabel = "x [m]",  ylabel = "y [m]",    aspect = DataAspect(), title = "Pt [Pa]"     )
    ax2 = Axis(gl2[1, 1], xlabel = "nt",     ylabel = "τII [Pa]", aspect = Lx_Ly,        title = "av(τII) [Pa]", ytickformat = "{:.2e}")
    ax3 = Axis(gl3[1, 1], xlabel = "x [m]",  ylabel = "y [m]",    aspect = DataAspect(), title = "τII [Pa]"    )
    ax4 = Axis(gl4[1, 1], xlabel = "x [m]",  ylabel = "y [m]",    aspect = DataAspect(), title = "C [Pa]"     )
    ax5 = Axis(gl5[1, 1], xlabel = "x [m]",  ylabel = "y [m]",    aspect = DataAspect(), title = "vx [m/s]"     )
    ax6 = Axis(gl6[1, 1], xlabel = "x [m]",  ylabel = "y [m]",    aspect = DataAspect(), title = "vy [m/s]"     )
    hm1 = heatmap!(ax1, xc, yc, p, colormap=:acton, colorrange = (pconf - pconf/2, pconf + pconf / 2))
    hm2 = heatmap!(ax3, xc, yc, τII, colormap=:hot, colorrange = (pconf - pconf/2, pconf + pconf / 2) )
    hm3 = heatmap!(ax4, xc, yc, Cc, colormap=:hot, colorrange = (C_r / 2.0, C_r) )
    hm4 = heatmap!(ax5, xv, yv, vx[:, 2:end-1],                     colormap=:roma )
    hm5 = heatmap!(ax6, xv, yv, vy[2:end-1, :],                     colormap=:roma )
    sc1 = scatter!(ax2, [0], [0], color = :black)
    cb1 = Colorbar(gl1[1, 2], hm1); cb1.tickformat = "{:.2e}"
    cb2 = Colorbar(gl3[1, 2], hm2); cb2.tickformat = "{:.2e}"
    cb3 = Colorbar(gl4[1, 2], hm3); cb3.tickformat = "{:.2e}"
    cb4 = Colorbar(gl5[1, 2], hm4); cb4.tickformat = "{:.2e}"
    cb5 = Colorbar(gl6[1, 2], hm5); cb5.tickformat = "{:.2e}"
    rowsize!(gl1, 1, Aspect(1, 1.0)); rowgap!(gl1, 10.0); colgap!(gl1, 10.0)
    rowsize!(gl2, 1, Aspect(1, 1.0)); rowgap!(gl2, 10.0); colgap!(gl2, 10.0)
    rowsize!(gl3, 1, Aspect(1, 1.0)); rowgap!(gl3, 10.0); colgap!(gl3, 10.0)
    rowsize!(gl4, 1, Aspect(1, 1.0)); rowgap!(gl4, 10.0); colgap!(gl4, 10.0)
    rowsize!(gl5, 1, Aspect(1, 1.0)); rowgap!(gl5, 10.0); colgap!(gl5, 10.0)
    rowsize!(gl6, 1, Aspect(1, 1.0)); rowgap!(gl6, 10.0); colgap!(gl6, 10.0)
    display(fg1)

    # TIME LOOP
    time = 0.0; τII_max_time = []; ntimes = []; τII_mean_time = []
    for idx_Time in 1:nt

        # Update time(
        # dt   = min(η_bg / G_r / 10.0, min(dx, dy) / max(maximum(abs.(vx)), maximum(abs.(vy))) / 4.1)
        # dt   = 5.0*min(η_bg, minimum(G.*tsc)) / maximum(G) / 10.0
        dt = tsc
        time = time + dt

        @printf("\ndt = %.5e\n", dt)

        # Update old physical fields
        τxx_old    .= τxx
        τyy_old    .= τyy
        τzz_old    .= τzz
        τxyv_old   .= τxyv
        τxy_old    .= @av_arr(τxyv)
        p_old      .= p
        Cc_old     .= Cc
        ε̇II_pl_old .= ε̇II_pl

        # Pseudo - transient solver loop
        max_Res_all = []; err = 2ϵ_tol
        λ̇crel .= 0.0
        idx_Time == 1 ? ηvep .= ηve : nothing
        # c2v(ηve, ηvev, AvNW,AvNE,AvSW, AvSE, wNW, wNE, wSW, wSE)
        for iter in 1:max_iter
            
            λ̇c_old .= λ̇crel
            compute_viscoelastic_rheology!(ηve, ηs, G.*dt)
            compute_viscoelastic_rheology!(ηvev, ηsv, Gv.*dt)
            compute_pseudotransient_parameters!(ηve, dτ_ρ, G̃dτ, Ṽ, max_LxLy, Re, r̃)
            compute_divergence_vs!(∇v, vx, vy, _dx, _dy)
            compute_vorticity!(ωxx, ωyy, ωxyv, ωyxv, vx, vy, _dx, _dy)
            v2c(ωxyv, ωxy)
            v2c(ωyxv, ωyx)
            compute_Jaumann_derivative!(Jxx, Jyy, Jxy, Jxyv, ωxx, ωyy, ωxy, ωyx, ωxyv, ωyxv, τxx, τyy, τxy, τxyv, τxx_old, τyy_old, τxy_old, τxyv_old, vx, vy, dt, _dx, _dy)
            compute_strainrates!(ε̇xx, ε̇yy, ε̇zz, ε̇xyv, vx, vy, ∇v, _dx, _dy)
            v2c(ε̇xyv, ε̇xy)
            update_stresses!(τxx, τyy, τzz, τxy, τxyv, ε̇xx, ε̇yy, ε̇zz, ε̇xy, ε̇xyv, Jxx, Jyy, Jzz, Jxy, Jxyv, ∇v, vx, vy, ηve, ηvev, G, Gv, G̃dτ, G̃dτv, _dx, _dy)
            c2v(Jxy, Jxyv, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
            compute_trial_stresses!(ε̇xx, ε̇yy, ε̇zz, ε̇xy, ε̇xyv, ε̇xx1, ε̇yy1, ε̇zz1, ε̇xy1, ε̇xyv1, ε̇II, τxx1, τyy1, τzz1, τxy1, ηve, G, Gv, Jxx, Jyy, Jzz, Jxy, Jxyv)
            compute_2ndInv!(τII, τxx1, τyy1, τzz1, τxy1)
            check_yield_1!(τII, p, Cc_old, Fcrel, sinφ, cosφ)
            if (iter % ncheck) == 0 max_Fcrel = maximum(Fcrel);  @printf("Ini Fc_rel = %1.3e \n", max_Fcrel / Psc) end
            compute_plastic_flow_potential!(Pla, Fcrel, Hc, K, ηve, sinφ, sinψ, cosφ, τII, τxx1, τyy1, τzz1, τxy1, λ̇c, λ̇crel, λ̇c_old, ∂Q∂τxx, ∂Q∂τyy,∂Q∂τzz, ∂Q∂τxy, ηvp, dt, rel)
            if iter % ncheck == 0
                correct_stress_pressure!(K, τxx, τyy, τzz, τxy, ηve, ε̇xx1, ε̇yy1, ε̇zz1, ε̇xy1, λ̇c, ∂Q∂τxx, ∂Q∂τyy, ∂Q∂τzz, ∂Q∂τxy, sinψ, pcorr, p, Hc, Cc, Cc_old, dt)
                compute_2ndInv!(τII, τxx, τyy, τzz, τxy)
                check_yield_2!(τII, Fc, pcorr, sinφ, cosφ, Cc, λ̇c, ηvp)
                max_Fc = maximum(Fc);  @printf("Check Fc_phys = %1.3e \n", max_Fc / Psc)
            end
            correct_stress_pressure!(K, τxx, τyy, τzz, τxy, ηve, ε̇xx1, ε̇yy1, ε̇zz1, ε̇xy1, λ̇crel, ∂Q∂τxx, ∂Q∂τyy, ∂Q∂τzz, ∂Q∂τxy, sinψ, pcorr, p, Hc, Cc, Cc_old, dt)
            compute_2ndInv!(τII, τxx, τyy, τzz, τxy)
            check_yield_2!(τII, Fcrel, pcorr, sinφ, cosφ, Cc, λ̇crel, ηvp)
            if (iter % ncheck) == 0 max_Fcrel = maximum(Fcrel);  @printf("Check Fc_rel = %1.3e \n", max_Fcrel / Psc) end
            compute_viscoelastoplastic_rheology!(ηvep, τII, ε̇II)
            interpolate_shearstress_vertices!(τxy, τxyv)
            compute_pseudotransient_parameters!(ηvep, dτ_ρ, G̃dτ, Ṽ, max_LxLy, Re, r̃)
            compute_pressure_change!(∇v, dpdτ, p, p_old, K, G̃dτ, dt, r̃)
            compute_velocity_change!(τxx, τyy, τxyv, pcorr, dvxdτ, dvydτ, dτ_ρ, _dx, _dy)
            update_pressure!(p, dpdτ)
            update_velocity!(vx, vy, dvxdτ, dvydτ)

            # Compute density
            compute_density!(ρ, p, K, ρ0, p0)

            # Residuals
            compute_pressure_residual!(pcorr, p_old, K, ∇v, Res_p, dt)
            compute_velocity_residual!(τxx, τyy, τxyv, pcorr, Res_vx, Res_vy, _dx, _dy)
            
            # Boundary conditions
            set_BC_velocity!(vx, vy, collect(xv), collect(yv), Lx, Ly, ε̇_bg, "ε̇_bg const.")

            # Monitor residuals 
            if iter % ncheck == 0
                @printf("\n")
                @printf("Time step = %d | time = %.2e \t Iteration = %d \n", idx_Time, time, iter)
                @printf("Residual p  = %.2e \n", (maximum(Res_p ) * dt) )
                @printf("Residual vx = %.2e \n", (maximum(Res_vx) / Psc * Lsc))
                @printf("Residual vy = %.2e \n", (maximum(Res_vy) / Psc * Lsc))
                @printf("Physics:\n")
                @printf "min(vx   ) = %.2e \t max(vx   ) = %.2e\n" minimum(vx    ) maximum(vx    )
                @printf "min(vy   ) = %.2e \t max(vy   ) = %.2e\n" minimum(vy    ) maximum(vy    )
                @printf "min(ηvep ) = %.2e \t max(ηvep ) = %.2e\n" minimum(ηvep  ) maximum(ηvep  )
                @printf "min(ηve  ) = %.2e \t max(ηve  ) = %.2e\n" minimum(ηve   ) maximum(ηve   )
                @printf "min(ηe   ) = %.2e \t max(ηe   ) = %.2e\n" minimum(G.*dt ) maximum(G.*dt )
                @printf "min(G    ) = %.2e \t max(G    ) = %.2e\n" minimum(G     ) maximum(G     )
                @printf "min(λ̇c   ) = %.2e \t max(λ̇c   ) = %.2e\n" minimum(λ̇c    ) maximum(λ̇c    )
                @printf "min(λ̇crel) = %.2e \t max(λ̇crel) = %.2e\n" minimum(λ̇crel ) maximum(λ̇crel )
                @printf "min(Cc   ) = %.2e \t max(Cc   ) = %.2e\n" minimum(Cc    ) maximum(Cc    )
                @printf "min(Fcrel) = %.2e \t max(Fcrel) = %.2e\n" minimum(Fcrel ) maximum(Fcrel )
                @printf "min(Fc   ) = %.2e \t max(Fc   ) = %.2e\n" minimum(Fc    ) maximum(Fc    )
                @printf "min(Pla  ) = %.2e \t max(Pla  ) = %.2e\n" minimum(Pla   ) maximum(Pla   )
                @printf "min(Hc   ) = %.2e \t max(Hc   ) = %.2e\n" minimum(Hc    ) maximum(Hc    )
                @printf "min(τxx  ) = %.2e \t max(τxx  ) = %.2e\n" minimum(τxx   ) maximum(τxx   )
                @printf "min(τyy  ) = %.2e \t max(τyy  ) = %.2e\n" minimum(τyy   ) maximum(τyy   )
                @printf "min(τzz  ) = %.2e \t max(τzz  ) = %.2e\n" minimum(τzz   ) maximum(τzz   )
                @printf "min(τxy  ) = %.2e \t max(τxy  ) = %.2e\n" minimum(τxy   ) maximum(τxy   )
                @printf "min(PsinC) = %.2e \t max(PsinC) = %.2e\n" minimum(p.*sinφ.-Cc_old.*cosφ ) maximum(p.*sinφ.-Cc_old.*cosφ )
                @printf "min(τII  ) = %.2e \t max(τII  ) = %.2e\n" minimum(τII   ) maximum(τII   )
                @printf "min(ε̇II  ) = %.2e \t max(ε̇II  ) = %.2e\n" minimum(ε̇II   ) maximum(ε̇II   )
                @printf "min(ε̇IIp ) = %.2e \t max(ε̇IIp ) = %.2e\n" minimum(ε̇II_pl) maximum(ε̇II_pl)
                err = maximum([maximum(Res_p ) * dt maximum(Res_vx) / Psc * Lsc maximum(Res_vy) / Psc * Lsc])
            end

            if err < ϵ_tol
                break
            end
        end
    
        # Store maximum shear stress to monitor build up
        push!(τII_max_time, maximum(τII))
        push!(τII_mean_time,   mean(τII))
        push!(ntimes, idx_Time)

        # Include plastic correction in converged fields
        update_fields_plasticity!(p, pcorr, ε̇II_pl, ε̇II_pl_old, λ̇crel, Hc, Cc, Cc_old, dt, C_r)

        # Update figure
        data     = [[ntimes[i], τII_mean_time[i]] for i in eachindex(ntimes)]
        hm1[3][] = p
        hm2[3][] = τII
        hm3[3][] = Cc
        hm4[3][] = vx
        hm5[3][] = vy
        sc1[1][] = data
        hm1.colorrange = (minimum(p), maximum(p))
        hm2.colorrange = (minimum(τII), maximum(τII))
        if minimum(Cc ) != maximum(Cc )
            hm3.colorrange = (minimum(Cc ), maximum(Cc ))
        end
        hm4.colorrange = (minimum(vx), maximum(vx))
        hm5.colorrange = (minimum(vy), maximum(vy))
        display(fg1)

        # Save image
        if save_figure
            fname = @sprintf("./png/Stokes2D_vevp_%05d.png", idx_Time)
            save(fname, fg1, px_per_unit = 1.5)
        end
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
function smooth_2DArray_diffusion!(A :: Matrix{Float64}, A2 :: Matrix{Float64}, nsteps :: Int64)
    nx, ny = size(A)
    for _ in 1:nsteps
        for idy in 1:ny
            for idx in 1:nx
                if idx > 1 && idx < nx && idy > 1 && idy < ny
                    A2[idx, idy] = A[idx, idy] + (1.0 / 4.1 / 5.0) * ( @d2_dx2(A) + @d2_dy2(A) )
                end
            end
        end
        A2[[1, end], :] .= A2[[2, end-1], :]
        A2[[1, end], :] .= A2[[2, end-1], :]
        A2[:, [1, end]] .= A2[:, [2, end-1]]
        A2[:, [1, end]] .= A2[:, [2, end-1]]
        A .= A2
    end

    # Return
    return nothing

end

# Compute sqrt of 2nd tensor invariant
@views function compute_2ndInv!(
    TII  :: Matrix{Float64},
    Txx  :: Matrix{Float64},
    Tyy  :: Matrix{Float64},
    Tzz  :: Matrix{Float64},
    Txy  :: Matrix{Float64},
)
    # Get size
    nx, ny = size(TII)

    # Compute sqrt of 2nd invariant
    for idy in 1:ny
        for idx in 1:nx
            TII[idx, idy] = sqrt( 0.5 * (Txx[idx, idy]^2 + Tyy[idx, idy]^2 + Tzz[idx, idy]^2 + 2.0 * Txy[idx, idy]^2))
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
    ηe  :: Matrix{Float64}
)
    # Get array sizes
    nx, ny = size(ηve)

    # Viscoelastic rheology
    for idy in 1:ny
        for idx in 1:nx
            ηve[idx, idy] = 1.0 / (1.0 / (ηe[idx, idy]) + 1.0 / (ηs[idx, idy]))
        end
    end

    # Return
    return nothing

end

# Trial stresses and strain rates
function compute_trial_stresses!(
    ε̇xx   :: Matrix{Float64},
    ε̇yy   :: Matrix{Float64},
    ε̇zz   :: Matrix{Float64},
    ε̇xy   :: Matrix{Float64},
    ε̇xyv  :: Matrix{Float64},
    ε̇xx1  :: Matrix{Float64},
    ε̇yy1  :: Matrix{Float64},
    ε̇zz1  :: Matrix{Float64},
    ε̇xy1  :: Matrix{Float64},
    ε̇xyv1 :: Matrix{Float64},
    ε̇II   :: Matrix{Float64},
    τxx1  :: Matrix{Float64},
    τyy1  :: Matrix{Float64},
    τzz1  :: Matrix{Float64},
    τxy1  :: Matrix{Float64},
    ηve   :: Matrix{Float64},
    G     :: Matrix{Float64},
    Gv    :: Matrix{Float64},
    Jxx   :: Matrix{Float64},
    Jyy   :: Matrix{Float64},
    Jzz   :: Matrix{Float64},
    Jxy   :: Matrix{Float64},
    Jxyv  :: Matrix{Float64}
)
    # Get sizes
    nxc, nyc = size(G)
    nxv, nyv = size(Gv)

    # Compute cell centers
    for idy in 1:nyc
        for idx in 1:nxc
            ε̇xx1[idx, idy]  = ε̇xx[idx, idy] - Jxx[idx, idy] / 2.0 / G[idx, idy]
            ε̇yy1[idx, idy]  = ε̇yy[idx, idy] - Jyy[idx, idy] / 2.0 / G[idx, idy]
            ε̇zz1[idx, idy]  = ε̇zz[idx, idy] - Jzz[idx, idy] / 2.0 / G[idx, idy]
            ε̇xy1[idx, idy]  = ε̇xy[idx, idy] - Jxy[idx, idy] / 2.0 / G[idx, idy]
            ε̇II[idx, idy]   = sqrt(0.5 * (ε̇xx1[idx, idy]^2 + ε̇yy1[idx, idy]^2 + ε̇zz1[idx, idy]^2 + 2.0 * (ε̇xy1[idx, idy])^2))
            τxx1[idx, idy] = 2.0 * ηve[idx, idy] * ε̇xx1[idx, idy]
            τyy1[idx, idy] = 2.0 * ηve[idx, idy] * ε̇yy1[idx, idy]
            τzz1[idx, idy] = 2.0 * ηve[idx, idy] * ε̇zz1[idx, idy]
            τxy1[idx, idy] = 2.0 * ηve[idx, idy] * ε̇xy1[idx, idy]
        end
    end

    # Compute vertices
    for idy in 1:nyv
        for idx in 1:nxv
            ε̇xyv1[idx, idy] = ε̇xyv[idx, idy] - Jxyv[idx, idy] / 2.0 / Gv[idx, idy]
        end
    end

    # Return
    return nothing
end

# Check yielding - part I
function check_yield_1!(
    τII    :: Matrix{Float64},
    p      :: Matrix{Float64},
    Cc_old :: Matrix{Float64},
    F      :: Matrix{Float64},
    sinφ   :: Matrix{Float64},
    cosφ   :: Matrix{Float64}
)
    # Get sizes
    nx, ny = size(p)

    # Evaluate yield
    for idy in 1:ny
        for idx in 1:nx
            F[idx, idy] = τII[idx, idy] - p[idx, idy] * sinφ[idx, idy] - Cc_old[idx, idy] * cosφ[idx, idy]
        end
    end
end

# Plastic multiplier and flow potential
function compute_plastic_flow_potential!(
    Pla    :: Matrix{Float64},
    F      :: Matrix{Float64},
    Hc     :: Matrix{Float64},
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
    λ̇c      :: Matrix{Float64},
    λ̇crel   :: Matrix{Float64},
    λ̇c_old  :: Matrix{Float64},
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
            λ̇c[idx, idy] = Pla[idx, idy] * (F[idx, idy] / (ηve[idx, idy] + ηvp + K[idx, idy]*dt*sinφ[idx, idy]*sinψ[idx, idy] + Hc[idx, idy]*cosφ[idx, idy]*dt))
            λ̇crel[idx, idy] = rel * Pla[idx, idy] * (F[idx, idy] / (ηve[idx, idy] + ηvp + K[idx, idy]*dt*sinφ[idx, idy]*sinψ[idx, idy] + Hc[idx, idy]*cosφ[idx, idy]*dt)) + (1.0 - rel) * λ̇c_old[idx, idy]
        end
    end

    # Derivative of flow potential
    for idy in 1:ny
        for idx in 1:nx
            ∂Q∂τxx[idx, idy] = 0.5 / τII[idx, idy] * τxx[idx, idy]
            ∂Q∂τyy[idx, idy] = 0.5 / τII[idx, idy] * τyy[idx, idy]
            ∂Q∂τzz[idx, idy] = 0.5 / τII[idx, idy] * τzz[idx, idy]
            ∂Q∂τxy[idx, idy] = 1.0 / τII[idx, idy] * τxy[idx, idy]
        end
    end

end

# Plastic corrections
function correct_stress_pressure!(
    K       :: Matrix{Float64},
    τxx     :: Matrix{Float64},
    τyy     :: Matrix{Float64},
    τzz     :: Matrix{Float64},
    τxy     :: Matrix{Float64},
    ηve     :: Matrix{Float64},
    ε̇xx1    :: Matrix{Float64},
    ε̇yy1    :: Matrix{Float64},
    ε̇zz1    :: Matrix{Float64},
    ε̇xy1    :: Matrix{Float64},
    λ̇c      :: Matrix{Float64},
    ∂Q∂τxx  :: Matrix{Float64},
    ∂Q∂τyy  :: Matrix{Float64},
    ∂Q∂τzz  :: Matrix{Float64},
    ∂Q∂τxy  :: Matrix{Float64},
    sinψ    :: Matrix{Float64},
    pcorr   :: Matrix{Float64},
    p       :: Matrix{Float64},
    Hc      :: Matrix{Float64},
    Cc      :: Matrix{Float64},
    Cc_old  :: Matrix{Float64},
    dt      :: Float64
)
    # Get sizes
    nx, ny = size(τxx)

    # Corrections
    for idy in 1:ny
        for idx in 1:nx
            τxx[idx, idy]   = 2.0 * ηve[idx, idy] * (ε̇xx1[idx, idy] -       λ̇c[idx, idy] * ∂Q∂τxx[idx, idy])
            τyy[idx, idy]   = 2.0 * ηve[idx, idy] * (ε̇yy1[idx, idy] -       λ̇c[idx, idy] * ∂Q∂τyy[idx, idy])
            τzz[idx, idy]   = 2.0 * ηve[idx, idy] * (ε̇zz1[idx, idy] -       λ̇c[idx, idy] * ∂Q∂τzz[idx, idy])
            τxy[idx, idy]   = 2.0 * ηve[idx, idy] * (ε̇xy1[idx, idy] - 0.5 * λ̇c[idx, idy] * ∂Q∂τxy[idx, idy])
            pcorr[idx, idy] = p[idx, idy] + K[idx, idy] * dt * λ̇c[idx, idy] * sinψ[idx, idy]
            Cc[idx, idy]    = Cc_old[idx, idy] + dt * sqrt(2/3) * λ̇c[idx, idy] * Hc[idx, idy]
        end
    end
end

# Check yield - part II
function check_yield_2!(
    τII   :: Matrix{Float64},
    F     :: Matrix{Float64},
    pcorr :: Matrix{Float64},
    sinφ  :: Matrix{Float64},
    cosφ  :: Matrix{Float64},
    Cc    :: Matrix{Float64},
    λ̇c    :: Matrix{Float64},
    ηvp   :: Float64
)
    # Get sizes
    nx, ny = size(τII)

    # Check yield
    for idy in 1:ny
        for idx in 1:nx
            F[idx, idy] = τII[idx, idy] - pcorr[idx, idy] * sinφ[idx, idy] - Cc[idx, idy] * cosφ[idx, idy] - ηvp * λ̇c[idx, idy]
        end
    end
end

# Include plastic corrections
function update_fields_plasticity!(
    p          :: Matrix{Float64},
    pcorr      :: Matrix{Float64},
    ε̇II_pl     :: Matrix{Float64},
    ε̇II_pl_old :: Matrix{Float64},
    λ̇c         :: Matrix{Float64},
    Hc         :: Matrix{Float64},
    Cc         :: Matrix{Float64},
    Cc_old     :: Matrix{Float64},
    dt         :: Float64,
    C_r        :: Float64
)
    # Get sizes
    ncx, ncy = size(p)

    # Include corrections
    for idy in 1:ncy
        for idx in 1:ncx
            p[idx, idy]      = pcorr[idx, idy]
            ε̇II_pl[idx, idy] = ε̇II_pl_old[idx, idy] + dt*sqrt(2.0/3.0)*λ̇c[idx, idy]
            Cc[idx, idy]     = Cc_old[idx, idy]     + dt*sqrt(2.0/3.0)*λ̇c[idx, idy]*Hc[idx, idy]
            # Limit softening
            Cc[idx, idy] < C_r/2.0 ? Hc[idx, idy] = 0.0 : nothing
        end
    end
    # Return
    return nothing
end

# Visco-elasto-plastic rheology
function compute_viscoelastoplastic_rheology!(
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
    Jxx      :: Matrix{Float64},
    Jyy      :: Matrix{Float64},
    Jxy      :: Matrix{Float64},
    Jxyv     :: Matrix{Float64},
    ωxx      :: Matrix{Float64},
    ωyy      :: Matrix{Float64},
    ωxy      :: Matrix{Float64},
    ωyx      :: Matrix{Float64},
    ωxyv     :: Matrix{Float64},
    ωyxv     :: Matrix{Float64},
    τxx      :: Matrix{Float64},
    τyy      :: Matrix{Float64},
    τxy      :: Matrix{Float64},
    τxyv     :: Matrix{Float64},
    τxx_old  :: Matrix{Float64},
    τyy_old  :: Matrix{Float64},
    τxy_old  :: Matrix{Float64},
    τxyv_old :: Matrix{Float64},
    vx       :: Matrix{Float64},
    vy       :: Matrix{Float64},
    dt       :: Float64,
    _dx      :: Float64,
    _dy      :: Float64
)
    # Get sizes
    nxc, nyc = size(Jxx)

    # Compute components
    for idy in 1:nyc
        for idx in 1:nxc

            # Old values
            Jxx[idx, idy] = -τxx_old[idx, idy] / dt
            Jyy[idx, idy] = -τyy_old[idx, idy] / dt
            Jxy[idx, idy] = -τxy_old[idx, idy] / dt

            # Advection
            if idx > 1 && idx < nxc && idy > 1 && idy < nyc
                Jxx[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τxx) * _dx
                Jyy[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τyy) * _dx
                Jxx[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τxx) * _dy
                Jyy[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τyy) * _dy
                Jxy[idx, idy] += max(0.0, vx[idx, idy]) * @d_dx(τxy) * _dx
                Jxy[idx, idy] += max(0.0, vy[idx, idy]) * @d_dy(τxy) * _dy
            end
            if idx < nxc && idy < nyc
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

    # Return
    return nothing
end

# Strain rates
function compute_strainrates!(
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
    τxx  :: Matrix{Float64},
    τyy  :: Matrix{Float64},
    τzz  :: Matrix{Float64},
    τxy  :: Matrix{Float64},
    τxyv :: Matrix{Float64},
    ε̇xx  :: Matrix{Float64},
    ε̇yy  :: Matrix{Float64},
    ε̇zz  :: Matrix{Float64},
    ε̇xy  :: Matrix{Float64},
    ε̇xyv :: Matrix{Float64},
    Jxx  :: Matrix{Float64},
    Jyy  :: Matrix{Float64},
    Jzz  :: Matrix{Float64},
    Jxy  :: Matrix{Float64},
    Jxyv :: Matrix{Float64},
    ∇v   :: Matrix{Float64},
    vx   :: Matrix{Float64},
    vy   :: Matrix{Float64},
    ηve  :: Matrix{Float64},
    ηvev :: Matrix{Float64},
    G    :: Matrix{Float64},
    Gv   :: Matrix{Float64},
    G̃dτ  :: Matrix{Float64},
    G̃dτv :: Matrix{Float64},
    _dx  :: Float64,
    _dy  :: Float64
)
    # Get sizes
    nx, ny = size(G)

    # Compute components
    for idy in 1:ny
        for idx in 1:nx
            τxx[idx, idy] = ( τxx[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jxx[idx, idy] + ε̇xx[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            τyy[idx, idy] = ( τyy[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jyy[idx, idy] + ε̇yy[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            τzz[idx, idy] = ( τzz[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jzz[idx, idy] - 1.0/3.0 * ∇v[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
            τxy[idx, idy] = ( τxy[idx, idy] + 2.0 * G̃dτ[idx, idy] * ( -1.0 / (2.0 * G[idx, idy]) * Jxy[idx, idy] + ε̇xy[idx, idy]) ) / (G̃dτ[idx, idy] / ηve[idx, idy] + 1.0)
        end
    end

    # Return
    return nothing
end

# Shear stress on vertices
function interpolate_shearstress_vertices!(
    τxy  :: Matrix{Float64},
    τxyv :: Matrix{Float64}
)
    # Get size
    nx, ny = size(τxy)

    # Interpolate
    for idy in 1:ny
        for idx in 1:nx
            if idx < nx && idy < ny
                τxyv[idx+1, idy+1] = @av(τxy)
            end
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
                dvxdτ[idx, idy] = @av_xa(dτ_ρ) * ( @d_dx(τxx) * _dx + @d_dyi(τxyv) * _dy - @d_dx(p) * _dx ) + (1.0 - 4 / nx) * dvxdτ[idx, idy]
            end
        end
    end
    for idy in 1:ny
        for idx in 1:nx
            if idy < ny
                dvydτ[idx, idy] = @av_ya(dτ_ρ) * ( @d_dy(τyy) * _dy + @d_dxi(τxyv) * _dx - @d_dy(p) * _dy ) + (1.0 - 4 / ny) * dvydτ[idx, idy]
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
    nx_vx, ny_vx = size(dvxdτ)[1], size(dvxdτ)[2]
    nx_vy, ny_vy = size(dvydτ)[1], size(dvydτ)[2]

    # @printf "nx_vx = %d \t ny_vx = %d \t nx_vy = %d \t ny_vy = %d\n" nx_vx nx_vy nx_vy ny_vy
    # error("Hello...")
    # Update velocity
    for idy in 1:ny_vx
        for idx in 1:nx_vx
            vx[idx + 1, idy + 1] += dvxdτ[idx, idy]
         end
    end
    for idy in 1:ny_vy
        for idx in 1:nx_vy
            vy[idx + 1, idy + 1] += dvydτ[idx, idy]
         end
    end

    # Return
    return nothing
end

# Boundary conditions
@views function set_BC_velocity!(
    vx     :: Matrix,
    vy     :: Matrix,
    xv     :: Vector{Float64},
    yv     :: Vector{Float64},
    Lx     :: Float64,
    Ly     :: Float64,
    ε̇_bg   :: Float64,
    BCtype :: String
)
    # Velocity in x
    # vx[:,   2:end-1] .= -0.5 .* (xv .- Lx / 2.0) .* ε̇_bg
    vx[:,   1] .= vx[:,     2]
    vx[:, end] .= vx[:, end-1]

    # Velocity in y
    # vy[2:end-1,   :] .= (yv' .- Ly) .* ε̇_bg
    vy[1,   :] .= vy[2,     :]
    vy[end, :] .= vy[end-1, :]

    # Return
    return nothing
end

# --------------------- #
#|      Run main       |#
# --------------------- #
compressible_vevp_stokes2D();