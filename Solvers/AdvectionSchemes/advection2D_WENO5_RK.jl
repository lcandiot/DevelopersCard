# This script showcases standard 1st order upwind advection and compare Runge-Kutta methods with Euler time integration
using CairoMakie, MathTeXEngine, ColorSchemes

# Set constants
const ϵ = 1e-5
const DatType = Float64

# Compute functions
function compute_qT!(
    T   :: Matrix{DatType},
    qT  :: Matrix{DatType},
    k   :: DatType,
    _dx :: DatType
)
    # Compute
    for idx_x in eachindex(T)
        if idx_x < size(T, 1)
            qT[idx_x + 1, 1] = k * (T[idx_x + 1] - T[idx_x]) * _dx
        end
    end

    # Return
    return nothing
end

function advect_WENOZ_plus!(
    T      :: AbstractArray,
    vxT    :: Matrix{DatType},
    vyT    :: Matrix{DatType},
    qTxa   :: Matrix{DatType},
    qTya   :: Matrix{DatType},
    ∂T∂t_a :: Matrix{DatType},
    α      :: Vector{DatType},
    β      :: Vector{DatType},
    d      :: Vector{DatType},
    q      :: Matrix{DatType},
    f      :: Vector{DatType},
    w      :: Matrix{DatType},
    vx     :: Matrix{DatType},
    vy     :: Matrix{DatType},
    _dx    :: DatType,
    _dy    :: DatType,
    λ      :: DatType
)
    # Reset
    ∂T∂t_a .= 0.0
    qTxa   .= 0.0
    qTya   .= 0.0
    vxT    .= 0.0
    vyT    .= 0.0

    # Assemble the fluxes vx*T , vy*T
    for idx_y in axes(vxT,2)
        for idx_x in axes(vxT, 1)
            if idx_x > 3 && idx_x < size(vxT, 1) - 2
                vx[idx_x - 3, idx_y] > 0.0 ? vxT[idx_x, idx_y] = vx[idx_x - 3, idx_y] * T[idx_x - 1, idx_y + 3] : vxT[idx_x, idx_y] = vx[idx_x - 3, idx_y] * T[idx_x, idx_y + 3]
            end
        end
    end
    for idx_y in axes(vyT,2)
        for idx_x in axes(vyT, 1)
            if idx_y > 3 && idx_y < size(vyT, 2) - 2
                vy[idx_x, idx_y - 3] > 0.0 ? vyT[idx_x, idx_y] = vy[idx_x, idx_y - 3] * T[idx_x + 3, idx_y - 1] : vyT[idx_x, idx_y] = vy[idx_x, idx_y - 3] * T[idx_x + 3, idx_y]
            end
        end
    end

    for idx_y in axes(vxT, 2)
        for idx_x in axes(vxT, 1)
            if idx_x > 3 && idx_x < size(vxT, 1) - 2
                # Left-biased stencil - positive flow direction
                f[1], f[2], f[3], f[4], f[5] = vxT[idx_x - 3, idx_y], vxT[idx_x - 2, idx_y], vxT[idx_x - 1, idx_y], vxT[idx_x, idx_y], vxT[idx_x + 1, idx_y]

                # Flux reconstructions
                q[1, 1] = ( 2.0 * f[1] - 7.0 * f[2] + 11.0 * f[3]) / 6.0
                q[2, 1] = (-1.0 * f[2] + 5.0 * f[3] +  2.0 * f[4]) / 6.0
                q[3, 1] = ( 2.0 * f[3] + 5.0 * f[4] -  1.0 * f[5]) / 6.0

                # Smoothness indicators
                β[1] = (13.0 / 12.0) * (f[1] - 2.0 * f[2] + f[3])^2 + 0.25 * (f[1] - 4.0 * f[2] + 3.0 * f[3])^2
                β[2] = (13.0 / 12.0) * (f[2] - 2.0 * f[3] + f[4])^2 + 0.25 * (f[2] - f[4])^2
                β[3] = (13.0 / 12.0) * (f[3] - 2.0 * f[4] + f[5])^2 + 0.25 * (3.0 * f[3] - 4.0 * f[4] + f[5])^2

                # Global smoothness indicator
                τ5 = abs(β[1] - β[3])

                # Nonlinear weights
                α[1]     = d[1] * (1.0 + (τ5 / (β[1] + ϵ))^2 + λ * (β[1] + ϵ) / (τ5 + ϵ))
                α[2]     = d[2] * (1.0 + (τ5 / (β[2] + ϵ))^2 + λ * (β[2] + ϵ) / (τ5 + ϵ))
                α[3]     = d[3] * (1.0 + (τ5 / (β[3] + ϵ))^2 + λ * (β[3] + ϵ) / (τ5 + ϵ))
                w[1, :] .= (α ./ sum(α))

                # Add to advective term
                qTxa[idx_x - 3, idx_y] += (w * q)[1]

                # Right-biased stencil - negative flow direction
                f[1], f[2], f[3], f[4], f[5] = vxT[idx_x + 3, idx_y], vxT[idx_x + 2, idx_y], vxT[idx_x + 1, idx_y], vxT[idx_x, idx_y], vxT[idx_x - 1, idx_y]
                
                # Flux reconstructions
                q[1, 1] = ( 2.0 * f[1] - 7.0 * f[2] + 11.0 * f[3]) / 6.0
                q[2, 1] = (-1.0 * f[2] + 5.0 * f[3] +  2.0 * f[4]) / 6.0
                q[3, 1] = ( 2.0 * f[3] + 5.0 * f[4] -  1.0 * f[5]) / 6.0

                # Smoothness indicators
                β[1] = (13.0 / 12.0) * (f[1] - 2.0 * f[2] + f[3])^2 + 0.25 * (f[1] - 4.0 * f[2] + 3.0 * f[3])^2
                β[2] = (13.0 / 12.0) * (f[2] - 2.0 * f[3] + f[4])^2 + 0.25 * (f[2] - f[4])^2
                β[3] = (13.0 / 12.0) * (f[3] - 2.0 * f[4] + f[5])^2 + 0.25 * (3.0 * f[3] - 4.0 * f[4] + f[5])^2

                # Global smoothness indicator
                τ5 = abs(β[1] - β[3])

                # Nonlinear weights
                α[1]     = d[1] * (1.0 + (τ5 / (β[1] + ϵ))^2)
                α[2]     = d[2] * (1.0 + (τ5 / (β[2] + ϵ))^2)
                α[3]     = d[3] * (1.0 + (τ5 / (β[3] + ϵ))^2)
                w[1, :] .= (α ./ sum(α))

                # Add to advective term
                qTxa[idx_x - 2, idx_y] += (w * q)[1]
            end
        end
    end
    for idx_y in axes(vyT, 2)
        for idx_x in axes(vyT, 1)
            if idx_y > 3 && idx_y < size(vyT, 2) - 2
                # Left-biased stencil - positive flow direction
                f[1], f[2], f[3], f[4], f[5] = vyT[idx_x, idx_y - 3], vyT[idx_x, idx_y - 2], vyT[idx_x, idx_y - 1], vyT[idx_x, idx_y], vyT[idx_x, idx_y + 1]

                # Flux reconstructions
                q[1, 1] = ( 2.0 * f[1] - 7.0 * f[2] + 11.0 * f[3]) / 6.0
                q[2, 1] = (-1.0 * f[2] + 5.0 * f[3] +  2.0 * f[4]) / 6.0
                q[3, 1] = ( 2.0 * f[3] + 5.0 * f[4] -  1.0 * f[5]) / 6.0

                # Smoothness indicators
                β[1] = (13.0 / 12.0) * (f[1] - 2.0 * f[2] + f[3])^2 + 0.25 * (f[1] - 4.0 * f[2] + 3.0 * f[3])^2
                β[2] = (13.0 / 12.0) * (f[2] - 2.0 * f[3] + f[4])^2 + 0.25 * (f[2] - f[4])^2
                β[3] = (13.0 / 12.0) * (f[3] - 2.0 * f[4] + f[5])^2 + 0.25 * (3.0 * f[3] - 4.0 * f[4] + f[5])^2

                # Global smoothness indicator
                τ5 = abs(β[1] - β[3])

                # Nonlinear weights
                α[1]     = d[1] * (1.0 + (τ5 / (β[1] + ϵ))^2 + λ * (β[1] + ϵ) / (τ5 + ϵ))
                α[2]     = d[2] * (1.0 + (τ5 / (β[2] + ϵ))^2 + λ * (β[2] + ϵ) / (τ5 + ϵ))
                α[3]     = d[3] * (1.0 + (τ5 / (β[3] + ϵ))^2 + λ * (β[3] + ϵ) / (τ5 + ϵ))
                w[1, :] .= (α ./ sum(α))

                # Add to advective term
                qTya[idx_x, idx_y - 3] += (w * q)[1]

                # Right-biased stencil - negative flow direction
                f[1], f[2], f[3], f[4], f[5] = vyT[idx_x, idx_y + 3], vyT[idx_x, idx_y + 2], vyT[idx_x, idx_y + 1], vyT[idx_x, idx_y], vyT[idx_x, idx_y - 1]
                
                # Flux reconstructions
                q[1, 1] = ( 2.0 * f[1] - 7.0 * f[2] + 11.0 * f[3]) / 6.0
                q[2, 1] = (-1.0 * f[2] + 5.0 * f[3] +  2.0 * f[4]) / 6.0
                q[3, 1] = ( 2.0 * f[3] + 5.0 * f[4] -  1.0 * f[5]) / 6.0

                # Smoothness indicators
                β[1] = (13.0 / 12.0) * (f[1] - 2.0 * f[2] + f[3])^2 + 0.25 * (f[1] - 4.0 * f[2] + 3.0 * f[3])^2
                β[2] = (13.0 / 12.0) * (f[2] - 2.0 * f[3] + f[4])^2 + 0.25 * (f[2] - f[4])^2
                β[3] = (13.0 / 12.0) * (f[3] - 2.0 * f[4] + f[5])^2 + 0.25 * (3.0 * f[3] - 4.0 * f[4] + f[5])^2

                # Global smoothness indicator
                τ5 = abs(β[1] - β[3])

                # Nonlinear weights
                α[1]     = d[1] * (1.0 + (τ5 / (β[1] + ϵ))^2 + λ * (β[1] + ϵ) / (τ5 + ϵ))
                α[2]     = d[2] * (1.0 + (τ5 / (β[2] + ϵ))^2 + λ * (β[2] + ϵ) / (τ5 + ϵ))
                α[3]     = d[3] * (1.0 + (τ5 / (β[3] + ϵ))^2 + λ * (β[3] + ϵ) / (τ5 + ϵ))
                w[1, :] .= (α ./ sum(α))

                # Add to advective term
                qTya[idx_x, idx_y - 2] += (w * q)[1]
            end
        end
    end

    # BC - Periodic
    qTxa[1, :]      .= qTxa[end-1, :]  # wrap last interface to start
    qTxa[end, :]    .= qTxa[2, :]      # wrap first real interface to end
    qTya[:, 1]      .= qTya[:, end-1]  # wrap last interface to start
    qTya[:, end]    .= qTya[:, 2]      # wrap first real interface to end

    # Compute advective change
    for idx_y in axes(∂T∂t_a, 2)
        for idx_x in axes(∂T∂t_a, 1)
            if idx_x > 3 && idx_x < size(∂T∂t_a, 1) - 2 && idx_y > 3 && idx_y < size(∂T∂t_a, 2) - 2
                ∂T∂t_a[idx_x, idx_y] = -(qTxa[idx_x - 1, idx_y - 3] - qTxa[idx_x - 2, idx_y - 3]) * _dx -(qTya[idx_x - 3, idx_y - 1] - qTya[idx_x - 3, idx_y - 2]) * _dy
            end
        end
    end

    # Return
    return nothing
end

# Time integration
function update_T_RK4_WENO!(
    T         :: AbstractArray,
    ∂T∂t_aRK1 :: Matrix{DatType},
    ∂T∂t_aRK2 :: Matrix{DatType},
    ∂T∂t_aRK3 :: Matrix{DatType},
    ∂T∂t_aRK4 :: Matrix{DatType},
    dt        :: DatType,
)
    # Update
    for idx_y in axes(T, 2)
        for idx_x in axes(T, 1)
            if idx_x > 3 && idx_x < size(T, 1) - 2 && idx_y > 3 && idx_y < size(T, 2) - 2
                T[idx_x, idx_y] += dt / 6.0 * (∂T∂t_aRK1[idx_x, idx_y] + 2.0 * ∂T∂t_aRK2[idx_x, idx_y] + 2.0 * ∂T∂t_aRK3[idx_x, idx_y] + ∂T∂t_aRK4[idx_x, idx_y])
            end
        end
    end

    # Return
    return nothing
end

# Define main function
@views function run_main()

    # Physics
    Lx   = 20.0
    Ly   = 14
    k    = 1.0
    Tmax = 1.0
    σh   = 0.5
    Vx   = 1.0
    Vy   = 1.0
    ω    = 1.0

    # Numerics
    ncx, ncy  = 101, 51
    nt   = 200
    dx   = Lx / ncx
    dy   = Ly / ncy
    _dx  = 1.0 / dx
    _dy  = 1.0 / dy
    CFL  = 0.4
    dt_d = 1.0 #min(dx^2, dy^2) / k  / 4.1
    nviz = 10

    # Initialize
    xc                       = collect(LinRange((-Lx+dx)/2.0, (Lx-dx)/2.0, ncx  ))
    yc                       = collect(LinRange((-Ly+dy)/2.0, (Ly-dy)/2.0, ncy  ))
    xv                       = collect(LinRange((-Lx   )/2.0, (Lx   )/2.0, ncx+1))
    yv                       = collect(LinRange((-Ly   )/2.0, (Ly   )/2.0, ncy+1))
    T_ini                    = zeros(DatType, ncx + 6, ncy + 6)
    for idx_y in eachindex(yc)
        for idx_x in eachindex(xc)
            T_ini[idx_x+3, idx_y+3] = Tmax * exp(-((xc[idx_x] + Lx/6.0)^2 + (yc[idx_y] + Ly/6.0)^2) / σh)
        end
    end
    T_WENOZp_RK = deepcopy(T_ini)
    qTxa_WENO   = zeros(DatType, ncx + 1, ncy    )
    qTya_WENO   = zeros(DatType, ncx    , ncy + 1)
    vxT         = zeros(DatType, ncx + 7, ncy)
    vyT         = zeros(DatType, ncx    , ncy + 7)
    ∂T∂t_aWRK1  = zeros(DatType, ncx + 6, ncy + 6)
    ∂T∂t_aWRK2  = zeros(DatType, ncx + 6, ncy + 6)
    ∂T∂t_aWRK3  = zeros(DatType, ncx + 6, ncy + 6)
    ∂T∂t_aWRK4  = zeros(DatType, ncx + 6, ncy + 6)
    vx          = zeros(DatType, ncx + 1, ncy    )
    vy          = zeros(DatType, ncx    , ncy + 1)
    for idx_y in axes(vx, 2)
        for idx_x in axes(vx, 1)
            vx[idx_x, idx_y] = -ω * yv[idx_y] * Vx
        end
    end
    for idx_y in axes(vy, 2)
        for idx_x in axes(vy, 1)
            vy[idx_x, idx_y] =  ω * xv[idx_x] * Vy
        end
    end
    α_WENO                   = [0.0, 0.0, 0.0]
    β_WENO                   = [0.0, 0.0, 0.0]
    q_WENO                   = zeros(DatType, 3, 1)
    f_WENO                   = [0.0, 0.0, 0.0, 0.0, 0.0]
    w_WENO                   = zeros(DatType, 1, 3)
    d_WENO                   = [0.1, 0.6, 0.3]
    λ_WENO                   = 0.0
    dt_a = min(dx,   dy  ) / max(abs(maximum(vx)), abs(maximum(vy))) / 4.1
    dt   = CFL * min(dt_d, dt_a)

    # Visualize initial configuration
    f = Figure()
    ax1 = Axis(f[1,1], xlabel = L"$$x []", ylabel = L"$$y []", aspect = DataAspect())
    ax2 = Axis(f[1,2], xlabel = L"$$x []", ylabel = L"$$y []", aspect = DataAspect())
    ax3 = Axis(f[2,1], xlabel = L"$$x []", ylabel = L"$$y []", aspect = DataAspect())
    hm1 = contourf!(ax1, xc, yc, T_WENOZp_RK[4:end-3, 4:end-3], colormap = :bilbao)
    hm2 = heatmap!(ax2, xv, yc, vx, colormap = :roma)
    hm3 = heatmap!(ax3, xc, yv, vy, colormap = :roma)
    cb  = Colorbar(f[2,2][1,1], hm1, label = L"$$T []")
    cb  = Colorbar(f[2,2][1,2], hm2, label = L"$$vx []")
    cb  = Colorbar(f[2,2][1,3], hm3, label = L"$$vy []")
    display(f)

    # Time loop
    for idx_t in 1:nt

        # WENO-Z+ - RK4 update
        ∂T∂t_aWRK1 .= 0.0
        ∂T∂t_aWRK2 .= 0.0
        ∂T∂t_aWRK3 .= 0.0
        ∂T∂t_aWRK4 .= 0.0
        advect_WENOZ_plus!(T_WENOZp_RK,                        vxT, vyT, qTxa_WENO, qTya_WENO, ∂T∂t_aWRK1, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, vy, _dx, _dy, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+ 0.5.*dt.*∂T∂t_aWRK1, vxT, vyT, qTxa_WENO, qTya_WENO, ∂T∂t_aWRK2, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, vy, _dx, _dy, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+ 0.5.*dt.*∂T∂t_aWRK2, vxT, vyT, qTxa_WENO, qTya_WENO, ∂T∂t_aWRK3, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, vy, _dx, _dy, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+      dt.*∂T∂t_aWRK3, vxT, vyT, qTxa_WENO, qTya_WENO, ∂T∂t_aWRK4, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, vy, _dx, _dy, λ_WENO)
        update_T_RK4_WENO!(T_WENOZp_RK, ∂T∂t_aWRK1, ∂T∂t_aWRK2, ∂T∂t_aWRK3, ∂T∂t_aWRK4, dt)

        # BC handling
        # T_WENOZp_RK[1:3, 1]       .= T_WENOZp_RK[end-5:end-3, 1]
        # T_WENOZp_RK[end-2:end, 1] .= T_WENOZp_RK[4:6, 1]

        # Update visualization
        if idx_t % nviz == 0
            # hm1[3][] .= T_WENOZp_RK[4:end-3, 4:end-3]
            hm1 = heatmap!(ax1, xc, yc, T_WENOZp_RK[4:end-3, 4:end-3], colormap = :bilbao)
            display(f)
        end

    end

    # Return
    return ln

end

# -----------------------------------------------------------------
# Run main function
ln = run_main();