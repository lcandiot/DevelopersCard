# This script showcases standard 1st order upwind advection and compare Runge-Kutta methods with Euler time integration
using CairoMakie, MathTeXEngine

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

# Advection functions
function advect_upwind!(
    T      :: AbstractArray,
    ∂T∂t_a :: AbstractArray,
    vx     :: DatType,
    _dx    :: DatType
)
    # Reset
    ∂T∂t_a .= 0.0

    # Advect
    for idx_x in eachindex(T)
        if idx_x > 1 && idx_x < size(T, 1)
            ∂T∂t_a[idx_x, 1] += -max(0.0, vx) * (T[idx_x,      1] - T[idx_x - 1, 1]) * _dx
            ∂T∂t_a[idx_x, 1] += -min(vx, 0.0) * (T[idx_x + 1 , 1] - T[idx_x,     1]) * _dx
        end
    end

    # Return
    return nothing
end

function advect_upwind_conservative!(
    T      :: AbstractArray,
    qTa    :: AbstractArray,
    ∂T∂t_a :: AbstractArray,
    vx     :: DatType,
    _dx    :: DatType
)
    # Reset
    ∂T∂t_a .= 0.0
    qTa    .= 0.0

    # Compute flux
    for idx_x in eachindex(qTa)
        if idx_x > 1 && idx_x < size(qTa, 1)
            qTa[idx_x] += max(0.0, vx) * T[idx_x - 1, 1]
            qTa[idx_x] += min(vx, 0.0) * T[idx_x    , 1]
        end
    end

    # Periodic BC handling
    vx > 0.0 ? qTa[end] = vx * T[end] : vx < 0.0 ? qTa[end] = vx * T[1  ] : qTa[end] = 0.0
    vx > 0.0 ? qTa[1  ] = vx * T[1  ] : vx < 0.0 ? qTa[1  ] = vx * T[end] : qTa[1  ] = 0.0

    # Advect
    for idx_x in eachindex(T)
        ∂T∂t_a[idx_x, 1] -= (qTa[idx_x + 1 , 1] - qTa[idx_x, 1]) * _dx
    end

    # Return
    return nothing
end

function advect_WENOZ!(
    T      :: AbstractArray,
    qTa    :: Matrix{DatType},
    ∂T∂t_a :: Matrix{DatType},
    α      :: Vector{DatType},
    β      :: Vector{DatType},
    d      :: Vector{DatType},
    q      :: Matrix{DatType},
    f      :: Vector{DatType},
    w      :: Matrix{DatType},
    vx     :: DatType,
    _dx    :: DatType
)

    # Reset
    ∂T∂t_a .= 0.0
    qTa    .= 0.0

    # Reconstruct flux in positive direction
    for idx_x in eachindex(T)
        if idx_x > 3 && idx_x < size(T, 1) - 2
            # Left-biased stencil - positive flow direction
            f[1], f[2], f[3], f[4], f[5] = T[idx_x - 3], T[idx_x - 2], T[idx_x - 1], T[idx_x], T[idx_x + 1]

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
            qTa[idx_x - 3] += max(0.0, vx) * (w * q)[1]

            # Right-biased stencil - negative flow direction
            f[1], f[2], f[3], f[4], f[5] = T[idx_x + 3], T[idx_x + 2], T[idx_x + 1], T[idx_x], T[idx_x - 1]
            
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
            qTa[idx_x - 2] += min(vx, 0.0) * (w * q)[1]
        end
    end

    # BC - Periodic
    qTa[1, 1]       = qTa[end-1, 1]  # wrap last interface to start
    qTa[end, 1]     = qTa[2, 1]      # wrap first real interface to end
    
    for idx_x in axes(∂T∂t_a, 1)
        if idx_x > 3 && idx_x < size(∂T∂t_a, 1) - 2
            ∂T∂t_a[idx_x] = -(qTa[idx_x - 2, 1] - qTa[idx_x - 3, 1]) * _dx
        end
    end

    # Return
    return nothing
end

function advect_WENOZ_plus!(
    T      :: AbstractArray,
    qTa    :: Matrix{DatType},
    ∂T∂t_a :: Matrix{DatType},
    α      :: Vector{DatType},
    β      :: Vector{DatType},
    d      :: Vector{DatType},
    q      :: Matrix{DatType},
    f      :: Vector{DatType},
    w      :: Matrix{DatType},
    vx     :: DatType,
    _dx    :: DatType,
    λ      :: DatType
)
    # Reset
    ∂T∂t_a .= 0.0
    qTa    .= 0.0

    # Reconstruct flux in positive direction
    for idx_x in eachindex(T)
        if idx_x > 3 && idx_x < size(T, 1) - 2
            # Left-biased stencil - positive flow direction
            f[1], f[2], f[3], f[4], f[5] = T[idx_x - 3], T[idx_x - 2], T[idx_x - 1], T[idx_x], T[idx_x + 1]

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
            qTa[idx_x - 3] += max(0.0, vx) * (w * q)[1]

            # Right-biased stencil - negative flow direction
            f[1], f[2], f[3], f[4], f[5] = T[idx_x + 3], T[idx_x + 2], T[idx_x + 1], T[idx_x], T[idx_x - 1]
            
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
            qTa[idx_x - 2] += min(vx, 0.0) * (w * q)[1]
        end
    end

    # BC - Periodic
    qTa[1, 1]       = qTa[end-1, 1]  # wrap last interface to start
    qTa[end, 1]     = qTa[2, 1]      # wrap first real interface to end
    
    for idx_x in axes(∂T∂t_a, 1)
        if idx_x > 3 && idx_x < size(∂T∂t_a, 1) - 2
            ∂T∂t_a[idx_x] = -(qTa[idx_x - 2, 1] - qTa[idx_x - 3, 1]) * _dx
        end
    end

    # Return
    return nothing
end

# Time integration
function update_T_Euler!(
    T      :: AbstractArray,
    ∂T∂t_a :: AbstractArray,
    dt     :: DatType
)
    # Update
    for idx_x in eachindex(T)
        if idx_x > 1 && idx_x < size(T, 1)
            T[idx_x, 1] += dt * ∂T∂t_a[idx_x, 1]
        end
    end

    # Return
    return nothing
end

function update_T_RK4!(
    T         :: AbstractArray,
    ∂T∂t_aRK1 :: Matrix{DatType},
    ∂T∂t_aRK2 :: Matrix{DatType},
    ∂T∂t_aRK3 :: Matrix{DatType},
    ∂T∂t_aRK4 :: Matrix{DatType},
    dt        :: DatType,
)
    # Update
    for idx_x in axes(T, 1)
        if idx_x > 1 && idx_x < size(T, 1)
            T[idx_x, 1] += dt / 6.0 * (∂T∂t_aRK1[idx_x, 1] + 2.0 * ∂T∂t_aRK2[idx_x, 1] + 2.0 * ∂T∂t_aRK3[idx_x, 1] + ∂T∂t_aRK4[idx_x, 1])
        end
    end

    # Return
    return nothing
end

function update_T_RK4_WENO!(
    T         :: AbstractArray,
    ∂T∂t_aRK1 :: Matrix{DatType},
    ∂T∂t_aRK2 :: Matrix{DatType},
    ∂T∂t_aRK3 :: Matrix{DatType},
    ∂T∂t_aRK4 :: Matrix{DatType},
    dt        :: DatType,
)
    # Update
    for idx_x in axes(T, 1)
        if idx_x > 3 && idx_x < size(T, 1) - 2
            T[idx_x, 1] += dt / 6.0 * (∂T∂t_aRK1[idx_x, 1] + 2.0 * ∂T∂t_aRK2[idx_x, 1] + 2.0 * ∂T∂t_aRK3[idx_x, 1] + ∂T∂t_aRK4[idx_x, 1])
        end
    end

    # Return
    return nothing
end

# Define main function
@views function run_main()

    # Physics
    Lx   = 20.0
    k    = 1.0
    Tmax = 1.0
    σh   = 0.5
    vx   = -1.0

    # Numerics
    ncx  = 101
    nt   = 10_000
    dx   = Lx / ncx
    _dx  = 1.0 / dx
    CFL  = 0.4
    dt_d = dx^2 / k  / 2.1
    dt_a = dx   / abs(vx) / 2.1
    dt   = CFL * min(dt_d, dt_a)
    nviz = 10

    # Initialize
    xc                       = collect(LinRange((-Lx+dx)/2.0, (Lx-dx)/2.0, ncx))
    T_ini                    = zeros(DatType, ncx + 2, 1)
    T_ini[2:end-1, :]       .= @. Tmax * exp(-((xc + Lx/3.0)^2) / σh)
    T_WENOZ_EU               = zeros(DatType, ncx + 6, 1)
    T_WENOZ_EU[4:end-3, :]  .= @. Tmax * exp(-((xc + Lx/3.0)^2) / σh)
    T_WENOZ_RK               = zeros(DatType, ncx + 6, 1)
    T_WENOZ_RK[4:end-3, :]  .= @. Tmax * exp(-((xc + Lx/3.0)^2) / σh)
    T_WENOZp_RK              = zeros(DatType, ncx + 6, 1)
    T_WENOZp_RK[4:end-3, :] .= @. Tmax * exp(-((xc + Lx/3.0)^2) / σh)
    T_EU                     = deepcopy(T_ini)
    T_RK                     = deepcopy(T_ini)
    qTd                      = zeros(DatType, ncx + 1, 1)
    qTa_WENO                 = zeros(DatType, ncx + 1, 1)
    ∂T∂t_a                   = zeros(DatType, ncx + 2, 1)
    ∂T∂t_aRK1                = zeros(DatType, ncx + 2, 1)
    ∂T∂t_aRK2                = zeros(DatType, ncx + 2, 1)
    ∂T∂t_aRK3                = zeros(DatType, ncx + 2, 1)
    ∂T∂t_aRK4                = zeros(DatType, ncx + 2, 1)
    ∂T∂t_aWRK1               = zeros(DatType, ncx + 6, 1)
    ∂T∂t_aWRK2               = zeros(DatType, ncx + 6, 1)
    ∂T∂t_aWRK3               = zeros(DatType, ncx + 6, 1)
    ∂T∂t_aWRK4               = zeros(DatType, ncx + 6, 1)
    α_WENO                   = [0.0, 0.0, 0.0]
    β_WENO                   = [0.0, 0.0, 0.0]
    q_WENO                   = zeros(DatType, 3, 1)
    f_WENO                   = [0.0, 0.0, 0.0, 0.0, 0.0]
    w_WENO                   = zeros(DatType, 1, 3)
    d_WENO                   = [0.1, 0.6, 0.3]
    λ_WENO                   = dx^(2/3)

    # Visualize initial configuration
    f = Figure()
    ax = Axis(f[1,1], xlabel = L"$$x []", ylabel = L"$$T []")
    ln0 = lines!(ax, xc, T_ini[2:end-1], color = :black, label = L"$$Initial")
    ln1 = lines!(ax, xc, T_EU[2:end-1], label = L"$$Euler")
    sc1 = scatter!(ax, xc, T_RK[2:end-1], label = L"$$RK4")
    ln2 = lines!(ax, xc, T_WENOZ_RK[4:end-3], color = :orange, label = L"$$WENO-Z")
    sc2 = scatter!(ax, xc, T_WENOZp_RK[4:end-3], color = :magenta, label = L"$$WENO-Z+")
    ln3 = lines!(ax, xc, T_WENOZp_RK[4:end-3], color = :magenta, label = L"$$WENO-Z+")
    axislegend(ax, merge = true)
    display(f)

    # Time loop
    for idx_t in 1:nt
        # Euler
        advect_upwind!(T_EU, ∂T∂t_a, vx, _dx)
        update_T_Euler!(T_EU, ∂T∂t_a, dt)

        # BC - Periodic
        T_EU[[1, end], 1] .= T_EU[[end - 1, 2], 1]

        # Runge Kutta 4
        ∂T∂t_aRK1 .= 0.0
        ∂T∂t_aRK2 .= 0.0
        ∂T∂t_aRK3 .= 0.0
        ∂T∂t_aRK4 .= 0.0
        advect_upwind!(T_RK, ∂T∂t_aRK1, vx, _dx)
        advect_upwind!(T_RK .+ 0.5.*dt.*∂T∂t_aRK1, ∂T∂t_aRK2, vx, _dx)
        advect_upwind!(T_RK .+ 0.5.*dt.*∂T∂t_aRK2, ∂T∂t_aRK3, vx, _dx)
        advect_upwind!(T_RK .+      dt.*∂T∂t_aRK3, ∂T∂t_aRK4, vx, _dx)
        update_T_RK4!(T_RK, ∂T∂t_aRK1, ∂T∂t_aRK2, ∂T∂t_aRK3, ∂T∂t_aRK4,dt)

        # BC - Periodic
        T_RK[[1, end], 1] .= T_RK[[end - 1, 2], 1]

        # WENO-Z - RK4 update
        ∂T∂t_aWRK1 .= 0.0
        ∂T∂t_aWRK2 .= 0.0
        ∂T∂t_aWRK3 .= 0.0
        ∂T∂t_aWRK4 .= 0.0
        advect_WENOZ!(T_WENOZ_RK,                         qTa_WENO, ∂T∂t_aWRK1, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx)
        advect_WENOZ!(T_WENOZ_RK .+ 0.5.*dt.*∂T∂t_aWRK1,  qTa_WENO, ∂T∂t_aWRK2, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx)
        advect_WENOZ!(T_WENOZ_RK .+ 0.5.*dt.*∂T∂t_aWRK2,  qTa_WENO, ∂T∂t_aWRK3, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx)
        advect_WENOZ!(T_WENOZ_RK .+      dt.*∂T∂t_aWRK3,  qTa_WENO, ∂T∂t_aWRK4, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx)
        update_T_RK4_WENO!(T_WENOZ_RK, ∂T∂t_aWRK1, ∂T∂t_aWRK2, ∂T∂t_aWRK3, ∂T∂t_aWRK4, dt)

        # BC handling
        # T_WENOZ_RK[1:3, 1]       .= T_WENOZ_RK[end-5:end-3, 1]
        # T_WENOZ_RK[end-2:end, 1] .= T_WENOZ_RK[4:6, 1]

        # WENO-Z+ - RK4 update
        ∂T∂t_aWRK1 .= 0.0
        ∂T∂t_aWRK2 .= 0.0
        ∂T∂t_aWRK3 .= 0.0
        ∂T∂t_aWRK4 .= 0.0
        advect_WENOZ_plus!(T_WENOZp_RK,                         qTa_WENO, ∂T∂t_aWRK1, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+ 0.5.*dt.*∂T∂t_aWRK1,  qTa_WENO, ∂T∂t_aWRK2, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+ 0.5.*dt.*∂T∂t_aWRK2,  qTa_WENO, ∂T∂t_aWRK3, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx, λ_WENO)
        advect_WENOZ_plus!(T_WENOZp_RK .+      dt.*∂T∂t_aWRK3,  qTa_WENO, ∂T∂t_aWRK4, α_WENO, β_WENO, d_WENO, q_WENO, f_WENO, w_WENO, vx, _dx, λ_WENO)
        update_T_RK4_WENO!(T_WENOZp_RK, ∂T∂t_aWRK1, ∂T∂t_aWRK2, ∂T∂t_aWRK3, ∂T∂t_aWRK4, dt)

        # BC handling
        # T_WENOZp_RK[1:3, 1]       .= T_WENOZp_RK[end-5:end-3, 1]
        # T_WENOZp_RK[end-2:end, 1] .= T_WENOZp_RK[4:6, 1]

        # Update visualization
        if idx_t % nviz == 0
            ln1[1][] .= [[pt[1], T_EU[i+1]   ] for (i, pt) in enumerate(ln1[1][])]
            sc1[1][] .= [[pt[1], T_RK[i+1]   ] for (i, pt) in enumerate(sc1[1][])]
            ln2[1][] .= [[pt[1], T_WENOZ_RK[i+3]] for (i, pt) in enumerate(ln2[1][])]
            sc2[1][] .= [[pt[1], T_WENOZp_RK[i+3]] for (i, pt) in enumerate(sc2[1][])]
            ln3[1][] .= [[pt[1], T_WENOZp_RK[i+3]] for (i, pt) in enumerate(ln3[1][])]
            display(f)
        end

    end

    # Return
    return ln

end

# -----------------------------------------------------------------
# Run main function
ln = run_main();