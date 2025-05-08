# This script showcases standard 1st order upwind advection and compare Runge-Kutta methods with Euler time integration
using CairoMakie, MathTeXEngine

# Set data type
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
        if idx_x < size(T, 1)
            ∂T∂t_a[idx_x + 1, 1] += -max(0.0, vx) * (T[idx_x + 1 , 1] - T[idx_x, 1]) * _dx
            ∂T∂t_a[idx_x    , 1] += -min(vx, 0.0) * (T[idx_x + 1 , 1] - T[idx_x, 1]) * _dx
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
        T[idx_x, 1] += dt * ∂T∂t_a[idx_x, 1]
    end

    # Return
    return nothing
end

@views function update_T_RK4!(
    T         :: AbstractArray,
    ∂T∂t_aRK1 :: Matrix{DatType},
    ∂T∂t_aRK2 :: Matrix{DatType},
    ∂T∂t_aRK3 :: Matrix{DatType},
    ∂T∂t_aRK4 :: Matrix{DatType},
    vx        :: DatType,
    dt        :: DatType,
    _dx       :: DatType
)
    # RK terms
    ∂T∂t_aRK1 .= 0.0
    ∂T∂t_aRK2 .= 0.0
    ∂T∂t_aRK3 .= 0.0
    advect_upwind!(T, ∂T∂t_aRK1, vx, _dx)
    advect_upwind!(T .+ 0.5.*dt.*∂T∂t_aRK1, ∂T∂t_aRK2, vx, _dx)
    advect_upwind!(T .+ 0.5.*dt.*∂T∂t_aRK2, ∂T∂t_aRK3, vx, _dx)
    advect_upwind!(T .+      dt.*∂T∂t_aRK3, ∂T∂t_aRK4, vx, _dx)

    # Update
    for idx_x in axes(T, 1)
        T[idx_x, 1] += dt / 6.0 * (∂T∂t_aRK1[idx_x, 1] + 2.0 * ∂T∂t_aRK2[idx_x, 1] + 2.0 * ∂T∂t_aRK3[idx_x, 1] + ∂T∂t_aRK4[idx_x, 1])
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
    vx   = 1.0

    # Numerics
    ncx  = 51
    nt   = 100
    dx   = Lx / ncx
    _dx  = 1.0 / dx
    CFL  = 1.0
    dt_d = dx^2 / k  / 2.1
    dt_a = dx   / abs(vx) / 2.1
    dt   = CFL * min(dt_d, dt_a)
    nviz = 25

    # Initialize
    xc                 = collect(LinRange((-Lx+dx)/2.0, (Lx-dx)/2.0, ncx))
    T_ini              = zeros(DatType, ncx + 2, 1)
    T_ini[2:end-1, :] .= @. Tmax * exp(-((xc + Lx/3.0)^2) / σh)
    T_EU               = deepcopy(T_ini)
    T_RK               = deepcopy(T_ini)
    qT                 = zeros(DatType, ncx+1, 1)
    ∂T∂t_a             = zeros(DatType, ncx  , 1)
    ∂T∂t_aRK1          = zeros(DatType, ncx  , 1)
    ∂T∂t_aRK2          = zeros(DatType, ncx  , 1)
    ∂T∂t_aRK3          = zeros(DatType, ncx  , 1)
    ∂T∂t_aRK4          = zeros(DatType, ncx  , 1)

    # Visualize initial configuration
    f = Figure()
    ax = Axis(f[1,1], xlabel = L"$$x []", ylabel = L"$$T []")
    ln0 = lines!(ax, xc, T_ini[2:end-1], color = :black, label = L"$$Initial")
    ln1 = lines!(ax, xc, T_EU[2:end-1], label = L"$$Euler")
    sc1 = scatter!(ax, xc, T_RK[2:end-1], label = L"$$RK4")
    axislegend(ax)
    display(f)

    # Time loop
    for idx_t in 1:nt
        # Euler
        advect_upwind!(T_EU[2:end-1, 1], ∂T∂t_a, vx, _dx)
        update_T_Euler!(T_EU[2:end-1, 1], ∂T∂t_a, dt)

        # Runge Kutta 4
        update_T_RK4!(T_RK[2:end-1, 1], ∂T∂t_aRK1, ∂T∂t_aRK2, ∂T∂t_aRK3, ∂T∂t_aRK4, vx, dt, _dx)

        # Update visualization
        if idx_t % nviz == 0
            ln1[1][] .= [[pt[1], T_EU[i+1]] for (i, pt) in enumerate(ln1[1][])]
            sc1[1][] .= [[pt[1], T_RK[i+1]] for (i, pt) in enumerate(sc1[1][])]
            display(f)
        end

    end

    # Return
    return ln

end

# -----------------------------------------------------------------
# Run main function
ln = run_main();