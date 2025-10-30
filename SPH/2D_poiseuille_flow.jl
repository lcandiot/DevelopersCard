using CairoMakie, Statistics, Random, GeoUtils, MathTeXEngine
using ColorSchemes, Printf

const DatType = Float64
const NEIGH_OFFSET = (
    ( 0,  0), # cell center
    (-1,  0), # left neighbour
    ( 1,  0), # right neighbour
    ( 0,  1), # top neighbour
    ( 0, -1), # bottom neighbour
    (-1,  1), # top left neighbour
    (-1, -1), # bottom left neighbour
    ( 1,  1), # top right neighbour
    ( 1, -1), # bottom right neighbour
)

# ====== Types and structs ====== #
abstract type SPH end
abstract type AbstractFlow     <: SPH end
abstract type AbstractKernels  <: SPH end
abstract type AbstractRheology <: SPH end

struct PoiseuilleFlow     <: AbstractFlow     end
struct CubicSpline        <: AbstractKernels  end
struct CubicSplineWithAcc <: AbstractKernels  end
struct MorrisFD           <: AbstractRheology end

struct Particles <: SPH
    x           :: Vector{DatType}
    y           :: Vector{DatType}
    vx          :: Vector{DatType}
    vy          :: Vector{DatType}
    ρ           :: Vector{DatType}
    p           :: Vector{DatType}
    m           :: Vector{DatType}
    ax          :: Vector{DatType}
    ay          :: Vector{DatType}
    wall_flag   :: Vector{Int}
    wall_position :: Vector{Int}
    idx_wall    :: UnitRange{Int}
    idx_fluid   :: UnitRange{Int}
end

struct PoiseuilleParams <: SPH
    # Geometry and fluid
    Lx      :: DatType
    Ly      :: DatType
    ρ0      :: DatType
    ν       :: DatType
    a_drive :: DatType
    # Discretization
    Δp      :: DatType
    h       :: DatType
    η_num   :: DatType
    # EOS
    γ       :: DatType
    c0      :: DatType
    B       :: DatType
    # Lattice and counts
    nxc     :: Int
    nyc     :: Int
    dx      :: DatType
    dy      :: DatType
    wall_rows_per_side :: Int
    np_fluid :: Int
    np_wall  :: Int
    np       :: Int
end

# ====== Constructors ====== #
function ConstructParams(::PoiseuilleFlow;
    Δp::DatType, Ly::DatType, Lx::DatType, ρ0::DatType, ν::DatType,
    a_drive::DatType, γ::DatType, c0::DatType)

    B   = (c0^2) * ρ0 / γ
    h   = 1.3 * Δp
    η_num = 0.01 * h^2

    # --- target hex aspect (ax, ay) and SNAP to Ly, then Lx ---
    ay_guess = DatType(√3/2) * Δp
    nyc      = max(1, round(Int, Ly / ay_guess))       # pick integer rows
    ay       = Ly / nyc                                # exact: Ly = nyc * ay
    ax       = 2ay / √3                                # keep perfect hex ratio

    nxc      = max(1, floor(Int, Lx / ax))             # columns that fit
    ax       = Lx / nxc                                # exact: Lx = nxc * ax   (periodic x looks cleaner)
    # NOTE: ax, ay are the *effective* lattice spacings now

    # store them in dx, dy (so init can use params.dx/params.dy)
    dx, dy = ax, ay

    wall_rows_per_side = 5
    npf = nxc * nyc
    npw = 2 * wall_rows_per_side * nxc
    np  = npf + npw

    return PoiseuilleParams(Lx, Ly, ρ0, ν, a_drive, Δp, h, η_num,
                            γ, c0, B, nxc, nyc, dx, dy,
                            wall_rows_per_side, npf, npw, np)
end

function InitializeParticles(::PoiseuilleFlow, params::PoiseuilleParams)
    # Geometry of the hex grid
    ax = params.Δp
    ay = DatType(sqrt(3)/2) * ax

    # Allocate arrays
    x  = Vector{DatType}(undef, params.np)
    y  = Vector{DatType}(undef, params.np)
    vx = fill(DatType(0), params.np)
    vy = fill(DatType(0), params.np)
    axv= fill(DatType(0), params.np)
    ayv= fill(DatType(0), params.np)
    ρ  = fill(params.ρ0, params.np)
    p  = params.B .* ((ρ ./ params.ρ0).^params.γ .- 1.0)
    wall_flag = fill(1, params.np)
    wall_position = fill(0, params.np)

    # Particle mass from area per particle in hex packing
    Ap = ax * ay
    m  = fill(params.ρ0 * Ap, params.np)

    # Indices
    idx_fluid = 1:params.np_fluid
    idx_wall  = params.np_fluid+1 : params.np

    # ---- Fluid positions (hex) ----
    c = 0
    for j in 0:params.nyc-1
        yj = (j + 0.5) * ay
        xshift = (isodd(j) ? ax/2 : 0.0)
        for i in 0:params.nxc-1
            c += 1
            xi = (i + 0.5) * ax + xshift
            # periodic wrap in x
            xi = reinsert_periodic_part(xi, params.Lx)
            x[c] = xi
            y[c] = yj
        end
    end
    @assert c == params.np_fluid

    # ---- Wall rows: continue the parity above/below ----
    # Bottom walls: j = -1, -2, ...
    nw = params.wall_rows_per_side
    for k in 1:nw
        jghost = -k                   # virtual row index just below j=0
        yj = (jghost + 0.5) * ay      # -> -0.5*ay, -1.5*ay, ...
        xshift = (isodd(jghost) ? ax/2 : 0.0)

        for i in 0:params.nxc-1
            c += 1
            xi = (i + 0.5) * ax + xshift
            xi = reinsert_periodic_part(xi, params.Lx)
            x[c] = xi
            y[c] = yj
            wall_flag[c] = 0
            wall_position[c] = -1                 # Flagging the boundary poisition: -1 = Bottom, 1 = Top
            vx[c] = 0.0; vy[c] = 0.0
            ρ[c] = params.ρ0; p[c] = 0.0
        end
    end

    # Top walls: continue after the last fluid row (j = nyc, nyc+1, ...)
    for k in 1:nw
        jghost = params.nyc - 1 + k   # just above the last fluid j = nyc-1
        yj = (jghost + 0.5) * ay      # -> Ly + 0.5*ay, Ly + 1.5*ay, ...
        xshift = (isodd(jghost) ? ax/2 : 0.0)

        for i in 0:params.nxc-1
            c += 1
            xi = (i + 0.5) * ax + xshift
            xi = reinsert_periodic_part(xi, params.Lx)
            x[c] = xi
            y[c] = yj
            wall_flag[c] = 0
            wall_position[c] = 1                # Flagging the boundary poisition: -1 = Bottom, 1 = Top
            vx[c] = 0.0; vy[c] = 0.0
            ρ[c] = params.ρ0; p[c] = 0.0
        end
    end

    @assert c == params.np  # sanity

    return Particles(x, y, vx, vy, ρ, p, m, axv, ayv, wall_flag, wall_position, idx_wall, idx_fluid)
end

# ====== Functions ====== #

function randomize_initial_particle_position!(particles::Particles, params::PoiseuilleParams; frac=0.10, rng=Random.default_rng())
    # frac is the jitter amplitude as a fraction of dx,dy (e.g. 0.10 = ±10%)
    jx = frac * params.dx
    jy = frac * params.dy
    @inbounds for ip in particles.idx_fluid
        particles.x[ip] += (2rand(rng) - 1) * jx
        particles.y[ip] += (2rand(rng) - 1) * jy
        # keep within domain for fluid (periodic in x; clamp in y)
        particles.x[ip] = reinsert_periodic_part(particles.x[ip], params.Lx)
        particles.y[ip] = clamp(particles.y[ip], 0.0, params.Ly)  # avoid penetrating walls at t=0
    end
    return nothing
end
@inline function compute_min_part_spacing(dx::Float64, L::Float64)      # This avoids that reinserted particles take the dp across the boundary and not the entire domain
    dx - L * round(dx / L)
end

@inline function reinsert_periodic_part(x::Float64, L::Float64)         # This reinserts particles that exit the domain
    x - floor(x / L) * L
end

@inline function W_kernel(
           ::CubicSpline,
    params :: PoiseuilleParams,
    Δx     :: DatType,
    Δy     :: DatType,
    α      :: DatType
)
    r2 = Δx*Δx + Δy*Δy
    r  = sqrt(r2)
    q  = r / params.h
    if q >= 2.0
        return 0.0
    elseif q < 1.0
        return α * (1 - 1.5*q*q + 0.75*q*q*q)
    else
        t = 2.0 - q
        return α * 0.25 * t*t*t
    end
end

@inline function ∇W_kernel(
           ::CubicSpline,
    params :: PoiseuilleParams,
    Δx     :: DatType,
    Δy     :: DatType,
    α      :: DatType
)
    r2 = Δx*Δx + Δy*Δy
    if r2 == 0.0
        return (0.0, 0.0)
    end
    r = sqrt(r2)
    q = r / params.h
    if q >= 2.0
        return (0.0, 0.0)
    end
    dWdq = (q < 1.0) ? (α * (-3.0*q + 2.25*q*q)) :
                       (α * (-0.75 * (2.0 - q) * (2.0 - q)))
    s = (dWdq / params.h) / r
    return (s*Δx, s*Δy)
end

@inline function ∇W_kernel(
           ::CubicSplineWithAcc,
    params :: PoiseuilleParams,
    Δx     :: DatType,
    Δy     :: DatType,
    α      :: DatType
)
    r2 = Δx*Δx + Δy*Δy
    if r2 == 0.0
        return (α, 0.0, 0.0, 0.0, 0.0)  # W(0)=αW, grad=0
    end
    r = sqrt(r2)
    q = r / params.h
    if q >= 2.0
        return (0.0, 0.0, 0.0, r2, r)
    end
    if q < 1.0
        Wv   = α * (1 - 1.5*q*q + 0.75*q*q*q)
        dWdq = α * (-3.0*q + 2.25*q*q)
    else
        t    = 2.0 - q
        Wv   = α * 0.25 * t*t*t
        dWdq = α * (-0.75 * t*t)
    end
    s = (dWdq / params.h) / r
    return (Wv, s*Δx, s*Δy, r2, r)
end

function create_CLL!(
    head      :: Vector{Int},
    next      :: Vector{Int},
    particles :: Particles,
    pb        :: NamedTuple,
    CLL       :: NamedTuple
)
    x0 = pb.xmin
    y0 = pb.ymin
    fill!(head, 0)
    fill!(next, 0)
    for ip in eachindex(particles.x)
        cx = 1 + mod(floor(Int, (particles.x[ip] - x0) / CLL.lcell), CLL.nx)
        cy = 1 +     floor(Int, (particles.y[ip] - y0) / CLL.lcell)
        cy = clamp(cy, 1, CLL.ny)
        cn = cx + (cy - 1) * CLL.nx
        next[ip] = head[cn]
        head[cn] = ip
    end

    # Sanity checks - Deactivate for production !
    # counts = zeros(Int, length(head))
    # for c in 1:length(head)
    #     p = head[c]
    #     while p != 0
    #         counts[c] += 1
    #         p = next[p]
    #     end
    # end
    # @assert sum(counts) == size(particles.x, 1)
end

@inline function reset_particle_arrays!(
    particles :: Particles,
    Arr       :: Vector{DatType},
    val       :: Union{DatType, Int}
)
    @inbounds for ip in particles.idx_fluid
        Arr[ip] = val
    end
    
    # Return
    return nothing
end

# ====== Physics functions ====== #
function compute_density!(
              ::PoiseuilleFlow,
    params    :: PoiseuilleParams,
    particles :: Particles,
    pb        :: NamedTuple,
    CLL       :: NamedTuple,
    head      :: Vector{Int},
    next      :: Vector{Int},
    α         :: DatType
)

    # Support
    supp2 = (2.0 * params.h)^2

    # Loop through particles
    @inbounds for ip in eachindex(particles.ρ)

        # Set numerator and denominator for Shepard Normalization
        num = 0.0
        den = 0.0

        # Get coords and initialize local density
        xp = particles.x[ip]; yp = particles.y[ip]

        # Determine cell of current particle
        cx = 1 + mod(floor(Int, (xp - pb.xmin) / CLL.lcell), CLL.nx)
        cy = 1 +     floor(Int, (yp - pb.ymin) / CLL.lcell)
        cy = clamp(cy, 1, CLL.ny)

        # Walk through neighbour cells
        for (off_x, off_y) in NEIGH_OFFSET
            # Determine neighbour cell (self is included here)
            cx2 = 1 + mod(cx - 1 + off_x, CLL.nx)   # stays in 1..nx by modulo (wraps)
            cy2 = clamp(cy + off_y, 1, CLL.ny)      # stays in 1..ny by clamping
            c   = cx2 + (cy2 - 1) * CLL.nx

            # Loop through neighbour particles
            j = head[c]
            while j != 0
                Δx = xp - particles.x[j]; Δy = yp - particles.y[j]
                Δx = compute_min_part_spacing(Δx, params.Lx)

                # Sum density in support zone
                r2 = Δx*Δx + Δy*Δy
                if r2 < supp2
                    Wij = W_kernel(CubicSpline(), params, Δx, Δy, α)
                    num += particles.m[j] * Wij
                    den += particles.m[j] / params.ρ0 * Wij
                end

                # Set next particle
                j = next[j]
            end
        end

        # Store calculated density
        particles.ρ[ip] = (den > 0.0) ? (num / den) : particles.ρ[ip]
    end
end

function density_to_pressure!(
    particles :: Particles,
    params    :: PoiseuilleParams
)
    @inbounds for ip in particles.idx_fluid
        particles.p[ip] = params.B * ((particles.ρ[ip]/params.ρ0)^params.γ - 1.0)
    end

    # Return
    return nothing
end

function compute_pressuregradient_acceleration!(
    params    ::PoiseuilleParams,
    particles :: Particles,
    pb        :: NamedTuple,
    CLL       :: NamedTuple,
    head      :: Vector{Int},
    next      :: Vector{Int},
    α         :: DatType
)
    # Support
    supp2 = (2.0 * params.h)^2

    # Loop over particles§
    @inbounds for ip in particles.idx_fluid
        xp = particles.x[ip]; yp = particles.y[ip]
        pp = particles.p[ip]; ρp = particles.ρ[ip]

        # Find lattice cell of current particle
        cx = 1 + mod(floor(Int, (xp - pb.xmin) / CLL.lcell), CLL.nx)
        cy = 1 +     floor(Int, (yp - pb.ymin) / CLL.lcell)
        cy = clamp(cy, 1, CLL.ny)

        # Loop through neighbour cells
        for (off_x, off_y) in NEIGH_OFFSET
            # Determine neighbour cell (self is included here)
            cx2 = 1 + mod(cx - 1 + off_x, CLL.nx)   # stays in 1..nx by modulo (wraps)
            cy2 = clamp(cy + off_y, 1, CLL.ny)      # stays in 1..ny by clamping
            c   = cx2 + (cy2 - 1) * CLL.nx

            # Loop through neighbouring particles
            j = head[c]
            while j != 0
                # 1: Exclude self -> Including itself would make the ∇Wij = 0 and only add computational cost without value
                if (ip == j); j = next[j]; continue end

                # 2: Exlude wall
                # if (particles.wall_flag[j] == 0); j = next[j]; continue end

                # 3: Now do the expensive stuff
                xj = particles.x[j]; yj = particles.y[j]
                Δx = xp - xj; Δy = yp - yj
                Δx = compute_min_part_spacing(Δx, params.Lx)
                r2 = Δx*Δx + Δy*Δy
                if r2 < supp2
                    pj = particles.p[j]
                    ρj = particles.ρ[j]
                    if particles.wall_flag[j] == 0
                        pj = pp
                        ρj = ρp
                    end
                    coeff    = -particles.m[j] * (pp / ρp^2 + pj/ρj^2)
                    dWx, dWy = ∇W_kernel(CubicSpline(), params, Δx, Δy, α)
                    particles.ax[ip] += coeff * dWx
                    particles.ay[ip] += coeff * dWy
                end
                j = next[j]
            end
        end
    end

    # Return
    return nothing
end

function compute_viscousterm_acceleration!(
              ::MorrisFD,
    params    :: PoiseuilleParams,
    particles :: Particles,
    pb        :: NamedTuple,
    CLL       :: NamedTuple,
    head      :: Vector{Int},
    next      :: Vector{Int},
    α         :: DatType
)
    # Support
    supp2 = (2.0 * params.h)^2
    η2    = params.η_num 

    # Loop over particles
    @inbounds for ip in particles.idx_fluid
        xp = particles.x[ip]; yp = particles.y[ip]
        vxp = particles.vx[ip]; vyp = particles.vy[ip]
        ρp  = particles.ρ[ip]

        # Find lattice cell of current particle
        cx = 1 + mod(floor(Int, (xp - pb.xmin) / CLL.lcell), CLL.nx)
        cy = 1 +     floor(Int, (yp - pb.ymin) / CLL.lcell)
        cy = clamp(cy, 1, CLL.ny)

        # Loop through neighbour cells
        for (off_x, off_y) in NEIGH_OFFSET
            # Determine neighbour cell (self is included here)
            cx2 = 1 + mod(cx - 1 + off_x, CLL.nx)   # stays in 1..nx by modulo (wraps)
            cy2 = clamp(cy + off_y, 1, CLL.ny)      # stays in 1..ny by clamping
            c   = cx2 + (cy2 - 1) * CLL.nx

            # Loop through neighbour particles
            j = head[c]
            while j != 0
                # 1: Exclude self -> Including itself would make the ∇Wij = 0 and only add computational cost without value
                if (ip == j); j = next[j]; continue end
                
                # 2: Now do the expensive stuff
                Δx = xp - particles.x[j]; Δy = yp - particles.y[j]
                Δx = compute_min_part_spacing(Δx, params.Lx)
                r2 = Δx*Δx + Δy*Δy
                if r2 >= supp2; j = next[j]; continue end
                
                # 3: kernel
                dWx, dWy = ∇W_kernel(CubicSpline(), params, Δx, Δy, α)
                
                # 4: Heavy stuff: If neighbour particle is fluid ...
                vxj = particles.vx[j]; vyj = particles.vy[j]
                xj = particles.x[j]; yj = particles.y[j]
                ρj = particles.ρ[j]

                if particles.wall_flag[j] == 0 # if neighbour particle is wall ...
                    # vx_pj, vy_pj = 0.0, 0.0
                    # if particles.wall_position[j] == -1
                    #     dp = abs(0.0 - yp)
                    #     dj = abs(0.0 - yj)
                    #     β = min(1.5, dj/dp)
                    #     vx_pj = -vxp
                    #     vy_pj = -vyp
                    # elseif particles.wall_position[j] == 1
                    #     dp = abs(params.Ly - yp)
                    #     dj = abs(params.Ly - yj)
                    #     β = min(1.5, dj/dp)
                    #     vx_pj = β * vxp
                    #     vy_pj = β * vyp
                    # end
                    vxj = -vxp
                    vyj = -vyp
                    ρj  = ρp
                end
                coeff_x = particles.m[j] * (params.ν * ρp + params.ν * ρj) / (r2 + η2) / (ρp * ρj)
                coeff_y = particles.m[j] * (params.ν * ρp + params.ν * ρj) / (r2 + η2) / (ρp * ρj)
                particles.ax[ip] += coeff_x * (Δx*dWx + Δy*dWy) * (vxp - vxj)
                particles.ay[ip] += coeff_y * (Δx*dWx + Δy*dWy) * (vyp - vyj)
                
                j = next[j]
            end
        end
    end

    # Return
    return nothing
end

function driving_force!(particles :: Particles, params :: PoiseuilleParams)
    @inbounds for ip in particles.idx_fluid
        particles.ax[ip] += params.a_drive
    end
end

function compute_stable_dt(
    ::PoiseuilleFlow,
    params::PoiseuilleParams,
    particles::Particles;
    CFL_cfl::DatType = 0.25,
    CFL_ν::DatType = 0.15,
    CFL_a::DatType = 0.25
)
    vmag = @views sqrt.(particles.vx .^ 2 .+ particles.vy .^ 2)
    amag = @views sqrt.(particles.ax .^ 2 .+ particles.ay .^ 2)
    vmax = maximum(vmag[particles.idx_fluid])
    amax = maximum(amag[particles.idx_fluid])

    dt_cfl = CFL_cfl * params.h / (params.c0 + vmax)
    dt_ν   = CFL_ν   * params.h^2 / params.ν
    dt_a   = CFL_a   * sqrt(params.h / max(amax, eps()))

    return min(dt_cfl, min(dt_ν, dt_a))
end

function integrate_symplectic!(
    particles::Particles,
    params::PoiseuilleParams,
    dt::DatType
)
    @inbounds for ip in particles.idx_fluid
        particles.vx[ip] += dt * particles.ax[ip]
        particles.vy[ip] += dt * particles.ay[ip]
        particles.x[ip]  += dt * particles.vx[ip]
        particles.y[ip]  += dt * particles.vy[ip]

        # Periodicity in x
        particles.x[ip] = reinsert_periodic_part(particles.x[ip], params.Lx)
    end
end

function sample_along_segment(particles::Particles, params::PoiseuilleParams, Q::AbstractVector;
                              p0::Tuple, p1::Tuple, n::Int=50, α::Float64=10.0/(7π*params.h^2))
    x0,y0 = p0; x1,y1 = p1
    xs = range(x0, x1; length=n) |> collect
    ys = range(y0, y1; length=n) |> collect
    Qs = similar(xs, Float64)

    supp2 = (2.0 * params.h)^2

    @inbounds for k in eachindex(xs)
        xr = xs[k]; yr = ys[k]
        num = 0.0; den = 0.0
        for j in eachindex(particles.x)
            dx = xr - particles.x[j]
            dy = yr - particles.y[j]
            dx = compute_min_part_spacing(dx, params.Lx)  # periodic x
            r2 = dx*dx + dy*dy
            if r2 < supp2
                Wij = W_kernel(CubicSpline(), params, dx, dy, α)
                wj  = particles.m[j] / particles.ρ[j]
                num += wj * Q[j] * Wij
                den += wj * Wij
            end
        end
        Qs[k] = (den > 0.0) ? num/den : 0.0
    end
    return xs, ys, Qs
end

function analytical_poiseuille(y, L, F, ν, t; N = 50)
    # F = body force per unit mass (acceleration)
    vx = @. (F/(2ν)) * y * (L - y)
    for n in 0:N-1
        vx .-= @. (4F*L^2)/(ν*π^3*(2n+1)^3) *
               sin((2n+1)*π*y/L) *
               exp(-((2n+1)^2*π^2*ν/L^2)*t)
    end
    return vx
end

# ====== MAIN ====== #
function main_2D_poiseuille_flow()

    # Model parameters
    Ly     = DatType(1.0e-3)                    # Height [m]
    Lx     = 1.0*Ly                    # Length [m]
    Δp     = DatType(2e-5)                      # Particle spacing [m]
    ρ0     = DatType(1000.0)                    # Density [kg/m3]
    ν      = DatType(1e-6)                      # Kinematic viscosity [m2/s]
    a_BC   = DatType(1e-2)                       # Driving acceleration at the boundary [m/s2]
    γ_Tait = DatType(7.0)                       # EOS for pressure calculation [-]
    c0     = DatType(30.0)                       # artificial speed of sound [m/s]
    nt     = 1_000_000                                # Number of time steps
    time   = 0.0                                # Phyiscal time
    nviz   = 10_000                                 # Plotting frequency
    ncheck = nviz                                 # Monitoring frequency

    # Initialize
    params    = ConstructParams(PoiseuilleFlow(), Δp = Δp, Ly = Ly, Lx = Lx, ρ0 = ρ0, ν = ν, a_drive = a_BC, γ = γ_Tait, c0 = c0)
    α         = 10.0 / (7.0 * π * params.h^2)      # Normalization factor for the weight kernel
    particles = InitializeParticles(PoiseuilleFlow(), params)
    # randomize_initial_particle_position!(particles, params; frac = 0.05)
    pb        = (                                  # Particle bounds used for the cell linked list -> move to particle struct later
        xmin = 0.0,
        xmax = params.Lx,
        ymin = minimum(particles.y),
        ymax = maximum(particles.y)
    )
    lcell = DatType(2.0 * params.h)
    CLL = (
        lcell = lcell,
        ny = ceil(Int, (pb.ymax - pb.ymin) / lcell),
        nx = ceil(Int, (pb.xmax - pb.xmin) / lcell),
    )
    head   = fill(0, CLL.nx*CLL.ny)
    next   = fill(0, params.np)
    dt     = 1e-6

    # Analytical solution
    _, ys, u_sph = sample_along_segment(particles, params, particles.vx; p0 = (params.Lx/2, 0), p1 = (params.Lx/2, params.Ly))
    u_ana = analytical_poiseuille(ys, params.Ly, params.a_drive, params.ν, 0.0001)
    u_max = a_BC * Ly^2 / 8.0 / ν
    
    # Initialize figure
    idxF = particles.idx_fluid[1:1:end]
    tname = @sprintf "Time: %.5f s" time
    f   = Figure(size = (750, 300), fontsize = 15)
    ax1 = Axis(f[1,1][1,1], limits = (0.0, params.Lx, 0.0, params.Ly), xlabel = L"$x$ [m]", ylabel =L"$y$ [m]", aspect = 1.0, title = tname)
    ax2 = Axis(f[1,2], xlabel = L"$y$ [m]", ylabel = L"$u_x$ [m.s$^{-1}$]", aspect = 1.0, limits = (0.0, Ly, 0.0, u_max))
    sc1 = scatter!(ax1, particles.x[idxF], particles.y[idxF], color = particles.vx[idxF], colormap=:batlow, colorrange = (0.0, u_max), markersize = 4)
    sc2 = scatter!(ax2, ys, u_sph, label = L"$$SPH")
    ln1 = lines!(ax2, ys, u_ana, color = :goldenrod, linewidth = 3, label = L"$$Analytical")
    cb  = Colorbar(f[1,1][1,2], sc1, tellheight = false, label = L"$u_x$ [m.s$^{-1}$]")
    axislegend(ax2, labelsize = 11, backgroundcolor = (:white, 0.0))
    display(f)

    # Interpolation settings - to get the Poiseuille parabola
    nsample   = 50
    u_sph_mat = Matrix{Float64}(undef, size(ys, 1), nsample)
    dx_spl = params.Lx / (nsample - 1)
    for itime in 1:nt
        # Create linked list
        create_CLL!(head, next, particles, pb, CLL)

        # Density and pressure EOS
        compute_density!(PoiseuilleFlow(), params, particles, pb, CLL, head, next, α)
        density_to_pressure!(particles, params)
        
        # Reset acceleration
        reset_particle_arrays!(particles, particles.ax, 0.0)
        reset_particle_arrays!(particles, particles.ay, 0.0)

        # Compute pressure gradient acceleration
        compute_pressuregradient_acceleration!(params, particles, pb, CLL, head, next, α)

        # Compute acceleration from viscous term
        compute_viscousterm_acceleration!(MorrisFD(),params, particles, pb, CLL, head, next, α)

        # Compute driving force
        driving_force!(particles, params)

        # Compute stable time step
        dt = compute_stable_dt(PoiseuilleFlow(), params, particles)

        # Update velocities and coordinates
        integrate_symplectic!(particles, params, dt)

        # Update time
        time += dt

        # Print to screen
        if itime % ncheck == 0
            println("vx min/max (fluid): ", extrema(particles.vx[particles.idx_fluid]))
            println("mean ρ (fluid): ", mean(particles.ρ[particles.idx_fluid]))
        end

        # Visualize
        if itime % nviz == 0
            for idx_spl in axes(u_sph_mat, 2)
                x_spl = 0.0 + idx_spl * dx_spl
                _, ys, u_sph = sample_along_segment(particles, params, particles.vx; p0 = (x_spl, 0), p1 = (x_spl, params.Ly))
                u_sph_mat[:, idx_spl] .= u_sph
            end
            u_sph_avg = mean(u_sph_mat, dims = 2)
            u_ana = analytical_poiseuille(ys, Ly, params.a_drive, params.ν, time)
            update_plot_Makie!(Plot1D(), sc1, particles.x[idxF], particles.y[idxF])
            update_plot_Makie!(Plot1D(), sc2, ys, u_sph_avg)
            update_plot_Makie!(Plot1D(), ln1, ys, u_ana)
            sc1.color = particles.vx[idxF]
            ax1.title = @sprintf "Time: %.5f s" time
            display(f)

            # Save the figure
            fname = @sprintf "./png/2D_poiseuille_%05d.png" itime
            save(fname, f, px_per_unit = 4)
        end
    end

    # Return
    return nothing
end

# Run main function
main_2D_poiseuille_flow();
