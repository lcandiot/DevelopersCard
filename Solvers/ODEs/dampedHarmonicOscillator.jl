# Solving the second order differential equation for a damped harmonic oscillator (F = m*a = m*d2x/dt2 = -k*x - c*dx/dt)
using OrdinaryDiffEq, CairoMakie, MathTeXEngine, Printf

# Define the harmonic oscillator function
function dampedHarmonicOscillator(d²udt², dudt, u, p, t) # Order of arguments matters
    d²udt² .= -p[1]^2 * u - 2.0 * p[2] * p[1] * dudt # d²xdt² = -ω^2 - 2 * ζ * ω * dxdt
end

function main_dampedHarmonicOscillator(; save_png :: Bool = false)

    # Physical parameters
    m     = 1.0                     # Mass [kg]
    k     = 1.0                     # Spring stiffness [N/m]
    c     = 0.4                    # Damping constant [N*s/m]
    tspan = (0.0, 10.0)             # Time span to solve [s]

    # Derive physical parameters
    ω = sqrt(k/m)
    ζ = c / (2.0 * sqrt(m*k))
    p = (ω, ζ)                      # Parameters stored in tuple

    # Initial conditions
    x₀   = [0.0]                    # Position array of the mass [m]
    dxdt = [1.0]                    # Velocity array of the mass [m/s]

    # Define problem
    prob = SecondOrderODEProblem(dampedHarmonicOscillator, dxdt, x₀, tspan, p)

    # Solve the problem and extract solutions
    sol = solve(prob, DPRKN6())
    x   = [sol.u[i][2] for i in eachindex(sol.t)]
    v   = [sol.u[i][1] for i in eachindex(sol.t)]
    t   = [sol.t[i]    for i in eachindex(sol.t)]

    # Visualize result
    fig = Figure(size = (600, 600))
    ax1 = Axis(fig[1,1], xlabel = L"$t$ [s]", ylabel = L"$x$ [m]", limits = (tspan[1], tspan[2], -max(x₀[1], dxdt[1]), max(x₀[1], dxdt[1])))
    for idx in eachindex(sol.t)
        pl1 = lines!(ax1, t[1:idx],  x[1:idx], color = :steelblue4, label = L"x")
        pl2 = lines!(ax1, t[1:idx],  v[1:idx], color = :goldenrod4, label = L"d$x$/d$t$")
        sc1 = scatter!(ax1, t[idx], x[idx],   color = :steelblue4)
        sc2 = scatter!(ax1, t[idx], v[idx],   color = :goldenrod4)
        if idx == 1
            axislegend(ax1)
        end
        display(fig)
        if save_png
            fname = @sprintf "./png/dampedHarmonicOscillator/dhOscillator_%03d.png" idx
            save(fname, fig, px_per_unit = 4)
        end
    end

    # Return
    return sol
end

# Run main function
sol = main_dampedHarmonicOscillator(; save_png = true);