# Solving the second order differential equation for a harmonic oscillator (F = m*a = m*d2x/dt2 = -k*x)
using OrdinaryDiffEq, CairoMakie, MathTeXEngine

# Define the harmonic oscillator function
function harmonicOscillator(d²udt², dudt, u, ω, t) # Order of arguments matters
    d²udt² .= -ω^2 * u
end

function main_harmonicOscillator()

    # Physical parameters
    m = 1.0                     # Mass [kg]
    k = 1.0                     # Spring stiffness [N/m]
    tspan = (0.0, 10.0)         # Time span to solve [s]

    # Derive physical parameters
    ω = sqrt(k/m)

    # Initial conditions
    x₀   = [0.0]              # Position array of the mass [m]
    dxdt = [1.0]              # Velocity array of the mass [m/s]

    # Define problem
    prob = SecondOrderODEProblem(harmonicOscillator, dxdt, x₀, tspan, ω)

    # Solve the problem
    sol = solve(prob, DPRKN6())

    # Visualize result
    fig = Figure(size = (500, 500))
    ax1 = Axis(fig[1,1], xlabel = L"$t$ [s]", ylabel = L"$x$ [m]")
    pl1 = plot!(ax1, sol)
    display(fig)

    # Return
    return nothing
end

# Run main function
main_harmonicOscillator();