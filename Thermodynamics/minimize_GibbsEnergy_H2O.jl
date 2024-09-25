# This script illustrates how to calculate stable phases of H2O by minimizing Gibbs energy
using Zygote, CairoMakie

# Define main function
function minimize_GibbsEnergy_H2O()
    
    # Constants
    R         = 8.314  # Gas constant (J/mol·K)
    T0        = 273.15  # Reference temperature (K)
    P0        = 1e5    # Reference pressure (Pa)
    ΔH_vap    = 44.01e3  # Heat of vaporization (J/mol)
    ΔH_fus    = 6.01e3   # Heat of fusion (J/mol)
    ΔV_liquid = 1e-5  # Molar volume of liquid water (m³/mol) at 298.15 K

    # Return
    return nothing

end

#--------------------------
# Run main
minimize_GibbsEnergy_H2O();

# Gibbs Free Energy functions
function gibbs_vapor(T, P, ΔH_vap)
    ΔG_0_vapor = -228.6e3  # Standard Gibbs free energy for vapor at 298.15 K (J/mol)
    return ΔG_0_vapor + ΔH_vap * ((T + 273.15) - T0) / T0 + R * (T + 273.15) * log(P / P0)
end

function gibbs_solid(T, P)
    ΔG_0_ice = -333.6e3  # Approximate Gibbs free energy for ice at 0°C and 1 bar (J/mol)
    ΔS_ice = 41.1        # Entropy change for ice (J/mol·K)
    return ΔG_0_ice + ΔS_ice * ((T + 273.15) - T0)
end

function gibbs_liquid(T, P)
    ΔG_0_liquid = -237.1e3  # Gibbs free energy for liquid at 298.15 K (J/mol)
    ΔH_liquid = ΔH_fus + ΔH_vap * (T + 273.15 - T0) / T0
    ΔV_liquid = 18.015 * (1e-6)  # Molar volume (m³/mol) adjusted for temperature
    return ΔG_0_liquid - ΔH_liquid * (T + 273.15 - T0) / T0 + ΔV_liquid * (P * 1e5 - P0)
end

# Function to calculate the stable phase using linear programming
function Gibbs_energy_total(T, P, X_vapor, X_liquid, X_solid)
    g_vapor  = gibbs_vapor(T, P)
    g_liquid = gibbs_liquid(T, P)
    g_solid  = gibbs_solid(T, P)

    g_total = g_vapor * X_vapor + g_liquid * X_liquid + g_solid * X_solid


end

# Define the temperature and pressure intervals
T_range = -4.0:1.0:120.0  # Temperature in °C
P_range = 0.0:0.1:3.0     # Pressure in bar

# Create a grid to store the stable phases
phases = Dict{Tuple{Float64, Float64}, Symbol}()

# Loop over the grid to calculate the stable phase
for T in T_range
    for P in P_range
        phases[(T, P)] = stable_phase(T, P)
    end
end

# Extract values for plotting
T_values = [key[1] for key in keys(phases)]
P_values = [key[2] for key in keys(phases)]
phase_values = [phases[key] for key in keys(phases)]

# Convert phase values to numeric values for coloring
phase_numeric = map(phase -> phase == :vapor ? 1 : phase == :ice ? 2 : 3, phase_values)

# Check for potential empty values
if isempty(phase_numeric)
    error("The phase_numeric array is empty, which indicates an issue in the phase calculation or storage process.")
end

# Create the phase diagram using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Temperature (°C)", ylabel = "Pressure (bar)", title = "H2O Phase Diagram")

scatter!(ax, T_values, P_values, color = phase_numeric, colormap = :viridis)

# Add a color legend
Legend(fig[1, 2], ax, ["Vapor", "Ice", "Liquid"], color = [:green, :blue, :yellow])

fig
