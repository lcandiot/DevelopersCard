# Linear regression for temperature dependent density. Turns out that Flux.jl works best for values close to 1.0. It is therefore important to normalise or scale all data before training the neural network. This is done in this example.
using Pkg
Pkg.instantiate()
using Flux, Statistics, CairoMakie, MathTeXEngine

function density_Tdependent()

    # Real data
    T_train_min = 400.0 
    T_train_max = 500.0 
    T_test_min  = 600.0 
    T_test_max  = 700.0 
    T_train     = hcat(T_train_min:1.0:T_train_max...)
    T_test      = hcat(T_test_min:1.0:T_test_max...)
    noise_train = rand(Float64, length(T_train)) .* 0.1    # Random noise 10 %
    noise_test  = rand(Float64, length(T_test)) .* 0.1    
    α           = 1e-4                                      # Coefficient of thermal expansion [K]
    ρ0          = 2800.0                                    # Reference density [kg / m^3]

    # Scaling
    Tsc         = maximum(T_train)                 # Temperature scale
    ρsc         = copy(ρ0)                                  # Density scale
    α          *= Tsc
    ρ0         /= ρsc
    T_train   ./= Tsc
    T_test    ./= Tsc

    # Generate data
    ρ_train     = ρ0 .- α .* (T_train) 
    ρ_test      = ρ0 .- α .* (T_test ) 
    T_train   .+= noise_train' .* rand(Float64, 1)          # Randomise T to resemble real data
    T_test    .+= noise_test' .* rand(Float64, 1)

    # Convert to Float 32 for the Flux.jl
    T_train, T_test = Float32.(T_train), Float32.(T_test)
    ρ_train, ρ_test = Float32.(ρ_train), Float32.(ρ_test)

    # Visualize initial data
    f   = Figure()
    ax1 = Axis(f[1,1])
    scatter!(ax1, vec(T_train .* Tsc), vec(ρ_train .* ρsc), color = :blue)
    scatter!(ax1, vec(T_test .* Tsc), vec(ρ_test .* ρsc), color = :orange)
    display(f)

    # Solver settings
    max_epc = 1e7
    εtol    = 1e-15
    ε̇       = 1e0 / length(T_train)
    ecnt    = 1

    # Build neural network model with 1 dense layer for 1 input and 1 output
    model = Dense(1 => 1)

    # Train the model
    err = 2εtol; loss_hist = [loss_func(model, T_train, ρ_train)]; epc_hist  = [ecnt]
    while err > εtol && ecnt <= max_epc
        ecnt += 1
        if (ecnt % 1000) == 0
            println("Iter = $(ecnt); err = $(err)")
        end
        train_model(loss_func, model, T_train, ρ_train, ε̇)
        err = loss_func(model, T_train, ρ_train)
        if abs(loss_hist[ecnt - 1] - err) < εtol
            println("Iter = $(ecnt); err = $(err)")
            break
        end
        push!(loss_hist, err)
        push!(epc_hist, ecnt)

    end

    # Make a prediction
    ρ_predict = model(T_train)

    # Visualize
    f   = Figure()
    ax1 = Axis(f[1,1], xlabel = "T [C]", ylabel = "ρ [kg.m⁻³]", title = "Fitting a line using neural networks")
    ax2 = Axis(f[2,1], yscale = log10, xlabel = "No. epochs [ ]", ylabel = "err [ ]")
    sc1 = scatter!(ax1, T_train[:] .* Tsc, ρ_train[:] .* ρsc, label = "Data")
    sc2 = lines!(ax1, T_train[:] .* Tsc, ρ_predict[:] .* ρsc, label = "Fit", color = :tomato)
    li1 = lines!(ax2, Float32.(epc_hist[:]), Float32.(loss_hist[:]), label = "Loss function")
    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)
    display(f)

    # Save Figure
    save("./doc/png/density_Tdependent.png", f, px_per_unit=3)

    # Return
    return nothing

end

# Loss function
function loss_func(model, x, y)
    ŷ = model(x)
    Flux.mse(ŷ, y)
end

# Train function
function train_model(loss_func, model, x, y, ε̇)
    dLdm, _, _      = gradient(loss_func, model, x, y)
    @. model.weight = model.weight - ε̇ * dLdm.weight
    @. model.bias   = model.bias - ε̇ * dLdm.bias
end

# Run main
density_Tdependent()