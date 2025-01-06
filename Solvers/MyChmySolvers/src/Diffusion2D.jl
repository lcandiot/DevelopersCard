using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using CairoMakie, KernelAbstractions

# Flux kernel
@kernel inbounds = true function compute_flux!(q, C, k, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O

    q.x[Idx...] = -k * ∂x(C, g, Idx...)
    q.y[Idx...] = -k * ∂y(C, g, Idx...)
end

# Update kernel
@kernel inbounds = true function update_C!(q, C, Δt, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O
    C[Idx...] -= Δt * divg(q, g, Idx...)
end

@views function main_diffusion2D(backend=CPU(); nxy = (126, 126))
    
    # Set architecture
    arch = Arch(backend)

    # Initialize grid
    grid   = UniformGrid(arch; origin = (-1, -1), extent = (2, 2), dims = nxy)
    launch = Launcher(arch, grid; outer_width = nothing)

    # Physics
    k = 1.0

    # Time stepping
    Δt = minimum(spacing(grid))^2 / k / ndims(grid) / 2.1

    # Plotting
    nviz = 10

    # Allocate fields
    C = Field(backend, grid, Center())
    q = VectorField(backend, grid)

    # Initial condition
    my_gaussian(x, y, x0, y0, σ) = exp(-((x - x0) / σ )^2 - ((y - y0) / σ )^2)
    set!(C, grid, my_gaussian; parameters = (x0 = 0.0, y0 = 0.0, σ = 0.1))
    bc!(arch, grid, C => Neumann())
    fg1 = Figure(size = (400, 400))
    ax1 = Axis(
        fg1[1,1][1,1], 
        xlabel = "x", 
        ylabel = "y",
        aspect = 1.0
    )
    hm1 = heatmap!(
        ax1,
        centers(grid)...,
        interior(C) |> Array,
    )
    Cb = Colorbar(fg1[1,1][1,2], hm1)
    display(fg1)

    # Time loop
    nt = 1000
    for idx_time in 1:nt

        # Compute physics
        launch(arch, grid, compute_flux! => (q, C, k, grid))
        launch(arch, grid, update_C! => (q, C, Δt, grid); bc = batch(grid, C => Neumann()) )

        # Visualize
        if nviz % idx_time == 0
            hm1[3] = interior(C) |> Array
            display(fg1)
            sleep(1)
        end
    end

    # Synchronize backend
    KernelAbstractions.synchronize(backend)


    # Return
    return nothing

end

n = 126

# Run main
main_diffusion2D(; nxy = (n, n) .- 2)