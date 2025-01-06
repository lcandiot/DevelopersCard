using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using Plots, KernelAbstractions

# Flux kernel
@kernel inbounds = true function compute_flux!(q, C, k, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O

    q.x[Idx...] = -k * ∂x(C, g, Idx...)
end

# Update kernel
@kernel inbounds = true function update_C!(q, C, Δt, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O
    C[Idx...] -= Δt * divg(q, g, Idx...)
end

@views function main_diffusion1D(backend=CPU(); nx = (126, ))
    
    # Set architecture
    arch = Arch(backend)

    # Initialize grid
    grid   = UniformGrid(arch; origin = (-1, ), extent = (2, ), dims = nx)
    launch = Launcher(arch, grid; outer_width = (4, ))

    # Physics
    k = 1.0

    # Time stepping
    Δt = spacing(grid)[1]^2 / k / ndims(grid) / 2.1

    # Allocate fields
    C = Field(backend, grid, Center())
    q = VectorField(backend, grid)

    # Initial condition
    my_gaussian(x, x0, σ) = exp(-((x - x0) / σ )^2 )
    set!(C, grid, my_gaussian; parameters = (x0 = 0.0, σ = 0.1))
    bc!(arch, grid, C => Neumann())
    Plots.display(Plots.plot(centers(grid), interior(C) |> Array))

    # Time loop
    nt = 100
    for idx_time in 1:nt

        # Compute physics
        launch(arch, grid, compute_flux! => (q, C, k, grid))
        launch(arch, grid, update_C! => (q, C, Δt, grid); bc = batch(grid, C => Neumann()) )

        # Visualize
        Plots.display(Plots.plot(centers(grid), interior(C) |> Array))
        sleep(0.1)
    end

    # Synchronize backend
    KernelAbstractions.synchronize(backend)


    # Return
    return nothing

end

n = 126

# Run main
main_diffusion1D(; nx = (n, ) .- 2)