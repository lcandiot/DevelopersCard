using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using CairoMakie, KernelAbstractions
using Enzyme
using Symbolics

# Flux kernel
@kernel inbounds = true function compute_flux!(q, C, k, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O

    q.x[Idx...] = -k * ∂x(C, g, Idx...)
    q.y[Idx...] = -k * ∂y(C, g, Idx...)
end

# Update kernel
# @kernel inbounds = true function update_C!(q, C, Δt, g :: StructuredGrid, O)
#     Idx = @index(Global, NTuple)
#     Idx += O
#     C[Idx...] -= Δt * divg(q, g, Idx...)
# end

@kernel inbounds = true function update_old!(C, C_old, O)
    Idx = @index(Global, NTuple)
    Idx += O

    C_old[Idx...] = C[Idx...]
end

@kernel inbounds = true function res_C!(q, C, C_old, R_C, Δt, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O
    R_C[Idx...] = (C[Idx...] - C_old[Idx...]) / Δt - divg(q, g, Idx...)
end

@kernel inbounds = true function JVP_C!(q, C, C_old, R_C, dC, dR_C, Δt, g :: StructuredGrid, O)
    Idx = @index(Global, NTuple)
    Idx += O
    F(C, C_old) = (C - C_old) / Δt - divg(q, g, Idx...)
    # To do: replace with autodiff_deferred but that hangs on the CPU
    @show F(C[Idx...], C_old[Idx...])
    @show Res = autodiff(set_runtime_activity(Forward), Const(F), Duplicated, Duplicated(C[Idx...], dC[Idx...]), Const(C_old[Idx...]))
    R_C[Idx...] = Res.val
    dR_C[Idx...] = Res.dval
    # R_C[Idx...] = (C[Idx...] - C_old[Idx...]) / Δt - divg(q, g, Idx...)
end

function res_C_wrap!(arch, grid, q, k, C, C_old, R_C, Δt, launch)
    launch(arch, grid, res_C! => (q, k, C, C_old, R_C, Δt, grid); bc = batch(grid, C => Neumann()) )
    return nothing
end

function my_assemble_Jacobian!(arch, grid, q, k, C, dC, C_old, R_C, dR_C, Δt, launch)
    autodiff(set_runtime_activity(Forward), res_C_wrap!,
        Const(arch), Const(grid), Const(q), Const(k),
        Duplicated(C, dC), Const(C_old), Duplicated(R_C, dR_C),
        Const(Δt),
        Const(launch)
    )

    # function inner(C)
    #     res_C_wrap(arch, grid, q, k, C, C_old, R_C, Δt)
    #     return R_C
    # end
    # jacobian(Enzyme.Forward, inner, vec(C))
end

# Time step function
function time_loop(C, C_old, R_C, q, Δt, k, arch, grid, launch)
    
    # Compute physics
    launch(arch, grid, compute_flux! => (q, C, k, grid))

    # Residual and Jacobian
    # launch(arch, grid, JVP_C! => (q, C, C_old, R_C, dC, dR_C, Δt, grid))
    # @show dR_C
    # @show dC

    launch(arch, grid, res_C! => (q, C, C_old, R_C, Δt, grid); bc = batch(grid, C => Neumann()) )
    # my_assemble_Jacobian!(arch, grid, q, k, C, dC, C_old, R_C, dR_C, Δt, launch)
    # @show dR_C |> Array
    # @show dC |> Array

    return nothing

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
    @show Δt = minimum(spacing(grid))^2 / k / ndims(grid) / 2.1

    # Plotting
    nviz = 10

    # Allocate fields
    C     = Field(backend, grid, Center())
    dC    = Field(backend, grid, Center())
    C_old = Field(backend, grid, Center())
    R_C   = Field(backend, grid, Center())
    dR_C  = Field(backend, grid, Center())
    q     = VectorField(backend, grid)

    # Initial condition
    my_gaussian(x, y, x0, y0, σ) = exp(-((x - x0)^2 / 2σ^2 ) - ((y - y0)^2 / 2σ^2 ))
    set!(C, grid, my_gaussian; parameters = (x0 = 0.0, y0 = 0.0, σ = 0.1))
    @show C |> Array
    set!(dC, grid, (x, y) -> 1.0) 
    set!(dR_C, grid, (x, y) -> 0.0) 
    set!(R_C, grid, (x, y) -> 0.0) 
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
    nt = 100
    # for idx_time in 1:nt

    #     # Compute physics
    #     launch(arch, grid, update_old! => (C, C_old))
    #     launch(arch, grid, compute_flux! => (q, C, k, grid))

    #     # Residual and Jacobian
    #     # launch(arch, grid, JVP_C! => (q, C, C_old, R_C, dC, dR_C, Δt, grid))
    #     # @show dR_C
    #     # @show dC

    #     # launch(arch, grid, res_C! => (q, k, C, C_old, R_C, Δt, grid); bc = batch(grid, C => Neumann()) )
    #     my_assemble_Jacobian!(arch, grid, q, k, C, dC, C_old, R_C, dR_C, Δt, launch)
    #     @show dR_C |> Array
    #     @show dC |> Array


    #     # Visualize
    #     if nviz % idx_time == 0
    #         hm1[3] = interior(C) |> Array
    #         display(fg1)
    #         sleep(1)
    #     end

    # end
    for idx_time in 1:nt
        @show idx_time
        launch(arch, grid, update_old! => (C, C_old))
        # time_loop(C, C_old, R_C, q, Δt, k, arch, grid, launch)
        # dC = Vector, dR_C = Jacobian-Vector-Product
        time_loop(C, C_old, R_C, q, Δt, k, arch, grid, launch)
        autodiff(
            set_runtime_activity(Forward), 
            time_loop, 
            Const, 
            Duplicated(C, R_C), 
            Const(C_old), 
            DuplicatedNoNeed(R_C, dR_C), 
            DuplicatedNoNeed(q, Enzyme.make_zero(q)),
            Const(Δt),
            Const(k),
            Const(arch),
            Const(grid),
            Const(launch)
        )
        @show dR_C |> Array
        # @show dC |> Array
        # @show C |> Array
        C .+= dR_C
        # @show C |> Array
        set!(dR_C, grid, (x, y) -> 0.0)
            #     # Visualize
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

n = 50

# Run main
main_diffusion2D(; nxy = (n, n) .- 2)