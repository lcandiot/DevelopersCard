using Enzyme
using LinearAlgebra, LinearOperators
using CairoMakie
using Krylov

function JVP_Finite_Diff(F, u, v)
    λ = 10e-6
    δ = λ * (λ + norm(u, Inf)/norm(v,Inf))

    (F(u + δ .* v) - F(u)) ./ δ
end

function compute_physics!(C, qx, Res, C_old, k, dx, dt)

    # Get size
    ncx = size(C)[1]

    # # Compute flux
    # for idx in 1:ncx-1
    #     qx[idx + 1] = -k * (C[idx + 1] - C[idx]) / dx
    # end

    # # Compute residual
    # for idx in eachindex(C)
    #     if idx > 1 && idx < ncx
    #         Res[idx] = -(C[idx] - C_old[idx]) / dt - (qx[idx + 1] - qx[idx]) / dx
    #     end
    # end

    # Compute residual
    for idx in eachindex(C)
        if idx > 1 && idx < ncx
            Res[idx] = -(C[idx] - C_old[idx]) / dt + k * (C[idx - 1] - 2C[idx] + C[idx + 1]) / dx^2
        else
            # Res[idx] = 0.0
            C[idx] = 0.0
        end
    end
end

function update_old!(A, A_old)
    for idx in eachindex(A)
        A_old[idx] = A[idx]
    end
end

# function JVP!(y, F, u, v) 
#     Enzyme.autodiff(
#         Forward, 
#         (temp, v) -> (temp .= F(v); nothing), 
#         Const, 
#         DuplicatedNoNeed(zero(y), y), # here we store the JVP
#         Duplicated(u, v)
#     )
#     return nothing
# end

function JVP!(C, ΔC, qx, dqx, Res, jvp, C_old, dC_old, k, dx, dt)
    Enzyme.make_zero!(dqx)
    Enzyme.make_zero!(jvp)
    Enzyme.make_zero!(dC_old)
    Enzyme.autodiff(
        Forward,
        compute_physics!,
        Const,
        DuplicatedNoNeed(C, ΔC),
        Const(qx),
        DuplicatedNoNeed(Res, jvp),
        Const(C_old),
        Const(k),
        Const(dx),
        Const(dt)
    )
	return nothing
end

"""
Calculate the Jacobian-Transpose Vector Product in-place by updating `y`.
"""
function JᵀVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt)
    # jvp .= 0 # Enzyme expects y to be zero
    Enzyme.make_zero!(dqx)
    Enzyme.make_zero!(jvp)
    Enzyme.make_zero!(dC_old)
    Enzyme.autodiff(
        Enzyme.Reverse,
        compute_physics!,
        Const,
        DuplicatedNoNeed(C, jvp),
        DuplicatedNoNeed(qx, dqx),
        DuplicatedNoNeed(Res, copy(ΔC)),
        DuplicatedNoNeed(C_old, dC_old),
        Const(k),
        Const(dx),
        Const(dt)
    )
	return nothing
end

function main_JVP()
    
    # Physics
    Lx = 1.0
    k  = 1.0

    # Numerics
    nx = 101
    dx = Lx / (nx - 1)
    dt = dx^2 / k / 2.1
    nt = 100

    # Initialisation
    x      = LinRange(-Lx / 2.0 + dx / 2.0, Lx / 2.0 - dx / 2.0, nx - 1)
    C      = exp.(-(x).^2 ./ 0.1.^2)
    C_old  = copy(C)
    dC     = [0.0 for _ in 1:nx-1]
    dC_old = [0.0 for _ in 1:nx-1]
    qx     = [0.0 for _ in 1:nx  ]
    dqx    = [0.0 for _ in 1:nx  ]
    Res    = [0.0 for _ in 1:nx-1]
    dRes   = [rand() for _ in 1:nx-1]
    jvp    = [0.0 for _ in 1:nx-1]

    # Check initial configuration
    C_ini = copy(C)
    fg1 = Figure(size=(300, 300))
    ax1 = Axis(fg1[1,1], xlabel = "x", ylabel = "C")
    ln1 = lines!(ax1, x, C_ini)
    display(fg1)

    # Time loop
    for idx_time in 1:nt
        # Update old
        update_old!(C, C_old)

        # # Get size
        # ncx = size(C)[1]

        # # Compute flux
        # for idx in 1:ncx-1
        #     qx[idx + 1] = -k * (C[idx + 1] - C[idx]) / dx
        # end

        # # Compute residual
        # for idx in eachindex(C)
        #     if idx > 1 && idx < ncx
        #         C[idx] -= dt * (qx[idx + 1] - qx[idx]) / dx
        #     end
        # end

        # Compute physics while calculating the JVP
        compute_physics!(C, qx, Res, C_old, k, dx, dt)

        Res_copy = -copy(Res)

        opJ = LinearOperator(Float64, nx-1, nx-1, false, false, 
            (jvp, ΔC) -> JVP!(C, ΔC, qx, dqx, Res, jvp, C_old, dC_old, k, dx, dt),
            (jvp, ΔC) -> JᵀVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt),
            (jvp, ΔC) -> JᵀVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt)
        )
        # opJ = LinearOperator(Float64, nx-1, nx-1, false, false, 
        #     (jvp, ΔC) -> JVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt),
        #     (jvp, ΔC) -> JᵀVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt),
        #     (jvp, ΔC) -> JᵀVP!(C, qx, dqx, Res, ΔC, jvp, C_old, dC_old, k, dx, dt)
        # )
        # Enzyme.autodiff(
        #     Forward,
        #     compute_physics!,
        #     Const,
        #     DuplicatedNoNeed(zero(jvp), jvp), # here we store the JVP
        #     Duplicated(C, dC),
        #     Duplicated(qx, dqx),
        #     Duplicated(Res, dRes),
        #     Duplicated(C_old, dC_old),
        #     Const(k),
        #     Const(dx),
        #     Const(dt)
        # )
        #@show Res
        # rand
        # Solve the system of equations
        # jvp = zeros(nx-1)
        # for _ in 1:3
        #     JVP!(C, dC, qx, dqx, Res, jvp, C_old, dC_old, k, dx, dt)
        #     @show norm(jvp)
        # end

        (incC, stats) = gmres(opJ, Res_copy, C)
        C .+= incC
        # @show C
        ln2 = lines!(ax1, x, C)
        display(fg1)
        # @show stats

    end
    # u = [4.0, 5.0, 7.0]
    # v = [1.0, 1.0, 1.0]
    # y = [0.0, 0.0, 0.0]

    # F(u) = u.^2
    # @show JVP_FD = JVP_Finite_Diff(F, u, v)


    # JVP_Enz = JVP!(y, F, u, v)

    # @show y
    # Return
    return nothing

end

# Run main
main_JVP();