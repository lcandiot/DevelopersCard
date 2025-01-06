using Enzyme
using LinearAlgebra

function main_JVP()
    
    u = [4.0, 5.0, 7.0]
    v = [1.0, 1.0, 1.0]

    y = [0.0, 0.0, 0.0]

    F(u) = u.^2
    @show JVP_FD = JVP_Finite_Diff(F, u, v)

    JVP!(y, F, u, v) = Enzyme.autodiff(Forward, 
    (temp, v) -> (temp .= F(v); nothing), Const, 
    DuplicatedNoNeed(zero(y), y), Duplicated(u, v))

    JVP_Enz = JVP!(y, F, u, v)

    @show y
    # Return
    return nothing

end

function JVP_Finite_Diff(F, u, v)
    λ = 10e-6
    δ = λ * (λ + norm(u, Inf)/norm(v,Inf))

    (F(u + δ .* v) - F(u)) ./ δ
end

# Run main
main_JVP();