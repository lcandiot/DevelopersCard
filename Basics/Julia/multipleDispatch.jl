using Random

# Define new main type
abstract type Operator end

# Define structs for different operations
struct MultiplyElements <: Operator end
struct DivideElements <: Operator end
struct MatMul <: Operator end

# Element wise multiplication operator
function operate_arrays(::MultiplyElements, A::AbstractArray, B::AbstractArray)
    return A .* B
end

# Elementwise division operator
function operate_arrays(::DivideElements, A::AbstractArray, B::AbstractArray)
    return A ./ B
end

# Matrix matrix multiplication operator
function operate_arrays(::MatMul, A::AbstractArray, B::AbstractArray)
    return A * B
end

# Define main
function run_main()
    # Randomness
    rng  = Random.MersenneTwister()
    seed = 123
    Random.seed!(rng, seed)

    # Initialise arrays
    A = rand(rng, 1:10, 3, 5)
    B = rand(rng, 1:10, 3, 5)
    C = rand(rng, 1:10, 5, 3)

    # Perform operations
    C_elmul  = operate_arrays(MultiplyElements(), A, B)
    C_eldiv  = operate_arrays(DivideElements(),   A, B)
    C_matmul = operate_arrays(MatMul(),           A, C)

    display(operate_arrays)
    display([A, B, C])
    display(C_elmul)
    display(C_eldiv)
    display(C_matmul)
end

# Run main
run_main()