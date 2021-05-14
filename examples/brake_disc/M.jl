module Mod

using LinearAlgebra
using SparseArrays

function ABAT!(C, A, B)
    M, N = size(A)
    @assert N == size(B, 1)
    @assert N == size(B, 2)
    @inbounds for k in 1:N, l in 1:N
        b = B[k, l]
        for c in 1:M
            for r in 1:M
                C[r, c] += A[r, k] * (b * A[c, l])
            end
        end
    end
    return C
end

function spABAT!(C, A, B)
    M, N = size(A)
    @assert N == size(B, 1)
    @assert N == size(B, 2)
    I, J, V = findnz(B)
    for v in 1:length(V)
        k, l = I[v], J[v]
        b = V[v]
        @inbounds for c in 1:M
            f = (b * A[c, l])
            @simd for r in 1:M
                C[r, c] += A[r, k] * f
            end
        end
    end
    return C
end

end

using BenchmarkTools
using LinearAlgebra
using SparseArrays
using .Mod

function test()
    M, N = 84, 130
    A = rand(M, N)
    B = rand(N, N)
    C = zeros(M, M)
    Mod.ABAT!(C, A, B)
    C2 = A * B * A'
    @show norm(C - C2)

    @btime Mod.ABAT!($C, $A, $B)
    @btime $C2 .= $A * $B * $A'
end

#test()


function sptest()
    M, N = 84, 130
    A = rand(M, N)
    B = sprand(N, N, 0.001)
    @show nnz(B)
    C = zeros(M, M)
    Mod.spABAT!(C, A, B)
    C2 = A * B * A'
    @show norm(C - C2)
    @show typeof(C2)

    @btime Mod.spABAT!($C, $A, $B)
    @btime $C2 .= $A * $B * $A'
end

#sptest()


function sptest2()
    N, M = (138720, 8064)
    A = sprand(M, N, 0.01)
    B = sprand(N, N, 0.0012/2)
    B = (B + B') / 2
    @show nnz(A) / prod(size(A))
    @show nnz(B) / prod(size(B))   
    @btime C2 = $A * $B * $A'
    @btime C2 = $A * ($B * $A')
    @btime C2 = ($A * $B) * $A'
    true
end

sptest2()


function sptest3()
    N, M = (138720, 8064)
       A = sprand(N, M, 0.01)
       B = sprand(N, N, 0.0012/2)
    B = (B + B') / 2
    @show nnz(A) / prod(size(A))
    @show nnz(B) / prod(size(B))    
    @btime C2 = $A' * $B * $A
    @btime C2 = $A' * ($B * $A)
    @btime C2 = ($A' * $B) * $A
    true
end

sptest3()
