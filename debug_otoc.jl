#!/usr/bin/env julia

using Printf
using PauliOperators
using LinearAlgebra
using DBF

include("./src/observables.jl")
include("./src/otoc_pauli.jl")
include("./src/otoc_exact.jl")

# Simple test case: 2 qubits, simple Hamiltonian, short time
N = 2
H = DBF.af_heisenberg(N, 1.0, 1.0, 1.0)
println("Hamiltonian:")
display(H)
println()

# Operators at sites 1 and 2
O1 = magz(N, 1)
O2 = magz(N, 2)

println("O1 (site 1):")
display(O1)
println()

println("O2 (site 2):")
display(O2)
println()

# At t=0, compute commutator
comm_0 = O2 * O1 - O1 * O2
println("Commutator [O2, O1] at t=0:")
display(comm_0)
println("Norm squared of commutator at t=0: ", inner_product(comm_0, comm_0))
println()

# Small time
t = 0.1

# Method 1: Trotter
otoc_trotter, _ = infinite_temp_otoc_pauli(O1, O2, H, t; dt=0.001)
println("OTOC (Trotter, dt=0.001) at t=$t: $otoc_trotter")

# Method 2: Exact
otoc_exact, _ = infinite_temp_otoc_exact(O1, O2, H, t)
println("OTOC (Exact) at t=$t: $otoc_exact")

println("\nDifference: ", abs(otoc_trotter - otoc_exact))
println("Relative difference: ", abs(otoc_trotter - otoc_exact) / otoc_exact)

# Let's also check the evolution of O2 directly
println("\n" * "="^70)
println("Checking O2 evolution...")
println("="^70)

# Exact evolution of O2
H_mat = Matrix(H)
O2_mat = Matrix(O2)
U = exp(-1im * t * H_mat)
O2t_exact_mat = U' * O2_mat * U

println("\nO2(t) using exact matrix exponentiation:")
println("(showing first few elements)")
println(O2t_exact_mat[1:min(4,size(O2t_exact_mat,1)), 1:min(4,size(O2t_exact_mat,2))])

# Trotter evolution of O2
O2t_trotter = deepcopy(O2)
n_steps = Int(round(t / 0.001))
for step in 1:n_steps
    global O2t_trotter
    for (p, c) in H
        θ = 2 * real(c) * 0.001  # Corrected: factor of 2 because evolve uses θ/2
        O2t_trotter = evolve(O2t_trotter, p, θ)
    end
end

println("\nO2(t) using Trotter (dt=0.001):")
O2t_trotter_mat = Matrix(O2t_trotter)
println("(showing first few elements)")
println(O2t_trotter_mat[1:min(4,size(O2t_trotter_mat,1)), 1:min(4,size(O2t_trotter_mat,2))])

println("\nDifference in O2(t) matrices:")
diff_mat = O2t_exact_mat - O2t_trotter_mat
println("Frobenius norm of difference: ", norm(diff_mat))
println("Relative error: ", norm(diff_mat) / norm(O2t_exact_mat))
