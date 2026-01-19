using DBF
using PauliOperators
using LinearAlgebra
using Printf
using Random

"""
Demonstration: 1D Hubbard model (4 sites) with half-filling and particle number preservation

For demonstration, we use a 1D system:
- 4 physical sites
- 2 spin orbitals per site (up and down)
- Total of 8 qubits

Half-filling means 2 electrons.
"""

println("="^80)
println("1D Fermi-Hubbard Model (4 sites) with Half-Filling")
println("Using Particle Number Preservation in DBF Ground State")
println("="^80)
println()

# System parameters
L = 4  # Number of sites
N_qubits = 2 * L  # 8 qubits (2 per site for spin up/down)

println("System parameters:")
println("  Number of sites: $(L)")
println("  Total qubits (spin up + down): $(N_qubits)")
println()

# Hubbard model parameters
t = 1.0   # Hopping amplitude
U = 2.0   # On-site interaction

println("Hubbard model parameters:")
println("  Hopping (t): $(t)")
println("  On-site interaction (U): $(U)")
println()

# Build the Hamiltonian
println("Building 1D Fermi-Hubbard Hamiltonian...")
H = DBF.hubbard_model_1D(L, t, U)
println("  Number of terms in Hamiltonian: ", length(H))
println("  Hamiltonian norm: ", norm(H))
println()

# Create the particle number operator
println("Creating particle number operator...")
N̂ = DBF.particle_number_operator(N_qubits)
println("  Number of terms in N̂: ", length(N̂))
println()

# Create an initial state - a simple product state with 2 particles
# Let's use |01010000> which has 2 particles
target_particles = 2
ψ = Ket{N_qubits}(0b00010001)  # particles at positions 1 and 5 (0-indexed bits 0 and 4)

println("Initial state:")
println("  State: |$(bitstring(Int128(ψ.v))[end-N_qubits+1:end])>")

# Verify particle number
particle_num = real(expectation_value(N̂, ψ))
println("  Particle number: $(particle_num)")
println("  Initial energy: ", real(expectation_value(H, ψ)))
println("  Initial variance: ", real(DBF.variance(H, ψ)))
println()

# Check if the Hamiltonian preserves particle number
println("Checking Hamiltonian properties...")
println("  Does H preserve particle number? ", DBF.preserves_particle_number(H))
println()

# Run DBF ground state optimization WITHOUT particle number preservation
println("="^80)
println("First: Running WITHOUT particle number preservation")
println("="^80)
println()

res_no_pn = DBF.dbf_groundstate(deepcopy(H), ψ,
                n_body=1,
                max_iter=10,
                conv_thresh=1e-6,
                evolve_coeff_thresh=1e-10,
                grad_coeff_thresh=1e-10,
                preserve_particle_number=false,
                verbose=1)

println()
println("="^80)
println("Now: Running WITH particle number preservation")
println("="^80)
println()

# Run DBF ground state optimization WITH particle number preservation
res_with_pn = DBF.dbf_groundstate(deepcopy(H), ψ,
                n_body=1,
                max_iter=10,
                conv_thresh=1e-6,
                evolve_coeff_thresh=1e-10,
                grad_coeff_thresh=1e-10,
                preserve_particle_number=true,
                verbose=1)

println()
println("="^80)
println("Comparison of Results:")
println("="^80)
println()

# Check particle number after transformation
final_particle_num_no_pn = real(expectation_value(N̂, ψ))
final_particle_num_with_pn = real(expectation_value(N̂, ψ))

println("WITHOUT Particle Number Preservation:")
println("  Number of generators: ", length(res_no_pn["generators"]))
println("  Final energy: ", real(res_no_pn["energies"][end]))
println("  Energy lowering: ", real(res_no_pn["energies"][1] - res_no_pn["energies"][end]))
println("  Final variance: ", real(res_no_pn["variances"][end]))
println("  Final H terms: ", length(res_no_pn["hamiltonian"]))
if length(res_no_pn["generators"]) > 0
    n_pn_preserving = sum(DBF.preserves_particle_number.(res_no_pn["generators"]))
    println("  Generators preserving particle number: $(n_pn_preserving)/$(length(res_no_pn["generators"]))")
end
println()

println("WITH Particle Number Preservation:")
println("  Number of generators: ", length(res_with_pn["generators"]))
println("  Final energy: ", real(res_with_pn["energies"][end]))
println("  Energy lowering: ", real(res_with_pn["energies"][1] - res_with_pn["energies"][end]))
println("  Final variance: ", real(res_with_pn["variances"][end]))
println("  Final H terms: ", length(res_with_pn["hamiltonian"]))
if length(res_with_pn["generators"]) > 0
    n_pn_preserving = sum(DBF.preserves_particle_number.(res_with_pn["generators"]))
    println("  All generators preserve particle number: ", all(DBF.preserves_particle_number.(res_with_pn["generators"])))
end
println()

# Show some example generators
if length(res_with_pn["generators"]) > 0
    println("Example particle-number-preserving generators:")
    for (i, gen) in enumerate(res_with_pn["generators"][1:min(5, length(res_with_pn["generators"]))])
        num_flips = count_ones(gen.x)
        println("  $(i): $(string(gen)) ($(num_flips) X/Y operators)")
    end
end

println()
println("="^80)
println("Demonstration complete!")
println("="^80)
