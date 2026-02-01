using PauliOperators
using LinearAlgebra
using Printf
using DBF

"""
    infinite_temp_otoc_pauli(O1::PauliSum{N}, O2::PauliSum{N}, H::PauliSum{N}, t::Real; dt::Real=0.01) where {N}

Calculate the infinite temperature Out-of-Time-Ordered Correlator (OTOC) using Pauli propagation.

At infinite temperature, the OTOC is defined as:
    F(t) = Tr([O2(t), O1]†[O2(t), O1]) / (2^N * Tr(O1†O1) * Tr(O2†O2))

where O2(t) = exp(iHt) O2 exp(-iHt) is the time-evolved operator.

# Arguments
- `O1::PauliSum{N}`: First observable operator (at time 0)
- `O2::PauliSum{N}`: Second observable operator (to be time-evolved)
- `H::PauliSum{N}`: Hamiltonian for time evolution
- `t::Real`: Time for evolution
- `dt::Real`: Time step for evolution (default: 0.01)

# Returns
- `Float64`: The OTOC value at time t
- `PauliSum{N}`: The time-evolved O2(t)
"""
function infinite_temp_otoc_pauli(O1::PauliSum{N}, O2::PauliSum{N}, H::PauliSum{N}, t::Real; dt::Real=0.01) where {N}
    # Time-evolve O2 using Trotterization
    O2t = deepcopy(O2)
    n_steps = Int(round(abs(t) / dt))
    
    for step in 1:n_steps
        # Evolve by each term in the Hamiltonian
        for (p, c) in H
            θ = -real(c) * dt  # Negative for forward time evolution
            O2t = evolve(O2t, p, θ)
        end
    end
    
    # Calculate commutator [O2(t), O1]
    comm = O2t * O1 - O1 * O2t
    
    # Calculate OTOC = Tr(comm† * comm) / (2^N * Tr(O1†O1) * Tr(O2†O2))
    # For Pauli operators, Tr(A†A) = 2^N * ||A||^2 where ||A||^2 is the squared Frobenius norm
    comm_norm_sq = inner_product(comm, comm)
    O1_norm_sq = inner_product(O1, O1)
    O2_norm_sq = inner_product(O2, O2)
    
    # OTOC (normalized)
    if abs(O1_norm_sq) > 0 && abs(O2_norm_sq) > 0
        otoc = real(comm_norm_sq / (O1_norm_sq * O2_norm_sq))
    else
        otoc = 0.0
    end
    
    return otoc, O2t
end

"""
    run_otoc_pauli(site1_list::Vector{Int}, site2::Int, N::Int, H::PauliSum, 
                   t_final::Real, dt::Real; evolution_dt::Real=0.01)

Run OTOC calculations for multiple sites using Pauli propagation.

# Arguments
- `site1_list::Vector{Int}`: List of site indices for O1 operators
- `site2::Int`: Site index for O2 operator
- `N::Int`: Number of qubits in the system
- `H::PauliSum`: Hamiltonian for time evolution
- `t_final::Real`: Final time for OTOC calculation
- `dt::Real`: Time step for OTOC sampling
- `evolution_dt::Real`: Time step for Hamiltonian evolution (default: 0.01)

# Returns
- `Vector{Float64}`: OTOC values for all sites and times (flattened)
"""
function run_otoc_pauli(site1_list::Vector{Int}, site2::Int, N::Int, H::PauliSum, 
                        t_final::Real, dt::Real; evolution_dt::Real=0.01)
    O2_initial = magz(N, site2)
    all_otocs = Float64[]
    
    println("OTOC calculation using Pauli propagation started")
    
    for t in 0.0:dt:t_final
        otoc_sites = Float64[]
        
        for site1 in site1_list
            O1 = magz(N, site1)
            otoc, _ = infinite_temp_otoc_pauli(O1, O2_initial, H, t; dt=evolution_dt)
            push!(otoc_sites, real(otoc))
        end
        
        append!(all_otocs, otoc_sites)
        @printf("Finished calculation for t = %.4f\n", t)
    end
    
    println("OTOC calculation completed")
    return all_otocs
end

"""
    single_timestep_otoc_pauli(site1::Int, site2::Int, N::Int, H::PauliSum, dt::Real)

Calculate OTOC for a single time step using Pauli propagation.

# Arguments
- `site1::Int`: Site index for O1 operator
- `site2::Int`: Site index for O2 operator  
- `N::Int`: Number of qubits in the system
- `H::PauliSum`: Hamiltonian for time evolution
- `dt::Real`: Time step

# Returns
- `Float64`: The OTOC value at time dt
"""
function single_timestep_otoc_pauli(site1::Int, site2::Int, N::Int, H::PauliSum, dt::Real)
    O1 = magz(N, site1)
    O2 = magz(N, site2)
    
    otoc, _ = infinite_temp_otoc_pauli(O1, O2, H, dt; dt=dt/10)
    return real(otoc)
end
