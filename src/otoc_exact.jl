# OTOC calculation using exact matrix exponentiation
# Assumes PauliOperators, LinearAlgebra, Printf are already loaded

"""
    infinite_temp_otoc_exact(O1::PauliSum{N}, O2::PauliSum{N}, H::PauliSum{N}, t::Real) where {N}

Calculate the infinite temperature OTOC using exact matrix exponentiation.

This method uses direct matrix exponentiation to evolve operators:
    O2(t) = exp(iHt) O2 exp(-iHt)

This is exact (within numerical precision) but computationally expensive.
Use for small systems to validate approximate methods.

# Arguments
- `O1::PauliSum{N}`: First observable operator (at time 0)
- `O2::PauliSum{N}`: Second observable operator (to be time-evolved)
- `H::PauliSum{N}`: Hamiltonian for time evolution
- `t::Real`: Time for evolution

# Returns
- `Float64`: The OTOC value at time t
- `PauliSum{N}`: The time-evolved O2(t)
"""
function infinite_temp_otoc_exact(O1::PauliSum{N}, O2::PauliSum{N}, H::PauliSum{N}, t::Real) where {N}
    # Convert to matrices
    H_mat = Matrix(H)
    O2_mat = Matrix(O2)
    O1_mat = Matrix(O1)
    
    # Compute evolution operator U = exp(-iHt)
    # Note: In quantum mechanics, U(t) = exp(-iHt), so U†(t) = exp(iHt)
    U = exp(-1im * t * H_mat)
    U_dag = U'
    
    # Time-evolve O2: O2(t) = U† O2 U = exp(iHt) O2 exp(-iHt)
    O2t_mat = U_dag * O2_mat * U
    
    # Calculate commutator [O2(t), O1]
    comm_mat = O2t_mat * O1_mat - O1_mat * O2t_mat
    
    # Calculate OTOC = Tr(comm† * comm) / (Tr(O1†O1) * Tr(O2†O2))
    comm_norm_sq = tr(comm_mat' * comm_mat)
    O1_norm_sq = tr(O1_mat' * O1_mat)
    O2_norm_sq = tr(O2_mat' * O2_mat)
    
    # OTOC (normalized)
    if abs(O1_norm_sq) > 0 && abs(O2_norm_sq) > 0
        otoc = real(comm_norm_sq / (O1_norm_sq * O2_norm_sq))
    else
        otoc = 0.0
    end
    
    # Convert O2(t) back to PauliSum (optional, can be expensive)
    # For now, just return the OTOC value
    O2t = O2  # Placeholder, actual conversion would be expensive
    
    return otoc, O2t
end

"""
    run_otoc_exact(site1_list::Vector{Int}, site2::Int, N::Int, H::PauliSum, 
                   t_final::Real, dt::Real)

Run OTOC calculations using exact matrix exponentiation method.

# Arguments
- `site1_list::Vector{Int}`: List of site indices for O1 operators
- `site2::Int`: Site index for O2 operator
- `N::Int`: Number of qubits in the system
- `H::PauliSum`: Hamiltonian for time evolution
- `t_final::Real`: Final time for OTOC calculation
- `dt::Real`: Time step for OTOC sampling

# Returns
- `Vector{Float64}`: OTOC values for all sites and times (flattened)
"""
function run_otoc_exact(site1_list::Vector{Int}, site2::Int, N::Int, H::PauliSum, 
                        t_final::Real, dt::Real)
    # Check system size - warn if too large
    dim = 2^N
    if dim > 1024
        @warn "System size 2^$N = $dim is large. Exact method may be slow and memory-intensive."
    end
    
    O2_initial = magz(N, site2)
    all_otocs = Float64[]
    
    println("OTOC calculation using exact matrix exponentiation started")
    
    for t in 0.0:dt:t_final
        otoc_sites = Float64[]
        
        for site1 in site1_list
            O1 = magz(N, site1)
            otoc, _ = infinite_temp_otoc_exact(O1, O2_initial, H, t)
            push!(otoc_sites, real(otoc))
        end
        
        append!(all_otocs, otoc_sites)
        @printf("Finished calculation for t = %.4f\n", t)
    end
    
    println("OTOC calculation completed")
    return all_otocs
end

"""
    single_timestep_otoc_exact(site1::Int, site2::Int, N::Int, H::PauliSum, dt::Real)

Calculate OTOC for a single time step using exact matrix exponentiation.

# Arguments
- `site1::Int`: Site index for O1 operator
- `site2::Int`: Site index for O2 operator  
- `N::Int`: Number of qubits in the system
- `H::PauliSum`: Hamiltonian for time evolution
- `dt::Real`: Time step

# Returns
- `Float64`: The OTOC value at time dt
"""
function single_timestep_otoc_exact(site1::Int, site2::Int, N::Int, H::PauliSum, dt::Real)
    O1 = magz(N, site1)
    O2 = magz(N, site2)
    
    otoc, _ = infinite_temp_otoc_exact(O1, O2, H, dt)
    return real(otoc)
end
