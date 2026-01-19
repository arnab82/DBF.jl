"""
    particle_number.jl

Functions for handling particle number preservation in Pauli operators.

In quantum computing with fermionic systems encoded in qubits, the particle number
operator counts the total number of particles (occupied states). For a system of N qubits,
the particle number operator is:

    N̂ = Σᵢ (I - Zᵢ)/2

where the sum runs over all qubits. This counts how many qubits are in the |1⟩ state.

A Pauli operator P preserves particle number if [N̂, P] = 0, i.e., if it commutes with
the particle number operator. For a single-qubit Pauli, this means:
- I and Z preserve particle number (they don't change occupancy)
- X and Y do not preserve particle number (they flip occupancy)

However, multi-qubit Paulis can preserve particle number even if they contain X or Y,
as long as the total number of occupancy changes sums to zero.
"""

"""
    particle_number_operator(N::Integer)

Create the particle number operator for N qubits.

The particle number operator is: N̂ = Σᵢ (I - Zᵢ)/2

This returns a `PauliSum` representing the number operator.

# Arguments
- `N::Integer`: Number of qubits

# Returns
- `PauliSum{N, ComplexF64}`: The particle number operator

# Example
```julia
N̂ = particle_number_operator(3)
```
"""
function particle_number_operator(N::Integer)
    result = PauliSum{N, ComplexF64}()
    
    # Add the identity term (N/2)
    identity_pauli = PauliBasis{N}(Int128(0), Int128(0))
    result[identity_pauli] = N/2
    
    # Subtract Zᵢ/2 for each qubit i
    for i in 1:N
        z_bits = Int128(1) << (i-1)
        z_pauli = PauliBasis{N}(z_bits, Int128(0))
        result[z_pauli] = get(result, z_pauli, 0.0) - 0.5
    end
    
    return result
end

"""
    preserves_particle_number(p::PauliBasis{N}) where N

Check if a Pauli operator preserves particle number.

A Pauli operator preserves particle number if it commutes with the particle number operator
N̂ = Σᵢ (I - Zᵢ)/2. This requires [Zᵢ, P] = 0 for all i.

Since Zᵢ anti-commutes with Xᵢ and Yᵢ, P preserves particle number if and only if it contains
NO X or Y operators on any qubit - it must be composed only of I and Z operators.

# Arguments
- `p::PauliBasis{N}`: A Pauli operator

# Returns
- `Bool`: true if the operator preserves particle number (contains only I and Z)

# Example
```julia
p1 = PauliBasis(Pauli(4, Z=[1,2]))  # preserves: only Z
p2 = PauliBasis(Pauli(4, X=[1]))    # does NOT preserve: has X
p3 = PauliBasis(Pauli(4, X=[1,2]))  # does NOT preserve: has X
```
"""
function preserves_particle_number(p::PauliBasis{N}) where N
    # In the symplectic representation, the x bit is set for both X and Y operators
    # P preserves particle number if and only if it has NO X or Y operators
    # i.e., p.x must be zero
    return p.x == 0
end

"""
    preserves_particle_number(p::Pauli{N}) where N

Check if a Pauli operator preserves particle number.

# Arguments
- `p::Pauli{N}`: A Pauli operator with coefficient

# Returns
- `Bool`: true if the operator preserves particle number
"""
preserves_particle_number(p::Pauli{N}) where N = preserves_particle_number(PauliBasis(p))

"""
    preserves_particle_number(ps::PauliSum{N,T}) where {N,T}

Check if all terms in a PauliSum preserve particle number.

# Arguments
- `ps::PauliSum{N,T}`: A sum of Pauli operators

# Returns
- `Bool`: true if all terms preserve particle number
"""
function preserves_particle_number(ps::PauliSum{N,T}) where {N,T}
    for (pauli, _) in ps
        if !preserves_particle_number(pauli)
            return false
        end
    end
    return true
end

"""
    filter_particle_number_preserving(ps::PauliSum{N,T}) where {N,T}

Filter a PauliSum to keep only terms that preserve particle number.

# Arguments
- `ps::PauliSum{N,T}`: A sum of Pauli operators

# Returns
- `PauliSum{N,T}`: A new PauliSum with only particle-number-preserving terms

# Example
```julia
ps = Pauli("X") + Pauli("Z") + Pauli("XX")
ps_filtered = filter_particle_number_preserving(ps)  # keeps only Z
```
"""
function filter_particle_number_preserving(ps::PauliSum{N,T}) where {N,T}
    result = PauliSum{N,T}()
    for (pauli, coeff) in ps
        if preserves_particle_number(pauli)
            result[pauli] = coeff
        end
    end
    return result
end

"""
    filter_particle_number_preserving(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}

Filter a PauliSum to keep only terms that commute with a given particle number operator.

This version checks actual commutation [N̂, p] = 0 rather than using a simple heuristic.
Use this when working in a transformed basis where the standard check may not apply.

# Arguments
- `ps::PauliSum{N,T}`: A sum of Pauli operators to filter
- `N̂::PauliSum{N}`: The particle number operator

# Returns
- `PauliSum{N,T}`: A new PauliSum with only terms that commute with N̂

# Example
```julia
N̂ = particle_number_operator(4)
ps = Pauli(4, X=[1]) + Pauli(4, Z=[1])
ps_filtered = filter_particle_number_preserving(ps, N̂)  # keeps only Z
```
"""
function filter_particle_number_preserving(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}
    result = PauliSum{N,T}()
    for (pauli, coeff) in ps
        # Check if [N̂, pauli] = 0
        p_op = Pauli(pauli)
        comm = N̂ * p_op - p_op * N̂
        coeff_clip!(comm, thresh=1e-12)
        
        if length(comm) == 0
            # Commutes with particle number operator
            result[pauli] = coeff
        end
    end
    return result
end
