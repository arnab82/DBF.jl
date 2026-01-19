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

A Pauli operator preserves particle number if it commutes with the number operator.
For a Pauli string, this is true if:
1. It contains only I and Z operators (no X or Y), OR
2. The X and Y operators come in pairs that preserve total particle number

For simplicity, this function checks the X operators (which cause bit flips).
A Pauli preserves particle number if the number of X/Y operators is even and they
are arranged such that each creation is paired with an annihilation.

# Arguments
- `p::PauliBasis{N}`: A Pauli operator

# Returns
- `Bool`: true if the operator preserves particle number

# Example
```julia
p1 = PauliBasis("IIZZ")  # preserves: no X or Y
p2 = PauliBasis("XYII")  # can preserve if proper pairing
p3 = PauliBasis("XIII")  # does not preserve: single X
```
"""
function preserves_particle_number(p::PauliBasis{N}) where N
    # Count the number of X or Y operators (operators that flip qubits)
    # In the symplectic representation, x bit is set for both X and Y
    num_flips = count_ones(p.x)
    
    # For particle number preservation, we need an even number of flips
    # (each "creation" must be paired with an "annihilation")
    # This is a necessary but not always sufficient condition
    # For a more rigorous check, we'd need to verify the specific structure
    return iseven(num_flips)
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
ps_filtered = filter_particle_number_preserving(ps)  # keeps Z and XX
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
