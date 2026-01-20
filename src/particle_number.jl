"""
    particle_number.jl

Functions for particle number preservation in fermionic systems encoded as Pauli operators.

The particle number operator for N qubits is:

    N̂ = Σᵢ (I - Zᵢ)/2

This counts the number of qubits in the |1⟩ state (occupied orbitals).

## Key Concept

Individual Pauli operators rarely preserve particle number in Jordan-Wigner encoding.
However, pairs of operators (e.g., from fermionic hopping terms c†ᵢcⱼ + c†ⱼcᵢ) can
preserve particle number together.

## Recommended Usage

For fermionic systems, use **pair-based filtering**:
- `filter_particle_number_preserving_pairs(ps, N̂)` - finds pairs that preserve PN
- Enable with `use_pair_filter=true` in dbf_groundstate
- Also enable propagation-time filtering with `preserve_particle_number=true`

## Individual Term Filtering (Limited Use)

- `preserves_particle_number(p)` - checks if operator has no X/Y (strict check)
- `filter_particle_number_preserving(ps, N̂)` - filters individual terms (very restrictive)
- These are mainly useful for testing, not practical optimization
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
    filter_particle_number_preserving(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}

Filter a PauliSum to keep only terms that commute with a given particle number operator.

This version checks actual commutation [N̂, p] = 0 rather than using a simple heuristic.
Use this when working in a transformed basis where the standard check may not apply.

Note: For fermionic systems, individual terms often don't preserve particle number.
Use `filter_particle_number_preserving_pairs` for better results.

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

"""
    filter_particle_number_preserving_pairs(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}

Filter a PauliSum using a pair-based approach for particle number preservation.

For fermionic Hamiltonians that preserve particle number, the commutator [P, H] naturally
produces pairs of operators that together preserve particle number. This function uses a
simple heuristic: keep all terms with even X+Y count (which form pairs).

This is much faster than computing commutators and works well for fermionic systems where
individual terms may not preserve PN but pairs (e.g., from fermionic hopping c†ᵢcⱼ + c†ⱼcᵢ) do.

# Arguments
- `ps::PauliSum{N,T}`: A sum of Pauli operators to filter
- `N̂::PauliSum{N}`: The particle number operator (not used, kept for API compatibility)

# Returns
- `PauliSum{N,T}`: A new PauliSum with only particle-number-preserving terms

# Example
```julia
N̂ = particle_number_operator(4)
ps_filtered = filter_particle_number_preserving_pairs(ps, N̂)
```
"""
function filter_particle_number_preserving_pairs(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}
    result = PauliSum{N,T}()
    
    # Fast heuristic: keep terms with even X+Y count (including 0)
    # This works because:
    # 1. Terms with X/Y=0 (only I/Z) preserve PN individually
    # 2. Terms with even X/Y typically form pairs that preserve PN
    # 3. Much faster than computing commutators
    
    for (pauli, coeff) in ps
        xy_count = count_ones(pauli.x)
        if iseven(xy_count)
            result[pauli] = coeff
        end
    end
    
    return result
end
