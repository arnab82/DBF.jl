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

This function identifies pairs of Pauli operators in the input PauliSum that together
preserve particle number. A pair is kept if the sum of the two operators commutes with
the particle number operator. This searches all pairwise combinations to find pairs
that preserve particle number when combined.

This is useful for fermionic systems where individual terms may not preserve particle number
but pairs (e.g., from fermionic hopping terms c†ᵢcⱼ + c†ⱼcᵢ) do.

# Arguments
- `ps::PauliSum{N,T}`: A sum of Pauli operators to filter
- `N̂::PauliSum{N}`: The particle number operator

# Returns
- `PauliSum{N,T}`: A new PauliSum with only pairs that preserve particle number

# Example
```julia
N̂ = particle_number_operator(4)
# Create gradient from Hubbard model
ps_filtered = filter_particle_number_preserving_pairs(ps, N̂)  # keeps pairs that preserve PN
```
"""
function filter_particle_number_preserving_pairs(ps::PauliSum{N,T}, N̂::PauliSum{N}) where {N,T}
    result = PauliSum{N,T}()
    processed = Set{PauliBasis{N}}()
    
    # Convert to array for easier iteration
    ps_array = collect(ps)
    
    # First pass: check individual terms
    for (pauli, coeff) in ps_array
        # Skip if already processed
        pauli in processed && continue
        
        # Check if individual term preserves particle number
        p_op = Pauli(pauli)
        comm = N̂ * p_op - p_op * N̂
        coeff_clip!(comm, thresh=1e-12)
        
        if length(comm) == 0
            # Individual term preserves PN - keep it
            result[pauli] = coeff
            push!(processed, pauli)
        end
    end
    
    # Second pass: find pairs that preserve PN
    for i in 1:length(ps_array)
        p1, c1 = ps_array[i]
        
        # Skip if already processed
        p1 in processed && continue
        
        # Try to find a partner
        for j in (i+1):length(ps_array)
            p2, c2 = ps_array[j]
            
            # Skip if already processed
            p2 in processed && continue
            
            # Check if the pair sum preserves particle number
            pair_sum = c1 * Pauli(p1) + c2 * Pauli(p2)
            comm = N̂ * pair_sum - pair_sum * N̂
            coeff_clip!(comm, thresh=1e-12)
            
            if length(comm) == 0
                # Pair preserves particle number - keep both terms
                result[p1] = c1
                result[p2] = c2
                push!(processed, p1)
                push!(processed, p2)
                break  # Found a partner, move to next term
            end
        end
    end
    
    return result
end
