# Observables for OTOC calculations
# Assumes PauliOperators is already loaded

"""
    magz(N::Int, site::Int)

Create a magnetization operator Sz for the specified site in an N-qubit system.
Returns a PauliSum representing the Z operator at the given site.

# Arguments
- `N::Int`: Total number of qubits in the system
- `site::Int`: Site index (1-based) where the Z operator is applied

# Returns
- `PauliSum{N}`: The magnetization operator
"""
function magz(N::Int, site::Int)
    return PauliSum(Pauli(N, Z=[site]))
end

"""
    magx(N::Int, site::Int)

Create a magnetization operator Sx for the specified site in an N-qubit system.
Returns a PauliSum representing the X operator at the given site.

# Arguments
- `N::Int`: Total number of qubits in the system
- `site::Int`: Site index (1-based) where the X operator is applied

# Returns
- `PauliSum{N}`: The magnetization operator
"""
function magx(N::Int, site::Int)
    return PauliSum(Pauli(N, X=[site]))
end

"""
    magy(N::Int, site::Int)

Create a magnetization operator Sy for the specified site in an N-qubit system.
Returns a PauliSum representing the Y operator at the given site.

# Arguments
- `N::Int`: Total number of qubits in the system
- `site::Int`: Site index (1-based) where the Y operator is applied

# Returns
- `PauliSum{N}`: The magnetization operator
"""
function magy(N::Int, site::Int)
    return PauliSum(Pauli(N, Y=[site]))
end
