# OTOC Calculation using Pauli Propagation

This directory contains tools for calculating Out-of-Time-Ordered Correlators (OTOCs) using Pauli operator propagation.

## Files

- `src/observables.jl`: Observable operators (magnetization in X, Y, Z directions)
- `src/otoc_pauli.jl`: OTOC calculation functions using Pauli propagation
- `otoc_comparison.jl`: Main script for running OTOC calculations and analysis

## Usage

Run the OTOC comparison script with:

```bash
julia --project=. otoc_comparison.jl [N] [t_final] [dt] [Jx] [Jy] [Jz]
```

### Arguments

- `N` - Number of qubits (default: 6)
- `t_final` - Final time for evolution (default: 1.0)
- `dt` - Time step for OTOC sampling (default: 0.1)
- `Jx` - X coupling strength (default: 1.0)
- `Jy` - Y coupling strength (default: 1.0)
- `Jz` - Z coupling strength (default: 1.0)

### Example

```bash
# Run with 4 qubits, evolve to t=0.5 with dt=0.1
julia --project=. otoc_comparison.jl 4 0.5 0.1

# Run with default parameters (6 qubits, t=1.0, dt=0.1)
julia --project=. otoc_comparison.jl
```

## Theory

### Out-of-Time-Ordered Correlator (OTOC)

The OTOC is a measure of quantum information scrambling and chaos in quantum systems. At infinite temperature, it is defined as:

```
F(t) = ⟨[W(t), V(0)]†[W(t), V(0)]⟩ / (⟨W†W⟩⟨V†V⟩)
```

where:
- `W(t) = exp(iHt) W exp(-iHt)` is the time-evolved operator
- `V(0)` is the initial operator
- `[A, B] = AB - BA` is the commutator

### Implementation

The OTOC calculation uses:

1. **Pauli Propagation**: Time evolution is performed using the Pauli operator evolution formula:
   ```
   O(θ) = exp(iθ/2 G) O exp(-iθ/2 G)
   ```
   where for non-commuting operators `[G, O] ≠ 0`:
   ```
   O(θ) = O cos(θ) - i sin(θ) G*O
   ```

2. **Trotterization**: The Hamiltonian evolution is approximated by breaking it into small time steps.

3. **Infinite Temperature**: The calculation assumes infinite temperature, simplifying the thermal averaging.

## Output

The script outputs:

1. **OTOC Evolution**: Values of OTOC for all sites at different times
2. **Single Time-Step**: OTOC value for one specific time step
3. **Analysis**: Mean and standard deviation of OTOC across sites at each time

### Example Output

```
============================================================
OTOC Results
============================================================
Time    Site 1  Site 2  Site 3  Site 4
------------------------------------------------------------
0.00    0.000000    0.000000    0.000000    0.000000    
0.10    0.004993    0.009982    0.004994    0.000001    
0.20    0.019892    0.039709    0.019895    0.000023    
0.30    0.044456    0.088535    0.044464    0.000119    
0.40    0.078292    0.155398    0.078311    0.000384    
0.50    0.120863    0.238854    0.120900    0.000954    

============================================================
OTOC Analysis
============================================================
OTOC spread over sites at different times:
  t = 0.00: mean = 0.000000, std = 0.000000
  t = 0.10: mean = 0.004992, std = 0.003529
  t = 0.20: mean = 0.019880, std = 0.014031
  t = 0.30: mean = 0.044394, std = 0.031260
  t = 0.40: mean = 0.078096, std = 0.054806
  t = 0.50: mean = 0.120393, std = 0.084112
```

## Functions

### `infinite_temp_otoc_pauli(O1, O2, H, t; dt=0.01)`

Calculate OTOC at infinite temperature for time `t`.

**Parameters:**
- `O1`: First observable (PauliSum)
- `O2`: Second observable to be time-evolved (PauliSum)
- `H`: Hamiltonian (PauliSum)
- `t`: Evolution time
- `dt`: Time step for Trotterization

**Returns:**
- OTOC value (Float64)
- Time-evolved O2(t) (PauliSum)

### `run_otoc_pauli(site1_list, site2, N, H, t_final, dt; evolution_dt=0.01)`

Run OTOC calculations for multiple sites over time.

**Parameters:**
- `site1_list`: List of site indices for O1 operators
- `site2`: Site index for O2 operator
- `N`: Number of qubits
- `H`: Hamiltonian
- `t_final`: Final time
- `dt`: Time step for sampling
- `evolution_dt`: Time step for evolution

**Returns:**
- Vector of OTOC values (flattened over sites and times)

### `single_timestep_otoc_pauli(site1, site2, N, H, dt)`

Calculate OTOC for a single time step.

**Parameters:**
- `site1`, `site2`: Site indices
- `N`: Number of qubits
- `H`: Hamiltonian
- `dt`: Time step

**Returns:**
- OTOC value (Float64)

## Physical Interpretation

- **OTOC ≈ 0**: Operators commute; information is localized
- **OTOC > 0**: Operators don't commute; information has spread
- **OTOC growth**: Indicates information scrambling and quantum chaos
- **Spread across sites**: Shows spatial propagation of quantum information

## Notes

- The implementation uses the existing DBF.jl framework for Pauli operator manipulation
- Time evolution uses Trotterization, so accuracy depends on the time step `dt`
- Smaller `dt` values give more accurate results but require more computation
