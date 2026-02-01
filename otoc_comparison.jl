#!/usr/bin/env julia

using Printf
using Dates
using PauliOperators
using LinearAlgebra
using DBF

include("./src/observables.jl")
include("./src/otoc_pauli.jl")

"""
Main script for comparing OTOC calculations using Pauli propagation.

This script demonstrates:
1. OTOC calculation using Pauli propagation for multiple sites over time
2. Single time-step OTOC comparison

Usage: julia otoc_comparison.jl [N] [t_final] [dt] [Jx] [Jy] [Jz]

Arguments:
    N       - Number of qubits (default: 6)
    t_final - Final time for evolution (default: 1.0)
    dt      - Time step for OTOC sampling (default: 0.1)
    Jx      - X coupling strength (default: 1.0)
    Jy      - Y coupling strength (default: 1.0)  
    Jz      - Z coupling strength (default: 1.0)
"""
function main()
    start_time = time()
    now = Dates.now()
    formatted = Dates.format(now, "yyyy-mm-dd HH:MM:SS")
    println("START : $formatted\n")

    # Parse command line arguments
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
    t_final = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.0
    dt = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.1
    Jx = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    Jy = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0
    Jz = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 1.0

    println("="^60)
    println("OTOC Calculation using Pauli Propagation")
    println("="^60)
    println("Parameters:")
    println("  Number of qubits (N):  $N")
    println("  Final time (t_final):  $t_final")
    println("  Time step (dt):        $dt")
    println("  Coupling (Jx, Jy, Jz): ($Jx, $Jy, $Jz)")
    println("="^60)
    println()

    # Create Heisenberg Hamiltonian
    println("Creating Heisenberg Hamiltonian...")
    H = DBF.af_heisenberg(N, Jx, Jy, Jz)
    println("Hamiltonian has $(length(H)) terms")
    println()

    # Define sites for OTOC calculation
    site2 = N ÷ 2  # Middle site for O2
    site1_list = collect(1:N)  # All sites for O1
    
    println("OTOC Configuration:")
    println("  O2 site: $site2")
    println("  O1 sites: $(site1_list)")
    println()

    # Calculate OTOC over time for all sites
    println("="^60)
    println("Computing OTOC using Pauli Propagation")
    println("="^60)
    otoc_results = run_otoc_pauli(site1_list, site2, N, H, t_final, dt; evolution_dt=dt/10)
    
    # Reshape results for analysis
    n_times = Int(round(t_final / dt)) + 1
    otoc_matrix = reshape(otoc_results, (N, n_times))
    
    println()
    println("="^60)
    println("OTOC Results")
    println("="^60)
    println("Time\t", join(["Site $i" for i in 1:N], "\t"))
    println("-"^60)
    
    for (i, t) in enumerate(0.0:dt:t_final)
        @printf("%.2f\t", t)
        for site in 1:N
            @printf("%.6f\t", otoc_matrix[site, i])
        end
        println()
    end
    println()

    # Single time-step comparison
    println("="^60)
    println("Single Time-Step OTOC Comparison")
    println("="^60)
    single_dt = dt
    site1_test = 1
    site2_test = site2
    
    println("Calculating OTOC for single time step:")
    println("  Site 1: $site1_test")
    println("  Site 2: $site2_test")
    println("  Time step: $single_dt")
    
    otoc_single = single_timestep_otoc_pauli(site1_test, site2_test, N, H, single_dt)
    @printf("  OTOC value: %.8f\n", otoc_single)
    println()

    # Calculate spread of OTOC (measure of information scrambling)
    println("="^60)
    println("OTOC Analysis")
    println("="^60)
    println("OTOC spread over sites at different times:")
    for (i, t) in enumerate(0.0:dt:t_final)
        otoc_values = otoc_matrix[:, i]
        avg_otoc = sum(otoc_values) / length(otoc_values)
        std_otoc = sqrt(sum((otoc_values .- avg_otoc).^2) / length(otoc_values))
        @printf("  t = %.2f: mean = %.6f, std = %.6f\n", t, avg_otoc, std_otoc)
    end
    println()

    # Print timing information
    end_time = time()
    elapsed = end_time - start_time
    println("="^60)
    @printf("Total computation time: %.2f seconds\n", elapsed)
    now = Dates.now()
    formatted = Dates.format(now, "yyyy-mm-dd HH:MM:SS")
    println("END : $formatted")
    println("="^60)
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
