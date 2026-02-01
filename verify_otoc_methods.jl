#!/usr/bin/env julia

"""
OTOC Method Comparison Script

This script compares two methods for calculating Out-of-Time-Ordered Correlators (OTOCs):
1. Trotterization-based method (approximate, fast)
2. Exact matrix exponentiation method (exact, slower)

The script verifies that both methods produce the same results within numerical precision.
"""

using Printf
using Dates
using PauliOperators
using LinearAlgebra
using DBF

include("./src/observables.jl")
include("./src/otoc_pauli.jl")
include("./src/otoc_exact.jl")

"""
Compare OTOC values from two methods and check if they agree within tolerance.
"""
function compare_otoc_values(otoc1::Vector{Float64}, otoc2::Vector{Float64}, 
                            method1_name::String, method2_name::String;
                            rtol::Float64=1e-6, atol::Float64=1e-8)
    if length(otoc1) != length(otoc2)
        error("OTOC vectors have different lengths: $(length(otoc1)) vs $(length(otoc2))")
    end
    
    max_diff = 0.0
    max_rel_diff = 0.0
    num_points = length(otoc1)
    num_agree = 0
    
    println("\n" * "="^70)
    println("Comparing $method1_name vs $method2_name")
    println("="^70)
    println("Number of comparison points: $num_points")
    println("Absolute tolerance: $atol")
    println("Relative tolerance: $rtol")
    println()
    
    for i in 1:num_points
        diff = abs(otoc1[i] - otoc2[i])
        rel_diff = abs(otoc1[i]) > atol ? diff / abs(otoc1[i]) : 0.0
        
        max_diff = max(max_diff, diff)
        max_rel_diff = max(max_rel_diff, rel_diff)
        
        # Check if values agree within tolerance
        if diff < atol || rel_diff < rtol
            num_agree += 1
        else
            @printf("  Point %3d: %s = %.8e, %s = %.8e, diff = %.3e, rel_diff = %.3e ❌\n",
                    i, method1_name, otoc1[i], method2_name, otoc2[i], diff, rel_diff)
        end
    end
    
    agreement_pct = 100.0 * num_agree / num_points
    
    println()
    println("Results:")
    println("  Points in agreement: $num_agree / $num_points ($(@sprintf("%.1f", agreement_pct))%)")
    @printf("  Maximum absolute difference: %.3e\n", max_diff)
    @printf("  Maximum relative difference: %.3e\n", max_rel_diff)
    println()
    
    if num_agree == num_points
        println("✓ SUCCESS: All values agree within tolerance!")
        return true
    else
        println("✗ FAILURE: Some values differ beyond tolerance")
        return false
    end
end

function main()
    start_time = time()
    now = Dates.now()
    formatted = Dates.format(now, "yyyy-MM-dd HH:mm:ss")
    println("START : $formatted\n")

    # Parse command line arguments
    N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4
    t_final = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.3
    dt = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.1
    Jx = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    Jy = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0
    Jz = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 1.0
    evolution_dt = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 0.01

    println("="^70)
    println("OTOC Method Comparison")
    println("="^70)
    println("Parameters:")
    println("  Number of qubits (N):     $N")
    println("  Final time (t_final):     $t_final")
    println("  Time step (dt):           $dt")
    println("  Coupling (Jx, Jy, Jz):    ($Jx, $Jy, $Jz)")
    println("  Evolution dt (Trotter):   $evolution_dt")
    println("="^70)
    println()

    # Warn if system is large
    if N > 6
        println("⚠ Warning: System size N=$N may be too large for exact method")
        println("  Matrix dimension: $(2^N)")
        println("  Consider using N ≤ 6 for comparison")
        println()
    end

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

    # Method 1: Trotterization-based evolution
    println("="^70)
    println("Method 1: Trotterization-based Evolution")
    println("="^70)
    t1_start = time()
    otoc_trotter = run_otoc_pauli(site1_list, site2, N, H, t_final, dt; evolution_dt=evolution_dt)
    t1_elapsed = time() - t1_start
    @printf("Method 1 completed in %.3f seconds\n\n", t1_elapsed)

    # Method 2: Exact matrix exponentiation
    println("="^70)
    println("Method 2: Exact Matrix Exponentiation")
    println("="^70)
    t2_start = time()
    otoc_exact = run_otoc_exact(site1_list, site2, N, H, t_final, dt)
    t2_elapsed = time() - t2_start
    @printf("Method 2 completed in %.3f seconds\n\n", t2_elapsed)

    # Compare results
    success = compare_otoc_values(otoc_trotter, otoc_exact, 
                                  "Trotter", "Exact";
                                  rtol=1e-4, atol=1e-6)

    # Display some example values
    n_times = Int(round(t_final / dt)) + 1
    otoc_trotter_matrix = reshape(otoc_trotter, (N, n_times))
    otoc_exact_matrix = reshape(otoc_exact, (N, n_times))
    
    println("\n" * "="^70)
    println("Sample OTOC Values (Site 1, different times)")
    println("="^70)
    println("Time\t\tTrotter\t\tExact\t\tDifference")
    println("-"^70)
    for (i, t) in enumerate(0.0:dt:t_final)
        val_trotter = otoc_trotter_matrix[1, i]
        val_exact = otoc_exact_matrix[1, i]
        diff = abs(val_trotter - val_exact)
        @printf("%.2f\t\t%.8f\t%.8f\t%.3e\n", t, val_trotter, val_exact, diff)
    end
    println()

    # Performance comparison
    println("="^70)
    println("Performance Comparison")
    println("="^70)
    @printf("Trotter method:  %.3f seconds (%.1fx)\n", t1_elapsed, 1.0)
    @printf("Exact method:    %.3f seconds (%.1fx)\n", t2_elapsed, t2_elapsed/t1_elapsed)
    println("="^70)
    println()

    # Print timing information
    end_time = time()
    elapsed = end_time - start_time
    println("="^70)
    @printf("Total computation time: %.2f seconds\n", elapsed)
    now = Dates.now()
    formatted = Dates.format(now, "yyyy-MM-dd HH:mm:ss")
    println("END : $formatted")
    println("="^70)
    
    # Exit with appropriate status
    if success
        println("\n✓ Validation PASSED: Both methods produce consistent results")
        exit(0)
    else
        println("\n✗ Validation FAILED: Methods produce different results")
        exit(1)
    end
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
