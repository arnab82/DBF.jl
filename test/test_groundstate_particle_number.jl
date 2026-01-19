using DBF
using PauliOperators
using LinearAlgebra
using Test
using Random

@testset "DBF Groundstate with Particle Number Preservation" begin
    @testset "Test particle number preservation in dbf_groundstate" begin
        N = 3
        Random.seed!(123)
        
        # Create a simple Hamiltonian
        H = DBF.heisenberg_1D(N, 1, 2, 3, z=0.1)
        DBF.coeff_clip!(H)
        H0 = deepcopy(H)
        
        # Find the lowest energy state
        kidx = argmin([real(expectation_value(H, Ket{N}(ψi))) for ψi in 1:2^N])
        ψ = Ket{N}(kidx)
        
        # Run DBF with particle number preservation
        res = DBF.dbf_groundstate(H, ψ,
                        max_iter=5, 
                        conv_thresh=1e-3,
                        evolve_coeff_thresh=1e-6,
                        grad_coeff_thresh=1e-6,
                        preserve_particle_number=true,
                        verbose=0)
        
        H_transformed = res["hamiltonian"]
        generators = res["generators"]
        angles = res["angles"]
        
        # Check that all generators preserve particle number
        for gen in generators
            @test DBF.preserves_particle_number(gen) == true
        end
        
        # Check that the energy is reasonable (should be close to original energy)
        e_original = real(expectation_value(H0, ψ))
        e_transformed = real(expectation_value(H_transformed, ψ))
        
        # The transformed Hamiltonian should have similar or lower energy
        @test e_transformed <= e_original + 1e-6
        
        println("Original energy: ", e_original)
        println("Transformed energy: ", e_transformed)
        println("Number of generators used: ", length(generators))
        println("All generators preserve particle number: ", all(DBF.preserves_particle_number.(generators)))
    end
    
    @testset "Compare with and without particle number preservation" begin
        N = 3
        Random.seed!(42)
        
        # Create a Hamiltonian
        H = DBF.heisenberg_1D(N, 1, 1, 1)
        DBF.coeff_clip!(H)
        
        # Find the lowest energy state
        kidx = argmin([real(expectation_value(H, Ket{N}(ψi))) for ψi in 1:2^N])
        ψ = Ket{N}(kidx)
        
        # Run without particle number preservation
        res_no_pn = DBF.dbf_groundstate(deepcopy(H), ψ,
                        max_iter=3,
                        conv_thresh=1e-3,
                        evolve_coeff_thresh=1e-6,
                        grad_coeff_thresh=1e-6,
                        preserve_particle_number=false,
                        verbose=0)
        
        # Run with particle number preservation
        res_with_pn = DBF.dbf_groundstate(deepcopy(H), ψ,
                        max_iter=3,
                        conv_thresh=1e-3,
                        evolve_coeff_thresh=1e-6,
                        grad_coeff_thresh=1e-6,
                        preserve_particle_number=true,
                        verbose=0)
        
        # Check that particle-number-preserving generators are indeed preserving
        for gen in res_with_pn["generators"]
            @test DBF.preserves_particle_number(gen) == true
        end
        
        # The non-preserving case may use some non-preserving generators
        # (not guaranteed, but likely)
        
        println("\nWithout particle number preservation:")
        println("  Number of generators: ", length(res_no_pn["generators"]))
        println("  Final energy: ", res_no_pn["energies"][end])
        
        println("\nWith particle number preservation:")
        println("  Number of generators: ", length(res_with_pn["generators"]))
        println("  Final energy: ", res_with_pn["energies"][end])
        println("  All generators preserve PN: ", all(DBF.preserves_particle_number.(res_with_pn["generators"])))
    end
end
