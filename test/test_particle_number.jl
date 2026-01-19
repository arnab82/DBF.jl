using DBF
using PauliOperators
using LinearAlgebra
using Test

@testset "Particle Number Preservation" begin
    @testset "particle_number_operator" begin
        # Test for 2 qubits
        N = 2
        N̂ = DBF.particle_number_operator(N)
        
        # Check that we have the right number of terms
        @test length(N̂) == N + 1  # Identity + N Z terms
        
        # Check identity term coefficient
        identity_pauli = PauliBasis{N}(Int128(0), Int128(0))
        @test N̂[identity_pauli] ≈ N/2
        
        # Check Z terms
        for i in 1:N
            z_bits = Int128(1) << (i-1)
            z_pauli = PauliBasis{N}(z_bits, Int128(0))
            @test N̂[z_pauli] ≈ -0.5
        end
    end
    
    @testset "preserves_particle_number - PauliBasis" begin
        N = 4
        
        # Identity preserves particle number (no X or Y)
        p_I = PauliBasis{N}(Int128(0), Int128(0))
        @test DBF.preserves_particle_number(p_I) == true
        
        # Single Z preserves particle number (no X or Y)
        p_Z = PauliBasis(Pauli(N, Z=[1]))
        @test DBF.preserves_particle_number(p_Z) == true
        
        # Multiple Z preserves particle number (no X or Y)
        p_ZZ = PauliBasis(Pauli(N, Z=[1, 2]))
        @test DBF.preserves_particle_number(p_ZZ) == true
        
        # Single X does NOT preserve particle number
        p_X = PauliBasis(Pauli(N, X=[1]))
        @test DBF.preserves_particle_number(p_X) == false
        
        # Single Y does NOT preserve particle number
        p_Y = PauliBasis(Pauli(N, Y=[1]))
        @test DBF.preserves_particle_number(p_Y) == false
        
        # XX does NOT preserve particle number (has X operators)
        p_XX = PauliBasis(Pauli(N, X=[1, 2]))
        @test DBF.preserves_particle_number(p_XX) == false
        
        # XY does NOT preserve particle number (has X and Y operators)
        p_XY = PauliBasis(Pauli(N, X=[1], Y=[2]))
        @test DBF.preserves_particle_number(p_XY) == false
        
        # XXX does NOT preserve particle number (has X operators)
        p_XXX = PauliBasis(Pauli(N, X=[1, 2, 3]))
        @test DBF.preserves_particle_number(p_XXX) == false
        
        # XXXX does NOT preserve particle number (has X operators)
        p_XXXX = PauliBasis(Pauli(N, X=[1, 2, 3, 4]))
        @test DBF.preserves_particle_number(p_XXXX) == false
        
        # XZ does NOT preserve particle number (has X operator)
        p_XZ = PauliBasis(Pauli(N, X=[1], Z=[2]))
        @test DBF.preserves_particle_number(p_XZ) == false
    end
    
    @testset "preserves_particle_number - Pauli" begin
        N = 4
        
        # Test with Pauli type (with coefficient)
        p_Z = Pauli(N, Z=[1])
        @test DBF.preserves_particle_number(p_Z) == true
        
        p_X = Pauli(N, X=[1])
        @test DBF.preserves_particle_number(p_X) == false
        
        # XX does NOT preserve (has X operators)
        p_XX = Pauli(N, X=[1, 2])
        @test DBF.preserves_particle_number(p_XX) == false
    end
    
    @testset "preserves_particle_number - PauliSum" begin
        N = 4
        
        # PauliSum with only particle-number-preserving terms (only I and Z)
        ps1 = PauliSum(N)
        ps1 += Pauli(N)  # Identity - preserves
        ps1 += Pauli(N, Z=[1])  # Z - preserves
        ps1 += Pauli(N, Z=[1, 2])  # ZZ - preserves
        @test DBF.preserves_particle_number(ps1) == true
        
        # PauliSum with a non-preserving term (has X)
        ps2 = PauliSum(N)
        ps2 += Pauli(N, Z=[1])
        ps2 += Pauli(N, X=[1])  # Single X - does NOT preserve
        @test DBF.preserves_particle_number(ps2) == false
    end
    
    @testset "filter_particle_number_preserving" begin
        N = 4
        
        # Create a PauliSum with mixed terms
        ps = PauliSum(N)
        ps += Pauli(N)  # Identity - preserves (no X/Y)
        ps += Pauli(N, Z=[1])  # Z - preserves (no X/Y)
        ps += Pauli(N, Z=[1, 2])  # ZZ - preserves (no X/Y)
        ps += Pauli(N, X=[1])  # X - does NOT preserve
        ps += Pauli(N, X=[1, 2])  # XX - does NOT preserve
        ps += Pauli(N, Y=[1])  # Y - does NOT preserve
        ps += Pauli(N, X=[1, 2, 3])  # XXX - does NOT preserve
        ps += Pauli(N, X=[1, 2, 3, 4])  # XXXX - does NOT preserve
        
        # Filter to keep only preserving terms
        ps_filtered = DBF.filter_particle_number_preserving(ps)
        
        # Should have 3 terms: Identity, Z, ZZ (only those with no X/Y)
        @test length(ps_filtered) == 3
        
        # Check that all remaining terms preserve particle number
        @test DBF.preserves_particle_number(ps_filtered) == true
        
        # Check that the filtered terms are correct
        @test haskey(ps_filtered, PauliBasis{N}(Int128(0), Int128(0)))  # Identity
        @test haskey(ps_filtered, PauliBasis(Pauli(N, Z=[1])))  # Z
        @test haskey(ps_filtered, PauliBasis(Pauli(N, Z=[1, 2])))  # ZZ
        
        # Check that non-preserving terms are not present
        @test !haskey(ps_filtered, PauliBasis(Pauli(N, X=[1])))  # X
        @test !haskey(ps_filtered, PauliBasis(Pauli(N, Y=[1])))  # Y
        @test !haskey(ps_filtered, PauliBasis(Pauli(N, X=[1, 2, 3])))  # XXX
    end
end
