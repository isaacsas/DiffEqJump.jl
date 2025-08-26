using JumpProcesses, DiffEqBase, Test, Statistics, Random
using StableRNGs
rng = StableRNG(12345)

@testset "SimpleSplitTauLeaping MassActionJump Tests" begin

    function run_split_ensemble(prob, n_traj, dt)
        sols = Vector{Any}(undef, n_traj)
        for i in 1:n_traj
            sols[i] = solve(prob, SimpleSplitTauLeaping(); dt)
        end
        return sols
    end
    
    # Test 1: Birth-death process via MassActionJumps
    @testset "Birth-Death MassActionJump" begin
        # ∅ → X (birth), X → ∅ (death)
        reactant_stoich = [Vector{Pair{Int,Int}}(),  # ∅ → X
                          [1 => 1]]                    # X → ∅
        net_stoich = [[1 => 1],                       # ∅ → X adds one
                      [1 => -1]]                       # X → ∅ removes one
        rates = [10.0, 0.1]  # birth rate, death rate
        
        ma_jumps = MassActionJump(rates, reactant_stoich, net_stoich)
        
        u0 = [100]
        tspan = (0.0, 100.0)
        dprob = DiscreteProblem(u0, tspan)
        
        # Run ensembles for statistics
        n_traj = 4000
        dt = .001

        # Direct method reference
        jprob_direct = JumpProblem(dprob, ma_jumps; rng, save_positions = (false, false))
        ensembleprob_direct = EnsembleProblem(jprob_direct)
        sol_direct = solve(ensembleprob_direct, SSAStepper(), 
                          EnsembleThreads(); trajectories=n_traj, saveat=tspan[2]/100)
        
        # SimpleSplitTauLeaping with small dt
        jprob_split = JumpProblem(dprob, PureLeaping(), ma_jumps; rng)
                
        sol_split = run_split_ensemble(jprob_split, n_traj, dt);
        
        # Extract final values
        direct_final = [sol.u[end][1] for sol in sol_direct]
        split_final = [sol.u[end][1] for sol in sol_split]
        
        # Test mean and variance (5% relative accuracy)
        @test mean(split_final) ≈ mean(direct_final) rtol=0.05
        @test var(split_final) ≈ var(direct_final) rtol=0.1
    end
    
    # Test 2: Simple reaction A + B → C
    @testset "A + B → C Reaction" begin
        reactant_stoich = [[1 => 1, 2 => 1]]  # A + B
        net_stoich = [[1 => -1, 2 => -1, 3 => 1]]  # -A -B +C
        rates = [0.001]
        
        ma_jumps = MassActionJump(rates, reactant_stoich, net_stoich)
        
        u0 = [100, 100, 0]
        tspan = (0.0, 10.0)
        dprob = DiscreteProblem(u0, tspan)
        
        n_traj = 1000
        
        # Direct reference
        jprob_direct = JumpProblem(dprob, ma_jumps; rng, save_positions = (false, false))
        ensembleprob_direct = EnsembleProblem(jprob_direct)
        sol_direct = solve(ensembleprob_direct, SSAStepper(),
                          EnsembleThreads(); trajectories=n_traj, saveat=tspan[2]/100)
        
        # SimpleSplitTauLeaping
        jprob_split = JumpProblem(dprob, PureLeaping(), ma_jumps; rng)
        sol_split = run_split_ensemble(jprob_split, n_traj, 0.0001)
        
        # Compare means of all species at final time
        for species in 1:3
            direct_vals = [sol.u[end][species] for sol in sol_direct]
            split_vals = [sol.u[end][species] for sol in sol_split]
            
            @test mean(split_vals) ≈ mean(direct_vals) rtol=0.05
            @test var(split_vals) ≈ var(direct_vals) rtol=0.1
        end
        
        # Check conservation: A + B → C means C = A(0) - A = B(0) - B
        for sol in sol_split
            @test sol.u[end][3] == u0[1] - sol.u[end][1]  # C = A(0) - A
            @test sol.u[end][3] == u0[2] - sol.u[end][2]  # C = B(0) - B
        end
    end
    
    # Test 3: Lotka-Volterra predator-prey
    @testset "Lotka-Volterra System" begin
        # X → 2X (prey birth)
        # X + Y → 2Y (predation) 
        # Y → ∅ (predator death)
        reactant_stoich = [[1 => 1],           # X
                          [1 => 1, 2 => 1],     # X + Y
                          [2 => 1]]             # Y
        net_stoich = [[1 => 1],                # X births
                      [1 => -1, 2 => 1],        # X dies, Y births
                      [2 => -1]]                # Y dies
        rates = [1.0, 0.001, .1]
        
        ma_jumps = MassActionJump(rates, reactant_stoich, net_stoich)
        
        u0 = [100, 100]  # Initial prey and predator
        tspan = (0.0, 20.0)
        dprob = DiscreteProblem(u0, tspan)
        
        n_traj = 100
        saveat = [5.0, 10.0, 15.0, 20.0]
        
        # Direct
        jprob_direct = JumpProblem(dprob, ma_jumps; rng, save_positions = (false, false))
        ensembleprob_direct = EnsembleProblem(jprob_direct)
        sol_direct = solve(ensembleprob_direct, SSAStepper(); trajectories=n_traj, saveat)
        
        # SimpleSplitTauLeaping with very small dt for accuracy
        jprob_split = JumpProblem(dprob, PureLeaping(), ma_jumps)
        sol_split = run_split_ensemble(jprob_split, n_traj, 0.0001)
        
        # Sample at multiple time points
        for test_t in saveat
            # Find closest time index
            t_idx_direct = findfirst(t -> t >= test_t, sol_direct.u[1].t)
            t_idx_split = findfirst(t -> t >= test_t, sol_split[1].t)
            
            for species in 1:2
                direct_vals = [sol(test_t; idxs = species) for sol in sol_direct]
                split_vals = [sol.u[t_idx_split][species] for sol in sol_split]
                @test mean(split_vals) ≈ mean(direct_vals) rtol=0.05
            end
        end
    end

    # Test 4: Validation tests
    @testset "Validation Tests" begin
        u0 = [10]
        tspan = (0.0, 10.0)
        dprob = DiscreteProblem(u0, tspan)
        
        # Test with non-PureLeaping aggregator
        reactant_stoich = [[1 => 1]]
        net_stoich = [[1 => -1]]
        rates = [0.1]
        ma_jumps = MassActionJump(rates, reactant_stoich, net_stoich)
        
        jprob_wrong = JumpProblem(dprob, Direct(), ma_jumps)
        @test_nowarn solve(jprob_wrong, SimpleSplitTauLeaping(), dt=0.1)
        
        # Test with no MassActionJumps (should fail)
        rj = RegularJump((out, u, p, t) -> out[1] = u[1], (du, u, p, t, counts, mark) -> du[1] = counts[1], 1)
        jprob_no_maj = JumpProblem(dprob, PureLeaping(), rj)
        @test_throws ErrorException solve(jprob_no_maj, SimpleSplitTauLeaping(), dt=0.1)
        
        # Test with ConstantRateJumps (should fail)
        crj = ConstantRateJump((u, p, t) -> u[1], integrator -> integrator.u[1] -= 1)
        jprob_crj = JumpProblem(dprob, PureLeaping(), crj)
        @test_throws ErrorException solve(jprob_crj, SimpleSplitTauLeaping(), dt=0.1)
    end

    # Test 5: Parameter requirements
    @testset "Parameter Requirements" begin
        u0 = [10]
        tspan = (0.0, 10.0)
        dprob = DiscreteProblem(u0, tspan)
        
        reactant_stoich = [[1 => 1]]
        net_stoich = [[1 => -1]]
        rates = [0.1]
        ma_jumps = MassActionJump(rates, reactant_stoich, net_stoich)
        jprob = JumpProblem(dprob, PureLeaping(), ma_jumps)
        
        # Test dt is required
        @test_throws ErrorException solve(jprob, SimpleSplitTauLeaping())
        
        # Test with dt provided
        @test_nowarn solve(jprob, SimpleSplitTauLeaping(), dt=0.1)
        
        # Test with seed
        sol1 = solve(jprob, SimpleSplitTauLeaping(), dt=0.1, seed=12345)
        sol2 = solve(jprob, SimpleSplitTauLeaping(), dt=0.1, seed=12345)
        @test sol1.u == sol2.u
    end
end