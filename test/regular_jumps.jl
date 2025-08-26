using JumpProcesses, DiffEqBase
using Test, LinearAlgebra
using StableRNGs
rng = StableRNG(12345)

function regular_rate(out, u, p, t)
    out[1] = (0.1 / 1000.0) * u[1] * u[2]
    out[2] = 0.01u[2]
end

const dc = zeros(3, 2)
dc[1, 1] = -1
dc[2, 1] = 1
dc[2, 2] = -1
dc[3, 2] = 1

function regular_c(du, u, p, t, counts, mark)
    mul!(du, dc, counts)
end

rj = RegularJump(regular_rate, regular_c, 2)
jumps = JumpSet(rj)
prob = DiscreteProblem([999, 1, 0], (0.0, 250.0))
jump_prob = JumpProblem(prob, PureLeaping(), rj; rng)
sol = solve(jump_prob, SimpleTauLeaping(); dt = 1.0)

# Test PureLeaping aggregator functionality
@testset "PureLeaping Aggregator Tests" begin
    # Test with MassActionJump
    u0 = [10, 5, 0]
    tspan = (0.0, 10.0)
    p = [0.1, 0.2]
    prob = DiscreteProblem(u0, tspan, p)
    
    # Create MassActionJump
    reactant_stoich = [[1 => 1], [1 => 2]]
    net_stoich = [[1 => -1, 2 => 1], [1 => -2, 3 => 1]]
    rates = [0.1, 0.05]
    maj = MassActionJump(rates, reactant_stoich, net_stoich)
    
    # Test PureLeaping JumpProblem creation
    jp_pure = JumpProblem(prob, PureLeaping(), JumpSet(maj))
    @test jp_pure.aggregator isa PureLeaping
    @test jp_pure.discrete_jump_aggregation === nothing
    @test jp_pure.massaction_jump !== nothing
    @test length(jp_pure.jump_callback.discrete_callbacks) == 0
    
    # Test with ConstantRateJump
    rate(u, p, t) = p[1] * u[1]
    affect!(integrator) = (integrator.u[1] -= 1; integrator.u[3] += 1)
    crj = ConstantRateJump(rate, affect!)
    
    jp_pure_crj = JumpProblem(prob, PureLeaping(), JumpSet(crj))
    @test jp_pure_crj.aggregator isa PureLeaping
    @test jp_pure_crj.discrete_jump_aggregation === nothing
    @test length(jp_pure_crj.constant_jumps) == 1
    
    # Test with VariableRateJump
    vrate(u, p, t) = t * p[1] * u[1]
    vaffect!(integrator) = (integrator.u[1] -= 1; integrator.u[3] += 1)
    vrj = VariableRateJump(vrate, vaffect!)
    
    jp_pure_vrj = JumpProblem(prob, PureLeaping(), JumpSet(vrj))
    @test jp_pure_vrj.aggregator isa PureLeaping
    @test jp_pure_vrj.discrete_jump_aggregation === nothing
    @test length(jp_pure_vrj.variable_jumps) == 1
    
    # Test with RegularJump
    function rj_rate(out, u, p, t)
        out[1] = p[1] * u[1]
    end
    
    rj_dc = zeros(3, 1)
    rj_dc[1, 1] = -1
    rj_dc[3, 1] = 1
    
    function rj_c(du, u, p, t, counts, mark)
        mul!(du, rj_dc, counts)
    end
    
    regj = RegularJump(rj_rate, rj_c, 1)
    
    jp_pure_regj = JumpProblem(prob, PureLeaping(), JumpSet(regj))
    @test jp_pure_regj.aggregator isa PureLeaping
    @test jp_pure_regj.discrete_jump_aggregation === nothing
    @test jp_pure_regj.regular_jump !== nothing
    
    # Test mixed jump types
    mixed_jumps = JumpSet(; massaction_jumps = maj, constant_jumps = (crj,), 
        variable_jumps = (vrj,), regular_jumps = regj)
    jp_pure_mixed = JumpProblem(prob, PureLeaping(), mixed_jumps)
    @test jp_pure_mixed.aggregator isa PureLeaping
    @test jp_pure_mixed.discrete_jump_aggregation === nothing
    @test jp_pure_mixed.massaction_jump !== nothing
    @test length(jp_pure_mixed.constant_jumps) == 1
    @test length(jp_pure_mixed.variable_jumps) == 1
    @test jp_pure_mixed.regular_jump !== nothing
    
    # Test spatial system error
    spatial_sys = CartesianGrid((2, 2))
    hopping_consts = [1.0]
    @test_throws ErrorException JumpProblem(prob, PureLeaping(), JumpSet(maj); 
                                          spatial_system = spatial_sys)
    @test_throws ErrorException JumpProblem(prob, PureLeaping(), JumpSet(maj); 
                                          hopping_constants = hopping_consts)
    
    # Test MassActionJump with parameter mapping
    maj_params = MassActionJump(reactant_stoich, net_stoich; param_idxs = [1, 2])
    jp_params = JumpProblem(prob, PureLeaping(), JumpSet(maj_params))
    scaled_rates = [p[1], p[2]/2]
    @test jp_params.massaction_jump.scaled_rates == scaled_rates
end

# Test SimpleTauLeaping saving controls
@testset "SimpleTauLeaping Saving Controls Tests" begin
    # Setup common test problem
    u0 = [100, 50, 0]
    tspan = (0.0, 10.0)
    prob = DiscreteProblem(u0, tspan)
    
    function test_rate(out, u, p, t)
        out[1] = 0.1 * u[1]  
        out[2] = 0.05 * u[2]
    end
    
    test_dc = zeros(3, 2)
    test_dc[1, 1] = -1
    test_dc[2, 1] = 1
    test_dc[2, 2] = -1  
    test_dc[3, 2] = 1
    
    function test_c(du, u, p, t, counts, mark)
        mul!(du, test_dc, counts)
    end
    
    rj = RegularJump(test_rate, test_c, 2)
    jump_prob = JumpProblem(prob, PureLeaping(), rj; rng=StableRNG(12345))
    dt = 0.5
    
    @testset "Basic saving controls" begin
        # Test save_everystep = false
        sol_no_steps = solve(jump_prob, SimpleTauLeaping(); dt=dt, save_everystep=false)
        @test length(sol_no_steps.t) == 2  # Only start and end
        @test sol_no_steps.t[1] == tspan[1]
        @test sol_no_steps.t[end] == tspan[2]
        
        # Test save_start = false 
        sol_no_start = solve(jump_prob, SimpleTauLeaping(); dt=dt, save_start=false)
        @test sol_no_start.t[1] != tspan[1] || length(sol_no_start.t) == 1
        
        # Test save_end = false
        sol_no_end = solve(jump_prob, SimpleTauLeaping(); dt=dt, save_end=false)
        @test sol_no_end.t[end] != tspan[2]
        
        # Test save_start = false, save_end = false
        sol_no_bounds = solve(jump_prob, SimpleTauLeaping(); dt=dt, 
                             save_start=false, save_end=false)
        expected_length = Int(ceil((tspan[2] - tspan[1]) / dt)) - 1
        @test length(sol_no_bounds.t) == expected_length
    end
    
    @testset "saveat functionality" begin
        # Test saveat as array
        saveat_times = [1.0, 3.0, 5.0, 7.0]
        sol_saveat = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=saveat_times, 
                          save_everystep=false)
        @test length(sol_saveat.t) >= length(saveat_times) + 2  # saveat + start + end
        for save_t in saveat_times
            @test any(abs.(sol_saveat.t .- save_t) .< dt/1e10)
        end
        
        # Test saveat as number (range)
        saveat_dt = 1.0
        sol_saveat_range = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=saveat_dt,
                                save_everystep=false) 
        expected_saves = collect(tspan[1]:saveat_dt:tspan[2])
        @test length(sol_saveat_range.t) >= length(expected_saves)
        
        # Test saveat with save_everystep=true (should include both)
        sol_both = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=[2.5, 7.5])
        regular_times = tspan[1]:dt:tspan[2]
        # Should have regular grid points plus saveat points
        @test length(sol_both.t) >= length(regular_times)
        @test any(abs.(sol_both.t .- 2.5) .< dt/1e10)
        @test any(abs.(sol_both.t .- 7.5) .< dt/1e10)
    end
    
    @testset "Edge cases and floating point handling" begin
        # Test saveat exactly at regular grid points
        regular_grid = collect(tspan[1]:dt:tspan[2])
        sol_exact = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=regular_grid[3:5])
        # Should not have duplicate times
        @test length(sol_exact.t) == length(unique(sol_exact.t))
        
        # Test saveat slightly offset from grid
        offset_times = regular_grid[2:4] .+ dt/1e12  # Very small offset
        sol_offset = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=offset_times)
        # Should merge with regular grid due to dtmin tolerance
        regular_sol = solve(jump_prob, SimpleTauLeaping(); dt=dt)
        @test length(sol_offset.t) == length(regular_sol.t)
        
        # Test custom dtmin
        large_dtmin = dt/100
        offset_times_large = regular_grid[2:4] .+ large_dtmin/2
        sol_custom_dtmin = solve(jump_prob, SimpleTauLeaping(); dt=dt, 
                                saveat=offset_times_large, dtmin=large_dtmin)
        # Should merge due to larger tolerance
        @test length(sol_custom_dtmin.t) == length(regular_sol.t)
        
        # Test saveat at boundaries
        boundary_times = [tspan[1], tspan[2]]
        sol_boundaries = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=boundary_times,
                              save_everystep=false)
        @test sol_boundaries.t[1] == tspan[1]
        @test sol_boundaries.t[end] == tspan[2]
        
        # Test saveat before start time (should save initial condition)
        early_times = [tspan[1] - 1.0, tspan[1] - 0.5, 2.0]
        sol_early = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=early_times,
                         save_everystep=false, save_start=false)
        @test length(sol_early.t) >= 3  # 2 early times + 1 valid + end
        @test any(sol_early.t .<= tspan[1])
        
        # Test empty saveat
        sol_empty = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=Float64[])
        regular_sol = solve(jump_prob, SimpleTauLeaping(); dt=dt)
        @test length(sol_empty.t) == length(regular_sol.t)
    end
    
    @testset "Solution correctness" begin
        # Test that solutions use same RNG sequence when settings are equivalent
        jump_prob1 = JumpProblem(prob, PureLeaping(), rj; rng=StableRNG(12345))
        jump_prob2 = JumpProblem(prob, PureLeaping(), rj; rng=StableRNG(12345))
        
        sol_default = solve(jump_prob1, SimpleTauLeaping(); dt=dt)
        sol_all_saves = solve(jump_prob2, SimpleTauLeaping(); dt=dt, 
                             saveat=collect(tspan[1]:dt:tspan[2]))
        
        # Should have same number of time points
        @test length(sol_default.t) == length(sol_all_saves.t)
        
        # Values at corresponding times should be identical with same RNG seed
        for i in 1:length(sol_default.t)
            @test isapprox(sol_default.t[i], sol_all_saves.t[i], atol=dt/1e10)
            @test sol_default.u[i] == sol_all_saves.u[i]  
        end
        
        # Test that no duplicate time points exist
        @test length(sol_default.t) == length(unique(sol_default.t))
        @test issorted(sol_default.t)
    end
    
    @testset "Combination tests" begin
        # Test all controls together
        sol_combined = solve(jump_prob, SimpleTauLeaping(); dt=dt,
                           save_start=false, save_end=false, 
                           save_everystep=false, saveat=[2.0, 4.0, 6.0])
        @test length(sol_combined.t) == 3  # Only saveat points
        @test all(t -> t in [2.0, 4.0, 6.0], sol_combined.t)
        
        # Test precedence: saveat should override regular grid when coincident
        saveat_at_grid = [2.0, 4.0]  # Assuming these align with regular grid
        sol_precedence = solve(jump_prob, SimpleTauLeaping(); dt=dt, saveat=saveat_at_grid)
        @test any(abs.(sol_precedence.t .- 2.0) .< dt/1e10)
        @test any(abs.(sol_precedence.t .- 4.0) .< dt/1e10)
        @test length(sol_precedence.t) == length(unique(sol_precedence.t))
    end
end
