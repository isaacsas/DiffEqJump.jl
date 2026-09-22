using JumpProcesses, OrdinaryDiffEq, StochasticDiffEq, SciMLBase
using Random, StableRNGs, Test

function stochastic_jump_problem(kind; seed = 0, variable = false, callback = nothing)
    f!(du, u, p, t) = (du .= -0.1 .* u)
    g!(du, u, p, t) = (du .= 0.2)
    f_rode!(du, u, p, t, W) = (du .= -0.1 .* u .+ 0.2 .* W)
    prob = kind === :SDE ?
        SDEProblem(f!, g!, [10.0], (0.0, 0.5), zeros(Int, 2); seed) :
        RODEProblem(f_rode!, [10.0], (0.0, 0.5), zeros(Int, 2); seed)
    rate(u, p, t) = 10.0
    affect!(integrator) = (integrator.u[1] += 1)
    jump = variable ? VariableRateJump(rate, affect!) : ConstantRateJump(rate, affect!)
    JumpProblem(prob, Direct(), jump; callback)
end

stochastic_rng_options(kind) =
    kind === :default ? (;) :
    kind === :seed ? (; seed = 456) :
    kind === :zero_seed ? (; seed = 0) :
    kind === :task_rng_seed ? (; rng = Random.default_rng(), seed = 456) :
    kind === :task_rng ? (; rng = Random.default_rng()) :
    (; rng = StableRNG(789), seed = 456)

@testset "Stochastic solve/init use the same RNG policy" begin
    @testset "$kind stored seed=$stored_seed options=$options_kind variable=$variable" for
            (kind, alg, variable) in ((:SDE, EM(), false), (:SDE, EM(), true),
                (:RODE, RandomEM(), false)),
            stored_seed in (0, 123),
            options_kind in (:default, :seed, :zero_seed, :task_rng_seed, :task_rng,
                :explicit_rng_seed)
        Random.seed!(999)
        sol = solve(stochastic_jump_problem(kind; seed = stored_seed, variable), alg;
            dt = 0.01, stochastic_rng_options(options_kind)...)

        Random.seed!(999)
        options = stochastic_rng_options(options_kind)
        integrator = init(stochastic_jump_problem(kind; seed = stored_seed, variable), alg;
            dt = 0.01, options...)
        if options_kind === :explicit_rng_seed
            @test SciMLBase.get_rng(integrator) === options.rng
        else
            @test SciMLBase.get_rng(integrator) isa Xoshiro
        end
        solve!(integrator)
        @test successful_retcode(sol)
        @test successful_retcode(integrator.sol)
        @test sol.t == integrator.sol.t
        @test sol.u == integrator.sol.u
        @test sol.seed == integrator.sol.seed
        expected_seed = options_kind === :explicit_rng_seed ? 0 :
            options_kind in (:seed, :task_rng_seed) ? 456 : stored_seed
        if expected_seed != 0 || options_kind === :explicit_rng_seed
            @test sol.seed == expected_seed
        end
    end
end

@testset "Stochastic initialization preserves callback merging" begin
    @testset "merge_callbacks=$merge_callbacks entry=$entry" for
            merge_callbacks in (true, false), entry in (:solve, :init)
        cb1 = DiscreteCallback((u, t, integrator) -> t == 0.25,
            integrator -> (integrator.p[1] += 1); save_positions = (false, false))
        cb2 = DiscreteCallback((u, t, integrator) -> t == 0.25,
            integrator -> (integrator.p[2] += 1); save_positions = (false, false))
        jprob = stochastic_jump_problem(:SDE; callback = cb1)
        options = (; seed = 123, dt = 0.01, tstops = [0.25], callback = cb2, merge_callbacks)
        sol = if entry === :solve
            solve(jprob, EM(); options...)
        else
            integrator = init(jprob, EM(); options...)
            solve!(integrator)
            integrator.sol
        end
        @test successful_retcode(sol)
        @test sol.prob.p == [Int(merge_callbacks), 1]
    end
end

@testset "Stochastic initialization respects the backend jump-copy option" begin
    for alias_jumps in (true, false)
        jprob = stochastic_jump_problem(:SDE)
        integrator = init(jprob, EM(); dt = 0.01, seed = 123,
            alias = SciMLBase.SDEAliasSpecifier(; alias_jumps))
        aliases_aggregation = any(integrator.opts.callback.discrete_callbacks) do callback
            callback.condition === jprob.discrete_jump_aggregation
        end
        @test aliases_aggregation == alias_jumps
        solve!(integrator)
        @test successful_retcode(integrator.sol)
    end
end

@testset "Stochastic initialization preserves explicit legacy alias_jump" begin
    for (kind, alg) in ((:SDE, EM()), (:RODE, RandomEM())),
            alias_jump in (true, false), native_alias in (nothing, true, false, :specifier)
        jprob = stochastic_jump_problem(kind)
        constructor = kind === :SDE ? SciMLBase.SDEAliasSpecifier : SciMLBase.RODEAliasSpecifier
        alias = native_alias === :specifier ?
            constructor(; alias_jumps = !alias_jump, alias_u0 = false) : native_alias
        integrator = init(jprob, alg; dt = 0.01, seed = 123, alias_jump, alias)
        aliases_aggregation = any(integrator.opts.callback.discrete_callbacks) do callback
            callback.condition === jprob.discrete_jump_aggregation
        end
        @test aliases_aggregation == alias_jump
        if native_alias in (false, :specifier)
            @test integrator.u !== jprob.prob.u0
        end
        solve!(integrator)
        @test successful_retcode(integrator.sol)
    end
end
