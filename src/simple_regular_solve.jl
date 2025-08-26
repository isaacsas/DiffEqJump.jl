struct SimpleTauLeaping <: DiffEqBase.DEAlgorithm end
struct SimpleSplitTauLeaping <: DiffEqBase.DEAlgorithm end

function validate_pure_leaping_inputs(jump_prob::JumpProblem, alg)
    if !(jump_prob.aggregator isa PureLeaping)
        @warn "When using $alg, please pass PureLeaping() as the aggregator to the \
        JumpProblem, i.e. call JumpProblem(::DiscreteProblem, PureLeaping(),...). \
        Passing $(jump_prob.aggregator) is deprecated and will be removed in the next breaking release."
    end
    isempty(jump_prob.jump_callback.continuous_callbacks) &&
    isempty(jump_prob.jump_callback.discrete_callbacks) &&
    isempty(jump_prob.constant_jumps) &&
    isempty(jump_prob.variable_jumps) &&
    get_num_majumps(jump_prob.massaction_jump) == 0 &&
    jump_prob.regular_jump !== nothing    
end

function validate_massjump_splitting_inputs(jump_prob::JumpProblem, alg)
    if !(jump_prob.aggregator isa PureLeaping)
        @warn "When using $alg, please pass PureLeaping() as the aggregator..."
    end
    # Only MassActionJumps allowed
    isempty(jump_prob.jump_callback.continuous_callbacks) &&
    isempty(jump_prob.jump_callback.discrete_callbacks) &&
    isempty(jump_prob.constant_jumps) &&
    isempty(jump_prob.variable_jumps) &&
    jump_prob.regular_jump === nothing &&
    get_num_majumps(jump_prob.massaction_jump) > 0
end

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleTauLeaping;
        seed = nothing, dt = error("dt is required for SimpleTauLeaping."))
    validate_pure_leaping_inputs(jump_prob, alg) ||
        error("SimpleTauLeaping can only be used with PureLeaping JumpProblems with only non-RegularJumps.")
    
    @unpack prob, rng = jump_prob
    (seed !== nothing) && seed!(rng, seed)

    rj = jump_prob.regular_jump
    rate = rj.rate # rate function rate(out,u,p,t)
    numjumps = rj.numjumps # used for size information (# of jump processes)
    c = rj.c # matrix-free operator c(u_buffer, uprev, tprev, counts, p, mark)

    if !isnothing(rj.mark_dist) == nothing # https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/250
        error("Mark distributions are currently not supported in SimpleTauLeaping")
    end

    u0 = copy(prob.u0)
    du = similar(u0)
    rate_cache = zeros(float(eltype(u0)), numjumps)

    tspan = prob.tspan
    p = prob.p

    n = Int((tspan[2] - tspan[1]) / dt) + 1
    u = Vector{typeof(prob.u0)}(undef, n)
    u[1] = u0
    t = tspan[1]:dt:tspan[2]

    # iteration variables
    counts = zero(rate_cache) # counts for each variable

    for i in 2:n # iterate over dt-slices
        uprev = u[i - 1]
        tprev = t[i - 1]
        rate(rate_cache, uprev, p, tprev)
        rate_cache .*= dt # multiply by the width of the time interval
        counts .= pois_rand.((rng,), rate_cache) # set counts to the poisson arrivals with our given rates
        c(du, uprev, p, tprev, counts, mark)
        u[i] = du + uprev
    end

    sol = DiffEqBase.build_solution(prob, alg, t, u,
        calculate_error = false,
        interp = DiffEqBase.ConstantInterpolation(t, u))
end

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleSplitTauLeaping;
        seed = nothing,
        dt = error("dt is required for SimpleSplitTauLeaping."))
    
    validate_massjump_splitting_inputs(jump_prob, alg) ||
        error("SimpleSplitTauLeaping currently only supports MassActionJumps with PureLeaping")
    
    @unpack prob, rng = jump_prob
    (seed !== nothing) && seed!(rng, seed)
    
    # Extract MassActionJumps
    ma_jumps = jump_prob.massaction_jump
    num_jumps = get_num_majumps(ma_jumps)
    
    # Pre-allocate
    u0 = copy(prob.u0)
    u_work = similar(u0)  # Working state vector
    
    tspan = prob.tspan
    p = prob.p
    n = Int((tspan[2] - tspan[1]) / dt) + 1
    u = Vector{typeof(u0)}(undef, n)
    u[1] = u0
    t = tspan[1]:dt:tspan[2]
    
    # Main loop - operator splitting with individual jump execution
    for i in 2:n
        copy!(u_work, u[i-1])
        
        # Split tau-leaping: evaluate and execute each jump separately
        @inbounds for j in 1:num_jumps
            rate = evalrxrate(u_work, j, ma_jumps)
            num_firings = pois_rand(rng, rate * dt)
            
            # Execute this jump num_firings times
            for _ in 1:num_firings
                executerx!(u_work, j, ma_jumps)
            end
        end
        
        u[i] = copy(u_work)
    end
    
    sol = DiffEqBase.build_solution(prob, alg, t, u,
        calculate_error = false,
        interp = DiffEqBase.ConstantInterpolation(t, u))
end

struct EnsembleGPUKernel{Backend} <: SciMLBase.EnsembleAlgorithm
    backend::Backend
    cpu_offload::Float64
end

function EnsembleGPUKernel(backend)
    EnsembleGPUKernel(backend, 0.0)
end

function EnsembleGPUKernel()
    EnsembleGPUKernel(nothing, 0.0)
end
