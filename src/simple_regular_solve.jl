struct SimpleTauLeaping <: DiffEqBase.DEAlgorithm end

function validate_pure_leaping_inputs(jump_prob::JumpProblem, alg)
    if !(jump_prob.aggregator isa PureLeaping)
        @warn "When using $alg, please pass PureLeaping() as the aggregator to the \
        JumpProblem, i.e. call JumpProblem(::DiscreteProblem, PureLeaping(),...). \
        Passing $(jump_prob.aggregator) is deprecated and will be removed in the next breaking release."
    end
    isvalid = isempty(jump_prob.jump_callback.continuous_callbacks) &&
        isempty(jump_prob.jump_callback.discrete_callbacks) &&
        isempty(jump_prob.constant_jumps) &&
        isempty(jump_prob.variable_jumps) &&
        get_num_majumps(jump_prob.massaction_jump) == 0 &&
        jump_prob.regular_jump !== nothing    
    
    if isvalid    
        (jump_prob.regular_jump.mark_dist === nothing) || 
            error("Mark distributions are currently not supported in SimpleTauLeaping")
    end

    return isvalid
end

function get_next_time(t_mesh, save_data, tspan_end)
    next_time = tspan_end
    
    # Always check next regular grid point (for numerical accuracy)
    (t_mesh < tspan_end) && (next_time = t_mesh)
    
    # Check next saveat point if exists
    @unpack use_saveat, saveat_idx, saveat_times = save_data
    if use_saveat && saveat_idx <= length(saveat_times)
        next_time = min(next_time, saveat_times[saveat_idx])
    end
    
    return next_time
end

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleTauLeaping;
        seed = nothing,
        dt = error("dt is required for SimpleTauLeaping."),
        saveat = nothing,
        save_everystep = nothing,
        save_start = nothing,
        save_end = nothing,
        dtmin = max(dt / 10^10, eps(dt)) )

    # validation
    validate_pure_leaping_inputs(jump_prob, alg) ||
        error("SimpleTauLeaping can only be used with PureLeaping JumpProblems with only non-RegularJumps.")

    rj = jump_prob.regular_jump
    @unpack rate, numjumps = rj    
    affects! = rj.c
    @unpack prob, rng = jump_prob
    (seed !== nothing) && seed!(rng, seed)
    @unpack tspan, p = prob
    u0 = copy(prob.u0)
    u = copy(u0)
    du = similar(u0)
    rate_cache = zeros(typeof(dt), numjumps)

    # Initialize time variables
    t = tspan[1]
    t_mesh = t + dt

    # setup saving
    save_data = initialize_saving(prob, dt, dtmin, saveat, save_everystep, save_start, 
        save_end)
    initial_save!(save_data, u, t)

    # iteration variables
    counts = zero(rate_cache) # counts for each variable

    # Main time-stepping loop
    while t < tspan[2] - dtmin
        # Determine next time point
        tnext = get_next_time(t_mesh, save_data, tspan[2])
        current_dt = tnext - t
        
        # Skip negligible steps
        if current_dt < dtmin
            continue
        end
        
        # Perform tau-leaping with current_dt
        rate(rate_cache, u, p, t)
        rate_cache .*= current_dt # multiply by the width of the time interval
        counts .= pois_rand.((rng,), rate_cache) # set counts to the poisson arrivals with our given rates
        affects!(du, u, p, t, counts, nothing)
        u = du + u
        t = tnext
        
        # Update regular grid tracker (always advance regardless of saving)
        if abs(t - t_mesh) < dtmin
            t_mesh += dt
        end
        
        # Handle saving using new API
        check_and_save!(save_data, u, t, dt, dtmin, t_mesh, tspan[2])
    end

    # Handle final step if not reached exactly
    if t < tspan[2] - dtmin
        # Take final step to tspan[2]
        final_dt = tspan[2] - t
        rate(rate_cache, u, p, t)
        rate_cache .*= final_dt
        counts .= pois_rand.((rng,), rate_cache)
        affects!(du, u, p, t, counts, mark)
        u = du + u
        t = tspan[2]
    end
    
    # Finalize saving using new API
    finalize_saving!(save_data, u, t, dtmin)

    sol = DiffEqBase.build_solution(prob, alg, save_data.t_vals, save_data.u_vals,
        calculate_error = false,
        interp = DiffEqBase.ConstantInterpolation(save_data.t_vals, save_data.u_vals))
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
