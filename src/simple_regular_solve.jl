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

function get_next_time(t_current, dt, regular_grid_next, saveat_times, 
                       saveat_idx, tspan_end, save_everystep, dtmin)
    next_time = tspan_end
    
    # Check next regular grid point if saving every step
    if save_everystep && regular_grid_next < tspan_end
        next_time = min(next_time, regular_grid_next)
    end
    
    # Check next saveat point if exists
    if !isempty(saveat_times) && saveat_idx <= length(saveat_times)
        saveat_t = saveat_times[saveat_idx]
        if saveat_t > t_current && saveat_t <= tspan_end
            next_time = min(next_time, saveat_t)
        end
    end
    
    return next_time
end

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleTauLeaping;
        seed = nothing,
        dt = error("dt is required for SimpleTauLeaping."),
        saveat = typeof(dt)[],
        save_everystep = isempty(saveat),  # note isempty(scalar) == false in Julia!
        save_start = save_everystep || isempty(saveat) || saveat isa Number ? true : jump_prob.prob.tspan[1] in saveat,
        save_end = save_everystep || isempty(saveat) || saveat isa Number ? true : jump_prob.prob.tspan[2] in saveat,
        dtmin = nothing)

    # validation
    validate_pure_leaping_inputs(jump_prob, alg) ||
        error("SimpleTauLeaping can only be used with PureLeaping JumpProblems with only non-RegularJumps.")

    rj = jump_prob.regular_jump
    @unpack rate, numjumps = rj    
    affects! = rj.c
    prob = jump_prob.prob
    @unpack tspan, p = prob
    rng = jump_prob.rng
    (seed !== nothing) && seed!(rng, seed)
    u0 = copy(prob.u0)
    u = copy(u0)
    du = similar(u0)
    rate_cache = zeros(float(eltype(u0)), numjumps)

    # Set default dtmin if not provided
    if dtmin === nothing
        dtmin = dt / 1e10
    end

    # solution output arrays
    t_vals = typeof(tspan[1])[]
    u_vals = typeof(u0)[]
    
    # Initialize time variables
    t = tspan[1]
    regular_grid_next = t + dt

    # Preprocess and unify all save times
    if saveat isa Number
        saveat_times = collect(tspan[1]:saveat:tspan[2])
    else
        saveat_times = copy(saveat)
        !issorted(saveat_times) && sort!(saveat_times)
    end
    
    # Handle save_start and save_end precedence over saveat
    if save_start
        # Add start time if not already present
        if isempty(saveat_times) || saveat_times[1] > tspan[1] + dtmin
            pushfirst!(saveat_times, tspan[1])
        end
    else
        # Remove start time if present in saveat (save_start=false takes precedence)
        if !isempty(saveat_times) && abs(saveat_times[1] - tspan[1]) <= dtmin
            popfirst!(saveat_times)
        end
    end
    
    if save_end
        # Add end time if not already present  
        if isempty(saveat_times) || saveat_times[end] < tspan[2] - dtmin
            push!(saveat_times, tspan[2])
        end
    else
        # Remove end time if present in saveat (save_end=false takes precedence)
        if !isempty(saveat_times) && abs(saveat_times[end] - tspan[2]) <= dtmin
            pop!(saveat_times)
        end
    end
    
    # Remove duplicates within dtmin tolerance
    if length(saveat_times) > 1
        unique_times = [saveat_times[1]]
        for i in 2:length(saveat_times)
            if abs(saveat_times[i] - unique_times[end]) > dtmin
                push!(unique_times, saveat_times[i])
            end
        end
        saveat_times = unique_times
    end
    
    # Handle initial saves
    saveat_idx = 1
    
    # Save all saveat times at or before initial time (with initial condition)
    while !isempty(saveat_times) && saveat_idx <= length(saveat_times) &&
          saveat_times[saveat_idx] <= t + dtmin
        push!(t_vals, saveat_times[saveat_idx])
        push!(u_vals, copy(u))
        saveat_idx += 1
    end

    # iteration variables
    counts = zero(rate_cache) # counts for each variable
    mark = nothing  # marks are not supported in SimpleTauLeaping

    # Main time-stepping loop
    while t < tspan[2] - dtmin
        # Determine next time point
        tnext = get_next_time(t, dt, regular_grid_next, saveat_times, 
                             saveat_idx, tspan[2], save_everystep, dtmin)
        current_dt = tnext - t
        
        # Skip negligible steps
        if current_dt < dtmin
            continue
        end
        
        # Perform tau-leaping with current_dt
        rate(rate_cache, u, p, t)
        rate_cache .*= current_dt # multiply by the width of the time interval
        counts .= pois_rand.((rng,), rate_cache) # set counts to the poisson arrivals with our given rates
        affects!(du, u, p, t, counts, mark)
        u = du + u
        t = tnext
        
        # Unified saving logic
        should_save = false
        
        # Check if at regular grid point
        if save_everystep && abs(t - regular_grid_next) < dtmin
            # Don't save end time if save_end=false
            is_end_time = abs(t - tspan[2]) < dtmin
            if !(is_end_time && !save_end)
                should_save = true
            end
            regular_grid_next += dt
        end
        
        # Check if at saveat point
        if !isempty(saveat_times) && saveat_idx <= length(saveat_times) &&
           abs(t - saveat_times[saveat_idx]) < dtmin
            should_save = true
            saveat_idx += 1
        end
        
        if should_save
            push!(t_vals, t)
            push!(u_vals, copy(u))
        end
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
    
    # Save any remaining saveat times (including final time if save_end was added)
    while !isempty(saveat_times) && saveat_idx <= length(saveat_times)
        save_t = saveat_times[saveat_idx]
        if save_t >= t - dtmin  # Only save times at or after final time
            push!(t_vals, save_t)
            push!(u_vals, copy(u))  # Use final state for any remaining saves
        end
        saveat_idx += 1
    end

    sol = DiffEqBase.build_solution(prob, alg, t_vals, u_vals,
        calculate_error = false,
        interp = DiffEqBase.ConstantInterpolation(t_vals, u_vals))
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
