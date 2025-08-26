struct SimpleTauLeaping <: DiffEqBase.DEAlgorithm end

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

function get_next_time(t_current, dt, regular_grid_next, saveat_times, 
                       saveat_idx, tspan_end, save_everystep, dtmin, use_saveat)
    next_time = tspan_end
    
    # Check next regular grid point if saving every step
    if save_everystep && regular_grid_next < tspan_end
        next_time = min(next_time, regular_grid_next)
    end
    
    # Check next saveat point if exists
    if use_saveat && saveat_idx <= length(saveat_times)
        saveat_t = saveat_times[saveat_idx]
        if saveat_t > t_current && saveat_t <= tspan_end
            # Only add saveat time if it doesn't coincide with regular grid
            if !(save_everystep && abs(saveat_t - regular_grid_next) < dtmin)
                next_time = min(next_time, saveat_t)
            end
        end
    end
    
    return next_time
end

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleTauLeaping;
        seed = nothing,
        dt = error("dt is required for SimpleTauLeaping."),
        save_start = true,
        save_end = true,
        saveat = nothing,
        save_everystep = saveat === nothing,
        dtmin = nothing)
    validate_pure_leaping_inputs(jump_prob, alg) ||
        error("SimpleTauLeaping can only be used with PureLeaping JumpProblems with only non-RegularJumps.")
    prob = jump_prob.prob
    rng = DEFAULT_RNG
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
    
    # Set default dtmin if not provided
    if dtmin === nothing
        dtmin = dt / 1e10
    end

    # Preprocess saveat
    use_saveat = saveat !== nothing
    saveat_times = if use_saveat
        if saveat isa Number
            collect(tspan[1]:saveat:tspan[2])
        else
            saveat
        end
    else
        typeof(tspan[1])[]  # Empty array with same type as t
    end

    # Initialize solution arrays
    t_vals = typeof(tspan[1])[]
    u_vals = typeof(u0)[]
    
    # Initialize time variables
    t = tspan[1]
    u = copy(u0)
    regular_grid_next = t + dt
    saveat_idx = 1
    
    # Save early saveat points at initial condition
    if use_saveat
        while saveat_idx <= length(saveat_times) && 
              saveat_times[saveat_idx] <= t + dtmin
            push!(t_vals, saveat_times[saveat_idx])
            push!(u_vals, copy(u0))
            saveat_idx += 1
        end
    end
    
    # Save initial condition if not already saved via saveat
    if save_start && (isempty(t_vals) || abs(t_vals[end] - t) > dtmin)
        push!(t_vals, t)
        push!(u_vals, copy(u))
    end

    # iteration variables
    counts = zero(rate_cache) # counts for each variable
    mark = nothing  # marks are not supported in SimpleTauLeaping

    # Main time-stepping loop
    while t < tspan[2] - dtmin
        # Determine next time point
        tnext = get_next_time(t, dt, regular_grid_next, saveat_times, 
                             saveat_idx, tspan[2], save_everystep, dtmin, use_saveat)
        current_dt = tnext - t
        
        # Skip negligible steps
        if current_dt < dtmin
            continue
        end
        
        # Perform tau-leaping with current_dt
        rate(rate_cache, u, p, t)
        rate_cache .*= current_dt # multiply by the width of the time interval
        counts .= pois_rand.((rng,), rate_cache) # set counts to the poisson arrivals with our given rates
        c(du, u, p, t, counts, mark)
        u = du + u
        t = tnext
        
        # Determine if we should save
        should_save = false
        is_end_time = abs(t - tspan[2]) < dtmin
        
        # Check if at regular grid point
        if save_everystep && abs(t - regular_grid_next) < dtmin
            should_save = !is_end_time || save_end  # Only save end if save_end=true
            regular_grid_next += dt
            
            # Check for coinciding saveat
            if use_saveat && saveat_idx <= length(saveat_times) &&
               abs(t - saveat_times[saveat_idx]) < dtmin
                saveat_idx += 1
            end
        # Check if at saveat point (not coinciding)
        elseif use_saveat && saveat_idx <= length(saveat_times) &&
               abs(t - saveat_times[saveat_idx]) < dtmin
            should_save = !is_end_time || save_end  # Only save end if save_end=true
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
        c(du, u, p, t, counts, mark)
        u = du + u
        t = tspan[2]
    end
    
    # Always save final time point if requested and not already saved
    if save_end && (isempty(t_vals) || abs(t_vals[end] - tspan[2]) > dtmin)
        push!(t_vals, tspan[2])
        push!(u_vals, copy(u))
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
