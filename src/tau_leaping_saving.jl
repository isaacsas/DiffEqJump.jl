# Tau-leaping saving utilities
# This module provides a reusable API for saving functionality in tau-leaping solvers

"""
    TauLeapingSaveData{T, U}

Mutable struct to hold all saving-related data and state for tau-leaping algorithms.

# Fields
- `t_vals::Vector{T}`: Storage for time points
- `u_vals::Vector{U}`: Storage for state values  
- `saveat_times::Vector{T}`: Unified vector of all save times
- `save_everystep::Bool`: Whether to save every step
- `save_start::Bool`: Whether to save initial condition
- `save_end::Bool`: Whether to save final condition
- `saveat_idx::Int`: Current index in saveat_times
"""
mutable struct TauLeapingSaveData{T, U}
    # Solution storage
    t_vals::Vector{T}
    u_vals::Vector{U}
    
    # Save configuration
    saveat_times::Vector{T}
    save_everystep::Bool
    save_start::Bool  
    save_end::Bool
    
    # State tracking
    saveat_idx::Int
end

"""
    initialize_saving(prob, dt, dtmin; saveat=nothing, save_everystep=nothing, save_start=nothing, save_end=nothing)

Initialize all saving-related data structures and parameters.

# Arguments
- `prob`: The problem being solved
- `dt`: Time step size (can be nothing for adaptive stepping)
- `dtmin`: Minimum time step for tolerance comparisons
- `saveat`: Times to save at (defaults to empty)
- `save_everystep`: Whether to save every step (defaults based on dt and saveat)
- `save_start`: Whether to save initial condition (defaults to true)
- `save_end`: Whether to save final condition (defaults to true)

# Returns
- `TauLeapingSaveData`: Initialized saving data structure
"""
function initialize_saving(prob, dt, dtmin; saveat=nothing, save_everystep=nothing, save_start=nothing, save_end=nothing)
    tspan = prob.tspan
    u0 = prob.u0
    
    # Handle default kwargs
    if saveat === nothing
        saveat = typeof(tspan[1])[]
    end
    
    if save_everystep === nothing
        save_everystep = dt === nothing ? false : isempty(saveat)
    end
    
    if save_start === nothing
        save_start = save_everystep || isempty(saveat) || saveat isa Number ? true : tspan[1] in saveat
    end
    
    if save_end === nothing  
        save_end = save_everystep || isempty(saveat) || saveat isa Number ? true : tspan[2] in saveat
    end
    
    # Preprocess and unify all save times
    if saveat isa Number
        saveat_times = collect(tspan[1]:saveat:tspan[2])
    else
        saveat_times = copy(saveat)
        !issorted(saveat_times) && sort!(saveat_times)
        
        # Validate that all saveat times are within tspan bounds
        if !isempty(saveat_times) && (saveat_times[1] < tspan[1] || saveat_times[end] > tspan[2])
            error("All saveat times must be within the time span [$(tspan[1]), $(tspan[2])]")
        end
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
    
    # Initialize storage arrays
    t_vals = typeof(tspan[1])[]
    u_vals = typeof(u0)[]
    
    return TauLeapingSaveData(t_vals, u_vals, saveat_times, save_everystep, save_start, save_end, 1)
end

"""
    initial_save!(save_data, u0, t0)

Handle initial condition saving.

Saves all saveat times at or before the initial time using the initial condition u0.
Updates saveat_idx to point to the next unsaved time.

# Arguments
- `save_data`: TauLeapingSaveData struct
- `u0`: Initial state
- `t0`: Initial time
"""
function initial_save!(save_data::TauLeapingSaveData, u0, t0)
    dtmin = eps(typeof(t0)) * 100  # Small tolerance for time comparisons
    
    # Save all saveat times at or before initial time (with initial condition)
    while !isempty(save_data.saveat_times) && save_data.saveat_idx <= length(save_data.saveat_times) &&
          save_data.saveat_times[save_data.saveat_idx] <= t0 + dtmin
        push!(save_data.t_vals, save_data.saveat_times[save_data.saveat_idx])
        push!(save_data.u_vals, copy(u0))
        save_data.saveat_idx += 1
    end
    
    return nothing
end

"""
    check_and_save!(save_data, u, t, dt, dtmin, regular_grid_next, tspan_end)

Handle saving within the timestepping loop.

Checks if the current time point should be saved (regular grid or saveat point) and
performs saves when conditions are met. Updates saveat_idx when saveat points are saved.

# Arguments
- `save_data`: TauLeapingSaveData struct
- `u`: Current state
- `t`: Current time
- `dt`: Time step size
- `dtmin`: Minimum time step for tolerance comparisons
- `regular_grid_next`: Next regular grid point
- `tspan_end`: End time of integration

# Returns
- `Bool`: Whether a save occurred
"""
function check_and_save!(save_data::TauLeapingSaveData, u, t, dt, dtmin, regular_grid_next, tspan_end)
    saved = false
    
    # Check if we should save at this time point
    should_save_regular = save_data.save_everystep && abs(t - regular_grid_next + dt) <= dtmin
    should_save_saveat = false
    
    # Check if current time matches a saveat point
    if !isempty(save_data.saveat_times) && save_data.saveat_idx <= length(save_data.saveat_times)
        saveat_t = save_data.saveat_times[save_data.saveat_idx]
        if abs(t - saveat_t) <= dtmin
            should_save_saveat = true
        end
    end
    
    # Special handling for end time when save_end=false
    is_end_time = abs(t - tspan_end) <= dtmin
    if is_end_time && !save_data.save_end
        should_save_regular = false
        should_save_saveat = false
    end
    
    # Perform save if needed
    if should_save_regular || should_save_saveat
        push!(save_data.t_vals, t)
        push!(save_data.u_vals, copy(u))
        saved = true
        
        # Update saveat index if we saved a saveat point
        if should_save_saveat
            save_data.saveat_idx += 1
        end
    end
    
    return saved
end

"""
    finalize_saving!(save_data, u_final, t_final, dtmin)

Complete the saving process after timestepping.

Saves any remaining saveat times (including final time if save_end was added)
using the final state for all remaining save points.

# Arguments
- `save_data`: TauLeapingSaveData struct
- `u_final`: Final state
- `t_final`: Final time
- `dtmin`: Minimum time step for tolerance comparisons
"""
function finalize_saving!(save_data::TauLeapingSaveData, u_final, t_final, dtmin)
    # Save any remaining saveat times (including final time if save_end was added)
    while !isempty(save_data.saveat_times) && save_data.saveat_idx <= length(save_data.saveat_times)
        saveat_t = save_data.saveat_times[save_data.saveat_idx]
        # Only save if the saveat time is at or after current time
        if saveat_t >= t_final - dtmin
            push!(save_data.t_vals, saveat_t)
            push!(save_data.u_vals, copy(u_final))
        end
        save_data.saveat_idx += 1
    end
    
    return nothing
end