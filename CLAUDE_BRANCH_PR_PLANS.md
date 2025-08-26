# Implementation Plan: Add Saving Controls to SimpleTauLeaping Methods

## Overview
Add support for the integrator saving controls (`saveat`, `save_everystep`, `save_start`, and `save_end`) to both `SimpleTauLeaping` and `SimpleSplitTauLeaping` methods, consistent with their documented behavior via the common solver interface at https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/

## Saving Controls Behavior

1. **`saveat`**: Specifies specific times to save the solution
   - Can be an array of times or a number (which expands to a range)
   - Default is `nothing`

2. **`save_everystep`**: Saves the result at every regular dt step
   - Default is `true`

3. **`save_start`**: Determines whether the initial condition is included
   - Default is `true`

4. **`save_end`**: Forces the final timepoint to be saved
   - Default is `true`

## Implementation Details

### 1. Update Function Signatures

Update both methods in `src/simple_regular_solve.jl`:

```julia
function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleTauLeaping;
        seed = nothing, 
        dt = error("dt is required for SimpleTauLeaping."),
        save_start = true, 
        save_end = true, 
        save_everystep = true, 
        saveat = nothing,
        dtmin = nothing)

function DiffEqBase.solve(jump_prob::JumpProblem, alg::SimpleSplitTauLeaping;
        seed = nothing,
        dt = error("dt is required for SimpleSplitTauLeaping."),
        save_start = true,
        save_end = true,
        save_everystep = true,
        saveat = nothing,
        dtmin = nothing)
```

### 2. Adaptive Time-Stepping Strategy

Instead of fixed dt steps, dynamically determine the next time point based on:
- Regular dt grid points (if `save_everystep=true`)
- `saveat` points
- Final time (`tspan[2]`)

The algorithm steps directly to the next required time point with an appropriate step size, avoiding interpolation.

### 3. Floating Point Tolerance Handling

- Use `dtmin` parameter for floating point comparisons (default: `dt/10^10`)
- If `abs(saveat_time - regular_dt_time) < dtmin`, treat as the same time point
- Prevents duplicate saves at essentially the same time
- Handles edge cases where `saveat` values coincide with regular timesteps

### 4. Core Algorithm Structure

```julia
# Initialize
t = tspan[1]
regular_grid_next = t + dt
saveat_idx = 1
t_vals = Vector{typeof(t)}()
u_vals = Vector{typeof(u0)}()

# Preprocess saveat
if saveat isa Number
    saveat_times = collect(tspan[1]:saveat:tspan[2])
else
    saveat_times = saveat
end

# Skip saveat points at or before start
if saveat_times !== nothing
    while saveat_idx <= length(saveat_times) && 
          saveat_times[saveat_idx] <= t + dtmin
        saveat_idx += 1
    end
end

# Save initial condition
if save_start
    push!(t_vals, t)
    push!(u_vals, copy(u))
end

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
    # [tau-leaping logic here]
    
    t = tnext
    
    # Determine if we should save
    should_save = false
    
    # Check if at regular grid point
    if save_everystep && abs(t - regular_grid_next) < dtmin
        should_save = true
        regular_grid_next += dt
        
        # Check for coinciding saveat
        if saveat_times !== nothing && saveat_idx <= length(saveat_times) &&
           abs(t - saveat_times[saveat_idx]) < dtmin
            saveat_idx += 1
        end
    # Check if at saveat point (not coinciding)
    elseif saveat_times !== nothing && saveat_idx <= length(saveat_times) &&
           abs(t - saveat_times[saveat_idx]) < dtmin
        should_save = true
        saveat_idx += 1
    end
    
    if should_save
        push!(t_vals, t)
        push!(u_vals, copy(u))
    end
end

# Handle save_end
if save_end && (isempty(t_vals) || abs(t_vals[end] - tspan[2]) > dtmin)
    # Final step to tspan[2] if needed
    # [final step logic]
    push!(t_vals, tspan[2])
    push!(u_vals, copy(u))
end
```

### 5. Helper Function for Time Determination

```julia
function get_next_time(t_current, dt, regular_grid_next, saveat_times, 
                       saveat_idx, tspan_end, save_everystep, dtmin)
    candidates = [tspan_end]
    
    # Add next regular grid point if saving every step
    if save_everystep && regular_grid_next < tspan_end
        push!(candidates, regular_grid_next)
    end
    
    # Add next saveat point if exists
    if saveat_times !== nothing && saveat_idx <= length(saveat_times)
        saveat_t = saveat_times[saveat_idx]
        if saveat_t > t_current && saveat_t <= tspan_end
            # Check if saveat coincides with regular grid
            if save_everystep && abs(saveat_t - regular_grid_next) < dtmin
                # Use regular grid time instead
            else
                push!(candidates, saveat_t)
            end
        end
    end
    
    return minimum(candidates)
end
```

### 6. Method-Specific Implementations

#### SimpleTauLeaping
- Use `current_dt` to scale rates: `rate_cache .*= current_dt`
- Apply Poisson sampling with scaled rates
- Execute regular jump updates

#### SimpleSplitTauLeaping
- For each mass action jump:
  - Evaluate rate
  - Scale by `current_dt`
  - Sample Poisson for number of firings
  - Execute jump that many times

## Special Cases to Handle

1. `saveat` as a single number (convert to range)
2. Empty `saveat` array
3. `saveat` points outside tspan
4. `saveat` points very close to tspan boundaries
5. Multiple `saveat` points within `dtmin` of each other
6. Ensuring no duplicate time points in solution

## Testing Requirements

1. **Individual control tests**:
   - Test each saving control separately
   - Verify default behaviors

2. **Combination tests**:
   - Test various combinations of saving controls
   - Verify precedence and interaction

3. **Floating point tolerance tests**:
   - `saveat` exactly at regular grid points
   - `saveat` slightly offset from regular grid (< dtmin)
   - Custom `dtmin` values

4. **Edge cases**:
   - `saveat` at tspan boundaries
   - Empty arrays
   - Very small/large dt values

5. **Correctness verification**:
   - No duplicate time points
   - Correct Poisson rates with variable step sizes
   - Solution accuracy maintained

## Benefits

- Exact values at `saveat` times (no interpolation error)
- Maintains accuracy of tau-leaping method
- Flexible stepping adapts to user-specified save points
- Handles floating point precision issues robustly
- Compatible with DiffEq common solver interface
- User control via `dtmin` when needed

## Files to Modify

1. `src/simple_regular_solve.jl` - Main implementation
2. `test/regular_jumps.jl` - Add tests for new functionality
3. `test/split_tau_leaping.jl` - Add tests for SimpleSplitTauLeaping

## Backwards Compatibility

Default behavior remains unchanged - all saving controls have defaults that maintain current behavior when not specified.