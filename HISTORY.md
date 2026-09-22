# Breaking updates and feature summaries across releases

## 10.0 (Breaking)

  - **Breaking**: The `rng` keyword argument has been removed from
    `JumpProblem`. Pass `rng` to `solve` or `init` instead:
    ```julia
    # Before (no longer works):
    jprob = JumpProblem(dprob, Direct(), jump; rng = Xoshiro(1234))
    sol = solve(jprob, SSAStepper())

    # After:
    jprob = JumpProblem(dprob, Direct(), jump)
    sol = solve(jprob, SSAStepper(); rng = Xoshiro(1234))
    ```
  - RNG state is now owned by the integrator, not the aggregator. This
    eliminates data races when sharing a `JumpProblem` across threads and
    ensures a single, consistent RNG priority across all solver pathways:
    `rng` > `seed` > `Random.default_rng()`.
  - `rng` and `seed` kwargs are fully supported on `solve`/`init` for all
    solver pathways (SSAStepper, ODE, SDE, tau-leaping).
  - `SSAIntegrator` now supports the `SciMLBase` RNG interface (`has_rng`,
    `get_rng`, `set_rng!`).
  - **Breaking**: The `scale_rates` and `useiszero` keyword arguments have been
    removed from `JumpProblem`. Set them on the `MassActionJump` directly:
    ```julia
    # Before (no longer works):
    jprob = JumpProblem(dprob, Direct(), maj; scale_rates = false)

    # After:
    maj = MassActionJump(rates, reactant_stoch, net_stoch; scale_rates = false)
    jprob = JumpProblem(dprob, Direct(), maj)
    ```
  - **Breaking**: Parameterized `MassActionJump`s (those constructed with
    `param_idxs` or a custom `param_mapper`) are now immutable — rates are
    computed from parameters at aggregator initialization rather than being
    materialized into the jump at `JumpProblem` construction time. This means:
      - `update_parameters!` has been removed. Mass action rates are now
        automatically recomputed from the current parameter values whenever the
        aggregator reinitializes. After modifying parameters (e.g. in a
        callback), call `reset_aggregated_jumps!(integrator)` to trigger
        reinitialization with the updated parameter values.
      - Custom parameter mappers (e.g. ModelingToolkitBase's
        `JumpSysMajParamMapper`) must implement the 3-arg callable API:
        `(mapper)(dest::AbstractVector, maj::MassActionJump, params)`.
        See [`MassActionJumpParamMapper`](@ref) for details.
  - Scalar `param_idxs` (e.g. `param_idxs = 1`) is now internally converted to
    a one-element vector. The scalar form continues to work as before.

## 9.33.1

  - Fixed `SSAStepper` saving the wrong state at the final requested `saveat`
    time when it precedes a later jump. The observation now contains the state
    at the requested time, and output times remain sorted when jump times are
    also saved.

## 9.33

  - Added user-specified rate bounds for `ConstantRateJump`s, enabling rates that
    are non-monotonic or increase in some species and decrease in others to be
    used with `RSSA` and `RSSACR` ([#653](https://github.com/SciML/JumpProcesses.jl/pull/653)).
    Pass `bounds(ulow, uhigh, u, p, t)` through the `bounds` keyword when
    constructing the jump. Here `ulow` and `uhigh` define the current species
    population bracket, and `u`, `p`, and `t` are the current state, parameters,
    and time. The function returns `RateBounds(; lrate, urate)` with nonnegative,
    finite bounds satisfying `lrate <= rate(v, p, t) <= urate` for every state
    `v` in that bracket. These bounds must hold throughout the bracket, not just
    at the current state. As with all `ConstantRateJump`s, the rate must remain
    constant between jumps and must not explicitly depend on time.

    For example, the following jump converts one particle of species 1 into
    species 2 at a rate that increases with species 1 and decreases with species 2:

    ```julia
    using JumpProcesses

    rate(u, p, t) = p[1] * u[1] / (1 + u[2])
    function affect!(integrator)
        integrator.u[1] -= 1
        integrator.u[2] += 1
        nothing
    end
    function bounds(ulow, uhigh, u, p, t)
        RateBounds(
            lrate = p[1] * ulow[1] / (1 + uhigh[2]),
            urate = p[1] * uhigh[1] / (1 + ulow[2])
        )
    end

    jump = ConstantRateJump(rate, affect!; bounds)
    prob = DiscreteProblem([10, 0], (0.0, 10.0), [1.0])

    # Each species affects the rate of jump 1; jump 1 changes both species.
    vartojumps_map = [[1], [1]]
    jumptovars_map = [[1, 2]]
    jprob = JumpProblem(prob, RSSA(), jump; vartojumps_map, jumptovars_map)
    sol = solve(jprob, SSAStepper())
    ```

    The same example works with `RSSACR()` in place of `RSSA()`. Omitting
    `bounds` preserves the existing behavior: rate bounds are computed by
    evaluating the rate at `ulow` and `uhigh`. This is valid for rates that are
    nondecreasing in every species or nonincreasing in every species, but is not
    generally valid for mixed dependence such as the example above.

## 9.14

  - Added the constant complexity next reaction method (CCNRM).

## 9.13

  - Added a default aggregator selection algorithm based on the number of passed
    in jumps. i.e. the following now auto-selects an aggregator (`Direct` in this
    case):
    
    ```julia
    using JumpProcesses
    rate(u, p, t) = u[1]
    affect(integrator) = (integrator.u[1] -= 1; nothing)
    crj = ConstantRateJump(rate, affect)
    dprob = DiscreteProblem([10], (0.0, 10.0))
    jprob = JumpProblem(dprob, crj)
    sol = solve(jprob, SSAStepper())
    ```

  - For `JumpProblem`s over `DiscreteProblem`s that only have `MassActionJump`s,
    `ConstantRateJump`s, and bounded `VariableRateJump`s, one no longer needs to
    specify `SSAStepper()` when calling `solve`, i.e. the following now works for
    the previous example and is equivalent to manually passing `SSAStepper()`:
    
    ```julia
    sol = solve(jprob)
    ```
  - Plotting a solution generated with `save_positions = (false, false)` now uses
    piecewise linear plots between any saved time points specified via `saveat`
    instead (previously the plots appeared piecewise constant even though each
    jump was not being shown). Note that solution objects still use piecewise
    constant interpolation, see [the
    docs](https://docs.sciml.ai/JumpProcesses/stable/tutorials/discrete_stochastic_example/#save_positions_docs)
    for details.

## 9.7

  - `Coevolve` was updated to support use with coupled ODEs/SDEs. See the updated
    documentation for details, and note the comments there about one needing to ensure
    rate bounds hold however the ODE/SDE stepper could modify dependent variables during a timestep.

## 9.3

  - Support for "bounded" `VariableRateJump`s that can be used with the `Coevolve`
    aggregator for faster simulation of jump processes with time-dependent rates.
    In particular, if all `VariableRateJump`s in a pure-jump system are bounded one
    can use `Coevolve` with `SSAStepper` for better performance. See the
    documentation, particularly the first and second tutorials, for details on
    defining and using bounded `VariableRateJump`s.
