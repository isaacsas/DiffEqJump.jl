# JumpProcesses.jl API

```@meta
CurrentModule = JumpProcesses
```

## Core Types

```@docs
ExtendedJumpArray
JumpProblem
PureLeaping
SSAStepper
SplitCoupledJumpProblem
reset_aggregated_jumps!
```

## Jump Types

```@docs
ConstantRateJump
MassActionJump
VariableRateJump
RegularJump
JumpSet
RateBounds
```

## Aggregator Types

Aggregators are the underlying algorithms used for sampling
[`ConstantRateJump`](@ref)s, [`MassActionJump`](@ref)s, and
[`VariableRateJump`](@ref)s.

```@docs
BracketData
CCNRM
Coevolve
Direct
DirectCR
DirectCRDirect
DirectFW
FRM
FRMFW
NRM
NSM
RDirect
RSSA
RSSACR
SortingDirect
get_num_majumps
needs_depgraph
needs_vartojumps_map
```

## Variable Rate Aggregators

```@docs
VariableRateAggregator
VR_Direct
VR_DirectFW
VR_FRM
```

## Tau-Leaping Algorithms

```@docs
EnsembleGPUKernel
SimpleAdaptiveTauLeaping
SimpleExplicitTauLeaping
SimpleImplicitTauLeaping
SimpleTauLeaping
SimpleTrapezoidalLeaping
```

## Spatial Jump APIs

```@docs
CartesianGrid
CartesianGridRej
SpatialMassActionJump
neighbors
num_sites
outdegree
```

## Reexported SciML common interface

`using JumpProcesses` also brings in the parts of the SciML common interface needed to
build the problem a [`JumpProblem`](@ref) wraps, solve it, drive the integrator from a
jump's `affect!`, and inspect the result -- so they do not have to be imported
separately. These names are owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/); JumpProcesses only re-exports
them:

  - Problems: [`DiscreteProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/discrete_types/),
    [`ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/),
    [`SDEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/sde_types/),
    [`EnsembleProblem`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/),
    `remake`, `NullParameters`
  - Functions: `DiscreteFunction`, `ODEFunction`, `SDEFunction`
  - Solutions: [`ODESolution`](https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/),
    `EnsembleSolution`, `EnsembleSummary`, and the `EnsembleAnalysis` module
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`
  - Solving: `solve`, `solve!`, `init`, `step!`
  - Integrator interface: `add_tstop!`, `add_saveat!`, `savevalues!`,
    `set_proposed_dt!`, `set_t!`, `set_u!`, `reinit!`, `terminate!`, `u_modified!`,
    `derivative_discontinuity!`
  - Return status: `ReturnCode`, `successful_retcode`
  - [Callbacks](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/):
    `DiscreteCallback`, `ContinuousCallback`, `VectorContinuousCallback`, `CallbackSet`

`DiscreteProblem` and `EnsembleProblem` in particular are what most downstream code
reaches for through JumpProcesses -- see SciML/MomentClosure.jl#111 for what happens
when they are not re-exported.

Note that [`SSAStepper`](@ref) only supports `DiscreteCallback`s;
`ContinuousCallback` and `VectorContinuousCallback` are re-exported for use with the
ODE/SDE integrators a `JumpProblem` can be paired with.

Anything else from SciMLBase -- the BVP, DAE, DDE, nonlinear and optimization problem
classes, the SciML operators, and the internals -- is not re-exported here; import it
from SciMLBase directly.

## Random Number Generator Control

JumpProcesses supports controlling the random number generator (RNG) used for
jump sampling via the `rng` and `seed` keyword arguments to `solve` or `init`.

### `rng` keyword argument

Pass any `AbstractRNG` to `solve` or `init`:

```julia
using Random, StableRNGs

# Using a StableRNG for cross-version reproducibility
sol = solve(jprob, SSAStepper(); rng = StableRNG(1234))

# Using Julia's built-in Xoshiro
sol = solve(jprob, Tsit5(); rng = Xoshiro(42))
```

### `seed` keyword argument

As a shorthand, pass an integer `seed` to create a `Xoshiro` generator:

```julia
sol = solve(jprob, SSAStepper(); seed = 1234)
# equivalent to: solve(jprob, SSAStepper(); rng = Xoshiro(1234))
```

### Resolution priority

When both `rng` and `seed` are passed to the same `solve`/`init` call, `rng`
takes priority:

| User provides | Result |
|---|---|
| `rng` via `solve`/`init` | Uses that `rng` |
| `seed` via `solve`/`init` | Creates `Xoshiro(seed)` |
| Nothing | Uses `Random.default_rng()` (SSAStepper, ODE, tau-leaping) or a randomly-seeded `Xoshiro` (SDE) |

### Behavior by solver pathway

| Solver | Default RNG (nothing passed) | `rng` / `seed` support |
|---|---|---|
| `SSAStepper` | `Random.default_rng()` | Full support via `solve`/`init` kwargs |
| ODE solvers (e.g., `Tsit5`) | `Random.default_rng()` | Full support via `solve`/`init` kwargs |
| SDE solvers (e.g., `SRIW1`) | Randomly-seeded `Xoshiro` | Full support; `TaskLocalRNG` is auto-converted to `Xoshiro` |
| `SimpleTauLeaping` | `Random.default_rng()` | Full support via `solve` kwargs |

!!! note
    For reproducible simulations, always pass an explicit `rng` or `seed`.
    The default RNG is shared global state and may produce different results
    depending on prior usage.

# Private / Developer API

```@docs
SSAIntegrator
```

## Internal Dispatch Pathways

The following table documents which code handles `solve`/`init` for each solver
type. This is relevant for developers working on JumpProcesses or its solver
backends.

| Solver type | `__solve` handled by | `__init` handled by | Uses `__jump_init`? |
|---|---|---|---|
| `SSAStepper` | JumpProcesses (`solve.jl`) | JumpProcesses (`SSA_stepper.jl`) | No |
| ODE (e.g., `Tsit5`) | JumpProcesses (`solve.jl`) | JumpProcesses (`solve.jl`) → OrdinaryDiffEq | Yes |
| SDE (e.g., `SRIW1`) | StochasticDiffEq | StochasticDiffEq | No |
| `SimpleTauLeaping` | JumpProcesses (`simple_regular_solve.jl`, custom `DiffEqBase.solve`) | N/A | No |

For **SSAStepper**, `rng` is resolved via `resolve_rng` in `SSA_stepper.jl`'s
`__init` and stored on the [`SSAIntegrator`](@ref).

For **ODE solvers**, `rng` is resolved via `resolve_rng` in `__jump_init`
(`solve.jl`) and forwarded to OrdinaryDiffEq's `init`, which stores it on the
`ODEIntegrator`.

For **SDE solvers**, StochasticDiffEq handles the full solve/init pathway
directly (JumpProcesses' ambiguity-fix `__solve` method is never dispatched to).
StochasticDiffEq has its own `_resolve_rng` that additionally handles
`TaskLocalRNG` conversion and the problem's stored seed.

For **tau-leaping**, JumpProcesses defines a custom `DiffEqBase.solve` that
bypasses the standard `__solve`/`__init` pathway. It calls `resolve_rng`
directly with the `rng` and `seed` kwargs from the `solve` call.
