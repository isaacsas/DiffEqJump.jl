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

For SSAStepper, OrdinaryDiffEq ODE/DAE solvers, SimpleTauLeaping, and
SimpleExplicitTauLeaping, an explicit `rng` takes priority over `seed`:

| User provides | Result |
|---|---|
| `rng` via `solve`/`init` | Uses that `rng` |
| `seed` via `solve`/`init` | Creates `Xoshiro(seed)` |
| Nothing | Uses `Random.default_rng()` |

StochasticDiffEq handles SDE/RODE RNGs for both `solve` and `init`. An explicit
RNG other than `TaskLocalRNG` takes priority over seeds. With no RNG, or with
`TaskLocalRNG`, a nonzero solve/init `seed` takes priority, followed by the
underlying problem's nonzero `seed`; otherwise StochasticDiffEq creates a
randomly seeded `Xoshiro`. It converts `TaskLocalRNG` to `Xoshiro` rather than
storing it on the stochastic integrator. In this pathway, `seed=0` means no
seed override, whereas SSAStepper and the ODE/DAE pathway use `Xoshiro(0)`.

### Behavior by solver pathway

| Solver | Default RNG (nothing passed) | `rng` / `seed` support |
|---|---|---|
| `SSAStepper` | `Random.default_rng()` | Full support via `solve`/`init` kwargs |
| OrdinaryDiffEq ODE/DAE solvers (e.g., `Tsit5`, `DFBDF`) | `Random.default_rng()` | Full support via `solve`/`init` kwargs |
| StochasticDiffEq SDE/RODE solvers (e.g., `SRIW1`, `RandomEM`) | `Xoshiro` from the stored problem seed, or a random seed | Full support; `TaskLocalRNG` follows the seed policy above and is converted to `Xoshiro` |
| `SimpleTauLeaping` | `Random.default_rng()` | Full support via `solve` kwargs |

!!! note
    For reproducible simulations, always pass an explicit `rng` or `seed`.
    Julia's default RNG is task-local; results may depend on prior draws in the task.

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
| OrdinaryDiffEq ODE/DAE (e.g., `Tsit5`, `DFBDF`) | JumpProcesses (`solve.jl`) | JumpProcesses' OrdinaryDiffEqCore extension → OrdinaryDiffEq | Yes |
| StochasticDiffEq SDE/RODE (e.g., `SRIW1`, `RandomEM`) | StochasticDiffEqCore | StochasticDiffEqCore's JumpProblem initializer | No |
| StochasticDiffEq jump algorithms (e.g., `TauLeaping`) | StochasticDiffEqCore | StochasticDiffEqCore's specialized jump-algorithm initializer | No |
| `SimpleTauLeaping` | JumpProcesses (`simple_regular_solve.jl`, custom `DiffEqBase.solve`) | N/A | No |

For **SSAStepper**, `rng` is resolved via `resolve_rng` in `SSA_stepper.jl`'s
`__init` and stored on the [`SSAIntegrator`](@ref).

For **OrdinaryDiffEq ODE/DAE solvers**, `rng` is resolved via `resolve_rng` in `__jump_init`
(`solve.jl`) and forwarded to OrdinaryDiffEq's `init`, which stores it on the
`ODEIntegrator`.

For **StochasticDiffEq SDE/RODE solvers**, the backend's more specific `__solve`
and `__init` methods handle the original JumpProblem and unchanged RNG/seed
inputs. Both routes use the backend's JumpProblem initializer, including its RNG
policy, jump-state copying, and callback setup. DiffEqBase merges stored problem
keywords and user callbacks before dispatching `init`.

Control stochastic jump-state copying through the backend's alias specifier:
`alias = SciMLBase.SDEAliasSpecifier(; alias_jumps = false)` for SDE problems,
or `SciMLBase.RODEAliasSpecifier` for RODE problems. Set the `alias_jumps` field
to `true` to reuse the original jump state. When that field is unspecified, the
backend aliases on thread 1 and copies on other threads. Both `solve` and `init`
follow this policy; `alias_jumps` is a field of the specifier, not a standalone
keyword argument.

For **tau-leaping**, JumpProcesses defines a custom `DiffEqBase.solve` that
bypasses the standard `__solve`/`__init` pathway. It calls `resolve_rng`
directly with the `rng` and `seed` kwargs from the `solve` call.
