module JumpProcessesStochasticDiffEqCoreExt

using JumpProcesses: JumpProblem
import SciMLBase
import OrdinaryDiffEqCore: StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm
import StochasticDiffEqCore

# Forward stochastic `init` of a `JumpProblem` to StochasticDiffEqCore's own initializer,
# `_sde_init`, which is exported for packages that need to bypass `__init` dispatch. On
# Core versions whose `__init` covers `Union{AbstractRODEProblem, JumpProblem}`, that
# method and the `OrdinaryDiffEqCore` extension's four-family method are each more
# specific in one argument, so without this method `init(jprob, sde_alg)` is ambiguous.
# Core versions with per-algorithm `JumpProblem` methods (OrdinaryDiffEq.jl#4614) are more
# specific than this union method and take precedence, and Core's jump-algorithm method
# (e.g. `TauLeaping`) is also more specific. Keep the signature distinct from Core's: do
# not split it into separate SDE and RODE methods.
#
# `init_call` has already merged the `JumpProblem`'s stored keywords, and `_sde_init`
# documents that they must not be merged again.
function SciMLBase.__init(prob::JumpProblem,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}; kwargs...)
    StochasticDiffEqCore._sde_init(prob, alg; kwargs...)
end

end
