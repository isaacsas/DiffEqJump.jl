module JumpProcessesOrdinaryDiffEqCoreExt

using JumpProcesses
import DiffEqBase
import SciMLBase
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, DAEAlgorithm,
    StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm

function _jump_init(_jump_prob, alg; merge_callbacks = true, kwargs...)
    kwargs = DiffEqBase.merge_problem_kwargs(_jump_prob; merge_callbacks, kwargs...)
    JumpProcesses.__jump_init(_jump_prob, alg; kwargs...)
end

function SciMLBase.__init(
        _jump_prob::JumpProcesses.JumpProblem{IIP, P},
        alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm,
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        kwargs...) where {IIP, P}
    # Cover the full intersection with OrdinaryDiffEqCore's initializer. The
    # narrower stochastic method below preserves its backend-owned setup.
    _jump_init(_jump_prob, alg; kwargs...)
end

function SciMLBase.__init(
        _jump_prob::JumpProcesses.JumpProblem{IIP, P},
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        alias_jump = nothing, alias = nothing, kwargs...) where {IIP, P}
    if alias_jump !== nothing
        # Preserve the legacy override previously consumed by __jump_init. The
        # backend already supports alias_jumps; retain its other alias settings.
        aliases = if alias isa Union{SciMLBase.SDEAliasSpecifier, SciMLBase.RODEAliasSpecifier}
            alias
        else
            constructor = alg isa StochasticDiffEqRODEAlgorithm ?
                SciMLBase.RODEAliasSpecifier : SciMLBase.SDEAliasSpecifier
            alias isa Bool ? constructor(; alias) : constructor()
        end
        fields = NamedTuple{fieldnames(typeof(aliases))}(
            ntuple(i -> getfield(aliases, i), fieldcount(typeof(aliases))))
        alias = typeof(aliases)(; fields..., alias_jumps = alias_jump)
    end
    # Use the stochastic backend's JumpProblem initialization, as its solve path does.
    # Keep rng and seed unresolved so it handles TaskLocalRNG conversion, problem
    # seeds, and seed metadata consistently. Keeping the wrapper also preserves its
    # jump aliasing and callback setup. The broader problem signature selects the
    # backend method without recursively dispatching to this ambiguity fix.
    # DiffEqBase.init_call already merged stored kwargs and user callbacks.
    invoke(SciMLBase.__init,
        Tuple{Union{SciMLBase.AbstractRODEProblem, JumpProcesses.JumpProblem},
            Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}},
        _jump_prob, alg; alias, kwargs...)
end

end
