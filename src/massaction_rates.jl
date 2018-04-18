###############################################################################
# Stochiometry for a given reaction is a vector of pairs mapping species id to
# stochiometric coefficient.
###############################################################################

# @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rateconst,
#                               stochmat::AbstractVector{Pair{S,V}}, 
#                               num_dep_specs=nothing)::typeof(rateconst) where {T,S,V}

@fastmath function evalrxrate(speciesvec::AbstractVector{T}, rxidx::S, 
                              majump::MassActionJump{U,V,W})::R where {T,S,R,U <: AbstractVector{R},V,W}
    val = one(T)

    rateconst = majump.scaled_rates[rxidx]
    stochmat  = majump.reactant_stoch[rxidx]

    @inbounds for specstoch in stochmat
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        @inbounds for k = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

    rateconst * val
end


# @inline @fastmath function executerx!(speciesvec::AbstractVector{T},
#                                       net_stoich::AbstractVector{Pair{S,V}}, 
#                                       num_dep_specs=nothing) where {T,S,V}

@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, rxidx::S, 
                                      majump::MassActionJump{U,V,W}) where {T,S,U,V,W}                                        

    net_stoich = majump.net_stoch[rxidx]                                      
    @inbounds for specstoch in net_stoich
        speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end


function scalerates!(unscaled_rates::AbstractVector{U}, stochmat::AbstractVector{V}, 
                     num_dep_specs_vec=nothing) where {U,S,T,W <: Pair{S,T}, V <: AbstractVector{W}}

    @inbounds for i in eachindex(unscaled_rates)
        coef = one(T)
        @inbounds for specstoch in stochmat[i]
            coef *= factorial(specstoch[2])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end

function scalerate(unscaled_rate::U, stochmat::AbstractVector{Pair{S,T}}, 
                   num_dep_specs=nothing) where {U <: Number, S, T}                   
    coef = one(T)
    @inbounds for specstoch in stochmat
        coef *= factorial(specstoch[2])
    end
    unscaled_rate /= coef
end


###############################################################################
# Stochiometry for a given reaction is a StaticArrays vector of pairs mapping 
# species id to stochiometric coefficient.
###############################################################################

using StaticArrays
# majump::MassActionJump{U,V,W})::R where {T,S,R,U<:AbstractVector{R},SA1<:SArray,SA2<:SArray,V<:AbstractVector{SA1},W<:AbstractVector{SA2}}
@fastmath function evalrxrate(speciesvec::AbstractVector{T}, rxidx::S, 
                              majump::MassActionJump{U,V,W})::R where {T,S,R,U <: AbstractVector{R},P1,P2,P<:Pair{P1,P2}, SA1 <:SArray{Tuple{3},P}, SA2 <:SArray{Tuple{6},P}, V <: AbstractVector{SA1}, W <: AbstractVector{SA2}}

    num_dep_specs = majump.num_dep_specs[rxidx]
    stochmat      = majump.reactant_stoch[rxidx]
    rateconst     = majump.scaled_rates[rxidx]

    val = one(T)
    @inbounds for idx = 1:num_dep_specs
        stoich  = stochmat[idx]
        specpop = speciesvec[stoich[1]]
        val    *= specpop
        @inbounds for k = 2:stoich[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

    rateconst * val
end

# @inline @fastmath function executerx!(speciesvec::AbstractVector{T},
#                                       net_stoich::SArray{Tuple{6},Pair{S,V}}, 
#                                       num_dep_specs) where {T,S,V}
#majump::MassActionJump{U,V,W})::R where {T,S,U,SA1<:SArray,SA2<:SArray,V<:AbstractVector{SA1},W<:AbstractVector{SA2}}
@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, rxidx::S,
                                      majump::MassActionJump{U,V,W}) where {T,S,U,P1,P2,P<:Pair{P1,P2}, SA1 <:SArray{Tuple{3},P}, SA2 <:SArray{Tuple{6},P}, V <: AbstractVector{SA1}, W <: AbstractVector{SA2}}
    
    num_chg_specs = majump.num_chg_specs[rxidx]
    net_stoich = majump.net_stoch[rxidx]

    @inbounds for idx = 1:num_chg_specs
        stoich = net_stoich[idx]
        speciesvec[stoich[1]] += stoich[2]
    end

    nothing
end

#                     num_dep_specs_vec::AbstractVector{W}) where {U,S,T,V <: SArray, W}
function scalerates!(unscaled_rates::AbstractVector{U}, stochmat::AbstractVector{V}, 
                     num_dep_specs_vec::AbstractVector{W}) where {U,S,T,R <: Pair{S,T}, V <: SArray{Tuple{3},R}, W}

    @inbounds for i in eachindex(unscaled_rates)
        coef = one(T)
        @inbounds for idx = 1:num_dep_specs_vec[i] 
            specstoch = stochmat[i][idx]
            coef     *= factorial(specstoch[2])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end

function scalerate(unscaled_rate::U, stochmat::SArray{Tuple{3},Pair{S,T}}, 
                   num_dep_specs) where {U <: Number, S, T}

    coef = one(T)
    @inbounds for idx = 1:num_dep_specs
        specstoch = stochmat[idx]
        coef     *= factorial(specstoch[2])
    end
    unscaled_rate /= coef
end
