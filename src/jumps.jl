struct ConstantRateJump{F1,F2} <: AbstractJump
  rate::F1
  affect!::F2
end

# when to switch from tuples for rates/affects to FunctionWrappers
const TUPLE_TO_FWRAPPER_CUTOFF = 10

function get_jump_info_tuples(constant_jumps)
  rates    = ((c.rate for c in constant_jumps)...)
  affects! = ((c.affect! for c in constant_jumps)...)

  return rates, affects!
end

using FunctionWrappers
function get_jump_info_fwrappers(u, p, t, constant_jumps)
  RateWrapper   = FunctionWrappers.FunctionWrapper{typeof(t),Tuple{typeof(u), typeof(p), typeof(t)}}
  rates         = [RateWrapper(c.rate) for c in constant_jumps]
  AffectWrapper = FunctionWrappers.FunctionWrapper{Void,Tuple{Any}}
  affects!      = [AffectWrapper(x->(c.affect!(x);nothing)) for c in constant_jumps]

  return rates, affects!
end


struct VariableRateJump{R,F,I,T,T2} <: AbstractJump
  rate::R
  affect!::F
  idxs::I
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
  abstol::T
  reltol::T2
end

VariableRateJump(rate,affect!;
                   idxs = nothing,
                   rootfind=true,
                   save_positions=(true,true),
                   interp_points=10,
                   abstol=1e-12,reltol=0) = VariableRateJump(
                              rate,affect!,idxs,
                              rootfind,interp_points,
                              save_positions,abstol,reltol)

struct RegularJump{R,C,DC,MD}
    rate::R
    c::C
    dc::DC
    mark_dist::MD
    constant_c::Bool
end

RegularJump(rate,c,dc;mark_dist = nothing,constant_c = false) =
            RegularJump(rate,c,dc,mark_dist,constant_c)

struct JumpSet{T1,T2,T3} <: AbstractJump
  variable_jumps::T1
  constant_jumps::T2
  regular_jump::T3
end

JumpSet(jump::ConstantRateJump) = JumpSet((),(jump,),nothing)
JumpSet(jump::VariableRateJump) = JumpSet((jump,),(),nothing)
JumpSet(jump::RegularJump) = JumpSet((),(),jump)
JumpSet() = JumpSet((),(),nothing)
JumpSet(jb::Void) = JumpSet()

# For Varargs, use recursion to make it type-stable

JumpSet(jumps::AbstractJump...) = JumpSet(split_jumps((), (), nothing, jumps...)...)

@inline split_jumps(vj, cj, rj) = vj, cj, rj
@inline split_jumps(vj, cj, rj, v::VariableRateJump, args...) = split_jumps((vj..., v), cj, rj, args...)
@inline split_jumps(vj, cj, rj, c::ConstantRateJump, args...) = split_jumps(vj, (cj..., c), rj, args...)
@inline split_jumps(vj, cj, rj, c::RegularJump, args...) = split_jumps(vj, cj, regular_jump_combine(rj,c), args...)
@inline split_jumps(vj, cj, rj, j::JumpSet, args...) = split_jumps((vj...,j.variable_jumps...), (cj..., j.constant_jumps...), regular_jump_combine(rj,j.regular_jump), args...)

regular_jump_combine(rj1::RegularJump,rj2::Void) = rj1
regular_jump_combine(rj1::Void,rj2::RegularJump) = rj2
regular_jump_combine(rj1::Void,rj2::Void) = rj1
regular_jump_combine(rj1::RegularJump,rj2::RegularJump) = error("Only one regular jump is allowed in a JumpSet")
