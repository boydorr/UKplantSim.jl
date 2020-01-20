using Simulation
import Simulation.AbstractTraitRelationship, AbstractTraits

mutable struct LCmatch{TR} <: AbstractTraitRelationship{TR}
end
function (::LCmatch{TR})(niche::TR, pref::Vector{TR}) where TR
  if niche in pref
    return 1.0
  else
    return 0.0
  end
end
iscontinuous(tr::LCmatch{TR}) where TR = false
function eltype(tr::LCmatch{TR}) where TR
    return TR
end


mutable struct LCtrait{D <: Number} <: AbstractTraits{D}
  vals::Array{Array{D, 1}, 1}
end

iscontinuous(trait::LCtrait{D}) where D = false
function eltype(trait::LCtrait{D}) where D
    return D
end

function LCtrait(vals::Array{Array{D, 1}, 1}) where D  <: AbstractFloat
    return LCtrait{typeof(1.0)}(vals)
end
