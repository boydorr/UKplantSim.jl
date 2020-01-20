using Simulation
import Simulation: AbstractTraitRelationship, AbstractTraits, iscontinuous, _traitfun, getpref, combineTR
import Base.eltype

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

iscontinuous(tr::multiplicativeTR3{TR1, TR2, TR3} ) where {TR1, TR2, TR3} =
    [iscontinuous(tr.tr1), iscontinuous(tr.tr2), iscontinuous(tr.tr3)]


function _traitfun(hab::HabitatCollection3, trts::TraitCollection3, rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
    res1 = _traitfun(hab.h1, trts.t1, rel.tr1, pos, spp)
    res2 = _traitfun(hab.h2, trts.t2, rel.tr2, pos, spp)
    res3 = _traitfun(hab.h3, trts.t3, rel.tr3, pos, spp)
    return combineTR(rel)(res1, res2, res3)
end

function _traitfun(hab::DiscreteHab, trts::LCtrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        vals = getpref(trts, spp)
    return rel(h, vals)
end

function getpref(traits::LCtrait, spp::Int64)
  return traits.vals[spp]
end
