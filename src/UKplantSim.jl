module UKplantSim

#include("extractHadUK.jl")
#export retrieve_HadUK

using Unitful
using Unitful.DefaultSymbols

include("units.jl")

include("ClimateTypes.jl")
export HadUK, LandCover, CropCover, Soils

include("Read.jl")
export readHadUK, readLC, readCHESS, readUKCP, readPlantATT, readNPMS, readCrop, readSoils

include("Tools.jl")
export OSGR_eastnorth, extractvalues, createRef, LC2015cats, CC2017cats,
    startingArray, combineLC, coarsenRef

include("abioticEnvUK.jl")
export lcAE, hadAE, soilAE

include("traitsUK.jl")
export LCmatch, LCtrait

end
