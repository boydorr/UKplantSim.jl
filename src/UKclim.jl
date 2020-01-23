module UKclim

include("extractHadUK.jl")
export retrieve_HadUK

include("ClimateTypes.jl")
export HadUK, LandCover, CropCover

include("Read.jl")
export readHadUK, readLC, readCHESS, readUKCP, readPlantATT, readNPMS, readCrop

include("Tools.jl")
export OSGR_eastnorth, extractvalues, createRef, LC2015cats, CC2017cats, startingArray, combineLC

include("abioticEnvUK.jl")
export lcAE, hadAE

include("traitsUK.jl")
export LCmatch, LCtrait

end
