module UKclim

include("extractHadUK.jl")
export retrieve_HadUK

include("ClimateTypes.jl")
export HadUK, CropCover

include("Read.jl")
export readHadUK, readLC, readPlantATT, readNPMS, readCrop

include("Tools.jl")
export OSGR_eastnorth, extractvalues, createRef, LC2015cats, CC2017cats, startingArray

include("abioticEnvUK.jl")
export lcAE, hadAE

end
