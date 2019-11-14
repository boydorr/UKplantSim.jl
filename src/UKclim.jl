module UKclim

include("extractHadUK.jl")
export retrieve_HadUK

include("ClimateTypes.jl")
export HadUK

include("Read.jl")
export readHadUK, readLC, readPlantATT, readNPMS

include("Tools.jl")
export OSGR_eastnorth, extractvalues, createRef, LC2015cats

include("abioticEnvUK.jl")
export lcAE, hadAE

end
