module UKclim

include("extractHadUK.jl")
export retrieve_HadUK

include("ClimateTypes.jl")
export HadUK

include("Read.jl")
export readHadUK, readLC

end
