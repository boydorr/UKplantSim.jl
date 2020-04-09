cd("/home/claireh/Documents/UK")
using UKclim
using JuliaDB
using JuliaDBMeta
using BritishNationalGrid
using Unitful
using Simulation.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using StatsBase
using DataStructures
using Simulation
using Distributions
using Diversity
using JLD
using Plots
gr()

bsbi = loadtable("plantdata/PlantData_87to99.txt")
species = loadtable("plantdata/PlantData_Species.txt")
squares = loadtable("plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.3e6m)
bsbi = @transform bsbi {refval = UKclim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

file = "CEH_landcover_2015.tif"
lc = readLC(file)
bsbi = @transform bsbi {lc = lc.array[:refval]}

function top_ten(x)
  count = collect(values(countmap(x)))
  species = collect(keys(countmap(x)))
  thresh = 0.9 * maximum(count)
  top_species = species[sortperm(count, rev=true)]
  return top_species[count[sortperm(count, rev=true)] .> thresh]
end
bsbi_lc = @groupby bsbi :lc {species = top_ten(:NAME)}

broadleavedwoods = filter(b -> b.lc == 1, bsbi)
filter(b -> b.NAME == "Quercus robur", broadleavedwoods)
spp = select(broadleavedwoods, :species)[1]
count = select(broadleavedwoods, :counts)[1]
spp[sortperm(count, rev=true)]
