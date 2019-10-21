using UKclim
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using Plots
plotlyjs()

retrieve_HadUK("tas", 1980, 1990, "/Users/claireh/Documents/PhD/GIT/UKclim/data/HadUK/tas")

dir = "data/HadUK/tas/"
times = collect(1980year:1month:1990year+11month)
tas1980 = readHadUK(dir, "tas", times)

plot(tas1980, 1980year + March, clim = (270, 290))

file = "data/CEH_landcover_2015.tif"
lc = readLC(file)

file = "data/CEH_landcover_2015_ireland.tif"
lci = readLC(file, 1.8e5, 3.7e5, 3e5, 4.6e5)

plot(lci)
