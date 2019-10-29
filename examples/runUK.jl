using UKclim
using JuliaDB
using JuliaDBMeta
using BritishNationalGrid
using Unitful
using MyUnitful
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using Simulation

# Read in landcover 2015 data
file = "data/CEH_landcover_2015.tif"
lc = readLC(file)
plot(lc)

dir = "data/HadUK/tas/10s/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/HadUK/sun/"
sun = readHadUK(dir, "sun", times)


lc.array[:, 1000.0m .. 1.25e6m]
