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
using Plots
plotlyjs()

uk = JuliaDB.load("data/Full_GBIF_africa")

# Read in landcover 2015 data
file = "data/CEH_landcover_2015.tif"
lc = readLC(file)
plot(lc)

# Import HadUK grid data
dir = "data/HadUK/tas/10s/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Create abiotic environment
solar = uconvert.(kJ, 1km^2 .* sun.array .* 1000.0*(W/m^2))
sol = SolarTimeBudget(solar, 1)

active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))
temp = hadAE(tas, sol, active)
rain = hadAE(rainfall, sol, active)
lcae = lcAE(lc, 1000.0kJ/km^2, 242_495km^2)

hab = HabitatCollection2(temp.habitat, rain.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(sol)}(hab, temp.active, sol, temp.names)
