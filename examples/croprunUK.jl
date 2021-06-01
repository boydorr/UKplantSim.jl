cd("/home/claireh/Documents/UK")
using UKclim
using JuliaDB
using JuliaDBMeta
using BritishNationalGrid
using Unitful
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using EcoSISTEM
using Distributions
using Diversity
using JLD
using Plots
gr()

##############
## BSBI RUN ##
##############

bsbi = loadtable("plantdata/PlantData_87to99.txt")
species = loadtable("plantdata/PlantData_Species.txt")
squares = loadtable("plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
bsbi = @transform bsbi {refval = UKclim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

traits = JuliaDB.load("BSBI_had_prefs_UK")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)
bsbi = filter(u-> u.NAME in JuliaDB.select(traits, :NAME), bsbi)
numSpecies = length(traits)
individuals = Int(1e9)

# Set up species requirements
sizes = abs.(rand(Normal(1.0, 0.1), numSpecies)) .* m^2
solarreq = collect(select(traits, :sun)) .* (kJ/km^2)
req1 = SolarRequirement(uconvert.(kJ, solarreq .* sizes))

waterreq = collect(select(traits, :rainfall)) .* (mm/km^2)
req2 = WaterRequirement(uconvert.(mm, waterreq .* sizes))

req = ReqCollection2(req1, req2)

tmean = collect(select(traits, :tas)) .* K
tsd = collect(select(traits, :tas_st)) .* K
tsd .+= 1e-3K
temp_traits = GaussTrait(tmean, tsd)

pmean = collect(select(traits, :rainfall)) .* mm
psd = collect(select(traits, :rain_st)) .* mm
psd .+= 1e-3mm
prec_traits = GaussTrait(pmean, psd)

lctrait = Array{Vector{Float64}, 1}(undef, numSpecies)
fill!(lctrait, collect(1:10))
lctraits = LCtrait(lctrait)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection3(temp_traits, prec_traits, lctraits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)


# Import HadUK grid data
dir = "HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Create abiotic environment
solar = uconvert.(kJ, 1km^2 .* sun.array .* 1000.0*(W/m^2))
sol = SolarTimeBudget(solar, 1)

rain = Array(rainfall.array)
water = WaterTimeBudget(rain, 1)

bud = BudgetCollection2(sol, water)

active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))
temp = hadAE(tas, sol, active)
rain = hadAE(rainfall, sol, active)

lc = readLC("CEH_landcover_2015.tif")
lcae = lcAE(lc, sol, active)

hab = HabitatCollection3(temp.habitat, rain.habitat, lcae.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, temp.active, bud, temp.names)

rel1 = Gauss{eltype(ae.habitat.h1)}()
rel2 = Gauss{eltype(ae.habitat.h2)}()
rel3 = LCmatch{eltype(ae.habitat.h3)}()
rel = multiplicativeTR3(rel1, rel2, rel3)
eco = Ecosystem(emptypopulate!, sppl, ae, rel)

#start = startingArray1(bsbi, length(traits), 10)
start = JLD.load("StartArray2.jld", "start")
eco.abundances.matrix .+= start

simulate!(eco, 20years, 1month)
JLD.save("BSBI_lc.jld", "abun")


# PLot results
cd("/home/claireh/Documents/UK")
using Diversity
using UKclim
using EcoSISTEM.Units
using JLD
using Plots
gr()
times = collect(2008year:1month:2017year+11month)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)
active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))

abun = JLD.load("BSBI_lc.jld", "abun")
abun = norm_sub_alpha(Metacommunity(abun), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("BSBI_lc.pdf")

bsbi = loadtable("plantdata/PlantData_87to99.txt")
species = loadtable("plantdata/PlantData_Species.txt")
squares = loadtable("plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
bsbi = @transform bsbi {refval = UKclim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

traits = JuliaDB.load("Crop_had_prefs_UK")
bsbi = filter(u-> u.NAME in JuliaDB.select(traits, :NAME), bsbi)
numSpecies = length(traits)
individuals = Int(1e9)

# Set up species requirements
sizes = abs.(rand(Normal(1.0, 0.1), numSpecies)) .* m^2
solarreq = collect(select(traits, :sun)) .* (kJ/km^2)
req1 = SolarRequirement(uconvert.(kJ, solarreq .* sizes))

waterreq = collect(select(traits, :rainfall)) .* (mm/km^2)
req2 = WaterRequirement(uconvert.(mm, waterreq .* sizes))

req = ReqCollection2(req1, req2)

tmean = collect(select(traits, :tas)) .* K
tsd = collect(select(traits, :tas_st)) .* K
tsd .+= 1e-3K
temp_traits = GaussTrait(tmean, tsd)

pmean = collect(select(traits, :rainfall)) .* mm
psd = collect(select(traits, :rain_st)) .* mm
psd .+= 1e-3mm
prec_traits = GaussTrait(pmean, psd)


lctrait = Array{Vector{Float64}, 1}(undef, numSpecies)
crop_names = ["Beta vulgaris", "Vicia faba", "Grass", "Zea mays", "Brassica napus", "Other", "Solanum tuberosum", "Hordeum vulgare", "Triticum aestivum", "Hordeum vulgare", "Triticum aestivum"]
crop_prefs = [22.0, 23.0, 25.0, 26.0, 28.0, [29.0, 31.0], [30.0, 32.0], [29.0, 31.0], [30.0, 32.0]]
cropids = vcat([findall(x .== select(traits, :NAME)) for x in crop_names]...)
fill!(lctrait, [1:2; 4:10])
map(cropids, crop_prefs) do id, p
    lctrait[id] = [p; 3.0; 27.0]
end
lctraits = LCtrait(lctrait)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection3(temp_traits, prec_traits, lctraits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)


# Import HadUK grid data
dir = "HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Create abiotic environment
solar = uconvert.(kJ, 1km^2 .* sun.array .* 1000.0*(W/m^2))
sol = SolarTimeBudget(solar, 1)

rain = Array(rainfall.array)
water = WaterTimeBudget(rain, 1)

bud = BudgetCollection2(sol, water)

active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))
temp = hadAE(tas, sol, active)
rain = hadAE(rainfall, sol, active)

lc = readLC("CEH_landcover_2015.tif")
crop = readCrop("Crop2017.tif")

newlc = combineLC(lc, crop)
newlc.array[newlc.array .== 24.0] /= 6.0
lcae = lcAE(newlc, sol, active)

hab = HabitatCollection3(temp.habitat, rain.habitat, lcae.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, temp.active, bud, temp.names)

rel1 = Gauss{eltype(ae.habitat.h1)}()
rel2 = Gauss{eltype(ae.habitat.h2)}()
rel3 = LCmatch{eltype(ae.habitat.h3)}()
rel = multiplicativeTR3(rel1, rel2, rel3)
eco = Ecosystem(emptypopulate!, sppl, ae, rel)

# start = startingArray1(bsbi, length(traits), 10)
# JLD.save("StartArrayCrop.jld", "start", start)

start = JLD.load("StartArrayCrop.jld", "start")
start = reshape(start, size(start, 1), 700, 1250)
ag = findall((newlc.array .== 3) .| (newlc.array .> 21))
start[:, ag] .= 0
for i in eachindex(cropids)
    locs = vcat([findall(newlc.array .== x) for x in lctrait[cropids[i]]]...)
    start[cropids[i], locs] .= 1e3
end
eco.abundances.matrix .+= reshape(start, size(start, 1), 700 * 1250)

simulate!(eco, 20years, 1month)
JLD.save("BSBI_abun_crop2.jld", "abun", eco.abundances.matrix)

# PLot results
cd("/home/claireh/Documents/UK")
using Diversity
using UKclim
using EcoSISTEM.Units
using JLD
using Plots
gr()
times = collect(2008year:1month:2017year+11month)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)
active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))

abun = JLD.load("BSBI_abun_crop2.jld", "abun")
abun = norm_sub_alpha(Metacommunity(abun), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("BSBI_crop.pdf")
