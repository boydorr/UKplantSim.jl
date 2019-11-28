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
using Distributions
using Diversity
using Plots
gr()

uk = JuliaDB.load("UKspecies")
uk = @transform uk {east = BNGPoint(lon = :decimallongitude, lat = :decimallatitude).e, north = BNGPoint(lon = :decimallongitude, lat = :decimallatitude).n}

# Create reference for UK grid
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
uk = @transform uk {refval = UKclim.extractvalues(:east * m, :north * m, ref)}

traits = load("GBIF_had_prefs_UK")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
uk = filter(u-> u.species in select(traits, :species), uk)
numSpecies = length(traits)
individuals = Int(1e9)

# Set up species requirements
sizes = abs.(rand(Normal(0.01, 0.01), numSpecies))
solarreq = collect(select(traits, :sun)) .* kJ
req1 = SolarRequirement(solarreq .* sizes)

waterreq = collect(select(traits, :rainfall)) .* mm
req2 = WaterRequirement(waterreq .* sizes)

req = ReqCollection2(req1, req2)

tmean = collect(select(traits, :tas)) .* K
tsd = collect(select(traits, :tas_st)) .* K
tsd .+= 1e-3K
temp_traits = GaussTrait(tmean, tsd)

pmean = collect(select(traits, :rainfall)) .* mm
psd = collect(select(traits, :rain_st)) .* mm
psd .+= 1e-3mm
prec_traits = GaussTrait(pmean, psd)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection2(temp_traits, prec_traits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)


# Read in landcover 2015 data
file = "CEH_landcover_2015.tif"
lc = readLC(file)
plot(lc)

# Import HadUK grid data
dir = "HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
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
#lcae = lcAE(lc, 1000.0kJ/km^2, 242_495km^2)

hab = HabitatCollection2(temp.habitat, rain.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, temp.active, bud, temp.names)

rel1 = Gauss{eltype(ae.habitat.h1)}()
rel2 = Gauss{eltype(ae.habitat.h2)}()
rel = multiplicativeTR2(rel1, rel2)
eco = Ecosystem(emptypopulate!, sppl, ae, rel)

start = startingArray(uk, numSpecies)
eco.abundances.matrix .+= start

abun = norm_sub_alpha(Metacommunity(start), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white,
grid = false, color = :algae, aspect_ratio = 1)

simulate!(eco, 1year, 1month)

abun = norm_sub_alpha(Metacommunity(abun), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
