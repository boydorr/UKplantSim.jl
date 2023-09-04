using Pkg; Pkg.instantiate()
using UKplantSim
using DataFrames
using DataFramesMeta
using BritishNationalGrid
using Unitful
using UKplantSim.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using EcoSISTEM
using Distributions
using Diversity
using JLD2
using Plots
using DataPipeline

handle = DataPipeline.initialise()

path = link_read!(handle, "UKModel/Had-prefs")
traits = JLD2.load_object(path)
removed_rows = findall(isnan.(traits.sun) .| isnan.(traits.rainfall) .| isnan.(traits.tas_st) .| isnan.(traits.rain_st) .| (traits.rain_st .<= 0) .| (traits.tas_st .<= 0))
keep_rows = setdiff(1:nrow(traits), removed_rows)
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)
numSpecies = nrow(traits)
individuals = Int(1e9)

# Set up species requirements
sizes = abs.(rand(Normal(0.01, 0.01), numSpecies))
solarreq = traits.sun .* kJ
req1 = SolarRequirement(solarreq .* sizes)

waterreq = traits.rainfall .* mm
req2 = WaterRequirement(waterreq .* sizes)

req = ReqCollection2(req1, req2)

tmean = traits.tas .* K
tsd = traits.tas_st .* K
tsd .+= 1e-3K
temp_traits = GaussTrait(tmean, tsd)

pmean = traits.rainfall .* mm
psd = traits.rain_st .* mm
psd .+= 1e-3mm
prec_traits = GaussTrait(pmean, psd)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection2(temp_traits, prec_traits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)


# Read in landcover 2015 data
file = link_read!(handle, "UKModel/LCM")
lc = readLC(file)

# Import HadUK grid data
dir = link_read!(handle, "UKModel/HadUK-tas")
times = collect(2021year:1month:2021year+11month)
tas = readHadUK(dir, "tas", times)
dir = link_read!(handle, "UKModel/HadUK-rain")
rainfall = readHadUK(dir, "rainfall", times)
dir = link_read!(handle, "UKModel/HadUK-sun")
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

path = link_read!(handle, "UKModel/StartArray")
start = JLD2.load_object(path)
eco.abundances.matrix .+= start[keep_rows, :]

abun = norm_sub_alpha(Metacommunity(start), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
# heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white,
grid = false, color = :algae, aspect_ratio = 1)

simulate!(eco, 1year, 1month)

abun = norm_sub_alpha(Metacommunity(abun), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
file = link_write!(handle, "GBIF-uk-run")
Plots.pdf(file)

DataPipeline.finalise(handle)
