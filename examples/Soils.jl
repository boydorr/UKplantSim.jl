using UKclim
using Simulation
using MyUnitful
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using JuliaDB
using JuliaDBMeta
using Distributions
using JLD

# using Plots
# gr()
# soils = readSoils("HuttonSoils.tif")
# heatmap(soils.array[:, 5e5m .. 1.25e6m]', aspect_ratio = 1, grid = false,
#     background_color =:lightblue, background_color_outside = :white)


# Import HadUK grid data
dir = "HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

soils = readSoils("HuttonSoils.tif")

# Create abiotic environment
solar = uconvert.(kJ, 1km^2 .* sun.array[:, 5e5m .. 1.25e6m, :] .* 1000.0*(W/m^2))
sol = SolarTimeBudget(solar, 1)

rain = Array(rainfall.array[:, 5e5m .. 1.25e6m, :])
water = WaterTimeBudget(rain, 1)

bud = BudgetCollection2(sol, water)

active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))

temp = HadUK(tas.array[:, 5e5m .. 1.25e6m, :])
temp = hadAE(temp, sol, active)
rain = HadUK(rainfall.array[:, 5e5m .. 1.25e6m, :])
rain = hadAE(rain, sol, active)
#lcae = lcAE(lc, 1000.0kJ/km^2, 242_495km^2)
soils = soilAE(soils, sol, active)

hab = HabitatCollection3(temp.habitat, rain.habitat, soils.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, temp.active, bud, temp.names)

traits = JuliaDB.load("BSBI_had_prefs_UK")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)
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

soil_traits = JuliaDB.load("BSBI_soil_prefs_uk")
soil_traits = LCtrait(collect(select(soil_traits, :soil)))

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection3(temp_traits, prec_traits, soil_traits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)

rel1 = Gauss{eltype(ae.habitat.h1)}()
rel2 = Gauss{eltype(ae.habitat.h2)}()
rel3 = LCmatch{eltype(ae.habitat.h3)}()
rel = multiplicativeTR3(rel1, rel2, rel3)
eco = Ecosystem(emptypopulate!, sppl, ae, rel)

#start = startingArray1(bsbi, length(traits), 10)
start = JLD.load("StartArray2.jld", "start")
start = reshape(start, 1626, 700, 1250)[:, :, 501:1250]
eco.abundances.matrix .+= reshape(start, 1626, 525000)

simulate!(eco, 20years, 1month)
JLD.save("BSBI_soil.jld", "abun", eco.abundances.matrix)
