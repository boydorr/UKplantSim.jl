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
using Plots

bsbi = loadtable("plantdata/PlantData_87to99.txt")
species = loadtable("plantdata/PlantData_Species.txt")
squares = loadtable("plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}
bsbi = distribute(bsbi, 12)

pa = readPlantATT("PLANTATT_19_Nov_08.csv")

spp_bsbi = unique(collect(select(bsbi, :NAME)))
spp_pa = collect(select(pa, :Taxon_name))
cross_species = spp_bsbi ∩ spp_pa
bsbi = filter(b -> b.NAME ∈ cross_species, bsbi)

ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
bsbi = @transform bsbi {refval = UKclim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

start = startingArray(bsbi, length(species), 10)

abun = norm_sub_alpha(Metacommunity(start), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white,
grid = false, color = :algae, aspect_ratio = 1)


bsbi = @transform bsbi {cr = coarsenRef(:refval, 700, 10)}

@everywhere namean(x) = mean(x[.!isnan.(x)])
@everywhere nastd(x) = std(x[.!isnan.(x)])
@everywhere using Statistics
dir = "HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
bsbi = @transform bsbi {tas = mean(meantas2015[:cr]), rainfall = mean(meanrainfall2015[:cr]), sun = mean(meansun2015[:cr])}

# Calculate averages per species and plot as histogram
bsbi_counts = collect(@groupby bsbi :NAME {tas = namean(:tas), rainfall = namean(:rainfall), sun = namean(:sun), tas_st = nastd(:tas), rain_st = nastd(:rainfall)})
save(bsbi_counts, "BSBI_had_prefs_UK")

lc = readLC("CEH_landcover_2015.tif")
bsbi = @transform bsbi {lc = lc.array[:refval]}
LC_counts = collect(@groupby bsbi :NAME {lc = [sum(:lc .== i) for i in 1:21]})
save(LC_counts, "BSBI_lc_prefs_uk")

save(bsbi, "BSBI_prefs_UK")

soils = readSoils("HuttonSoils.tif")
bsbi = @transform bsbi {soil = soils.array[:refval]}
soil_counts = collect(@groupby bsbi :NAME {soil = unique(:soil)})
soil_counts = @where soil_counts !all(isnan.(:soil))
save(soil_counts, "BSBI_soil_prefs_uk")
