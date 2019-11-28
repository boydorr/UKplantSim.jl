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

bsbi = loadtable("plantdata/PlantData_87to99.txt")
species = loadtable("plantdata/PlantData_Species.txt")
squares = loadtable("plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}

pa = readPlantATT("PLANTATT_19_Nov_08.csv")

spp_bsbi = unique(select(bsbi, :NAME))
spp_pa = select(pa, :Taxon_name)
cross_species = spp_bsbi ∩ spp_pa
bsbi = filter(b -> b.NAME ∈ cross_species, bsbi)

ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
bsbi = @transform bsbi {refval = UKclim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

bsbi = distribute(bsbi, 1)
start = startingArray(bsbi, length(species), 10)

abun = norm_sub_alpha(Metacommunity(start), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white,
grid = false, color = :algae, aspect_ratio = 1)


function coarsenRef(refval::Int64, width::Int64, sf::Int64)
    ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
    x, y = convert_coords(refval, size(ref.array, 1))
    xs = collect(x:(x + sf -1)); ys =  collect(y:(y + sf-1))
    xs = xs[xs .< 700]; ys = ys[ys .< 1250]
    return ref.array[xs, ys][1:end]
end

@everywhere namean(x) = mean(x[.!isnan.(x)])
@everywhere nastd(x) = std(x[.!isnan.(x)])
@everywhere using Statistics
dir = "HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
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
bsbi = @transform bsbi {tas = mean(meantas2015[coarsenRef(:refval, 700, 10)])}

# Calculate averages per species and plot as histogram
bsbi_counts = collect(@groupby bsbi :species {tas = namean(:tas), rainfall = namean(:rainfall), sun = namean(:sun), tas_st = nastd(:tas), rain_st = nastd(:rainfall)})
save(had_counts, "GBIF_had_prefs_UK")
