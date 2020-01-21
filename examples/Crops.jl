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

import UKclim.LandCover

cd("/home/claireh/Documents/UK/")

crop = readCrop("Crop2017.tif")

times = collect(2008year:1month:2017year+11month)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)
active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))

cc = Array(crop.array[:, 1000m .. 1.25e6m])
cc[isnan.(cc)] .= 0
cc[.!active] .= NaN
crop = CropCover(AxisArray(cc, Axis{:northing}(sun.array.axes[1].val), Axis{:easting}(sun.array.axes[2].val)))
plot(crop)
Plots.pdf("CropCover2017.pdf")


lc = readLC("CEH_landcover_2015.tif")
plot(lc)
Plots.pdf("LandCover2017.pdf")

function combineLC(lc::LandCover, cc::CropCover)
    lc = lc.array[:, 0m .. 1.25e6m]
    cropland = findall(cc.array .> 0)
    agland = findall(lc .== 3)
    inter = agland ∩ cropland
    lc[inter] += (cc.array[inter] .+ 18)
    return LandCover(lc)
end

newlc = combineLC(lc, crop)
plot(newlc)
Plots.pdf("LandCropCombo.pdf")

start = JLD.load("StartArray.jld", "start")
start = reshape(start, size(start, 1), 700, 1250)
ag = findall(lc.array .== 3)
start[:, ag] .= 0
crops = zeros(11, 700, 1250)
for i in 1:11
    locs = findall((newlc.array .== 3) .| (newlc.array .== (i+21)))
    crops[i, locs] .= 1e3
end
start = cat(start, crops, dims = 1)
abun = start[end, :, :]
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("StartAbun.pdf")

abun = reshape(start, size(start, 1), 700 * 1250)
abun = norm_sub_alpha(Metacommunity(abun), 1.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("StartCropDiv.pdf")

using JLD
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
namean(x) = mean(x[.!isnan.(x)])
nastd(x) = std(x[.!isnan.(x)])

trt_means = map(1:11) do i
    locs = findall(newlc.array .== (i+21))
    return [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
end
trt_means = hcat(trt_means...)
JLD.save("Crop_trait_means.jld", "trt_means", trt_means)

trt_stds = map(1:11) do i
    locs = findall(newlc.array .== (i+21))
    return [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]
end
trt_stds = hcat(trt_stds...)
JLD.save("Crop_trait_stds.jld", "trt_stds", trt_stds)
