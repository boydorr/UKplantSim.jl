using UKclim
using Plots
plotlyjs()

soils = readSoils("data/HuttonSoils.tif")
heatmap(soils.array)
