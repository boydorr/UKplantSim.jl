r = raster("Documents/PhD/Data/CEH/CEH_landcover_2015.tif")
r2 = raster("Documents/PhD/Data/CEH/CEH_landcover_2015_ireland.tif")

r = projectRaster(r, crs = crs(r2), res = c(1000, 1000))
r = extend(r, extent)
r2 = extend(r2, extent)
plot(r2)
plot(r, add = T)
merge(r, r2)
