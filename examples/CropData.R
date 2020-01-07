library(rgdal)
setwd("~/Downloads/")
fgdb <- "~/Downloads/Download_Crops2017_1399393/lccm-2017_3308728/lccm-2017_3308728.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)
fc <- readOGR(dsn=fgdb,layer="crop_map_2017")
fc@data
mode_func <- function(x, na.rm){
  as.factor(names(which.max(table(x))))
}

library(tidyverse)
fc_agg = aggregate(fc, by= "crop_code", dissolve = T)

library(raster)
r <- raster(extent(fc), res= c(1000,1000))
rr <- rasterize(fc, r, field = "crop_code", fun = mode_func)
