library(rgdal)
library(raster)
library(lidR)

# load point cloud
las <- readLAS('G:/DopeData/Data/Lidar/Kall/Test/Raw/3dm_32_322_5593_1_nw.laz')

# Metadata
summary(las)
las_check(las)

## Point Density
# The point density is an important characteristic  
# that determines the max. resolution of the dem

# create raster of the point density per m²
point_density <- grid_density(las, res = 1)
plot(point_density)

# Histogram of the points per pixel
hist(d, main = "Distribution of Points per m2",
     xlab = "Point Count", ylab = "Frequency")
summary(d)

## Generate DEMs

# DTM (bare surface)
dtm <- grid_terrain(las, res = 1, knnidw(k = 10, p = 2), keep_lowest = TRUE)
plot(dtm)

# DSM (Surface + objects)
dsm <- grid_canopy(las, res = 1, dsmtin(max_edge = 0))

# Canopy Height Model
chm <- dsm - dtm
plot(chm)
summary(chm)
chm[chm<0] <- 0

# Hillshade
slope <- terrain(dtm, opt='slope')
aspect <- terrain(dtm, opt='aspect')
hs <- hillShade(slope, aspect, angle=45, direction=315)
plot(hs)

## Extract buildings

# Class 21 = buildings, 17 = bridges
unique(las@data$Classification)
las_buildings <- filter_poi(las, Classification %in%  c(21, 17))
lidR::plot(las_buildings) # opens in new window

## Tree detection
ttops <- find_trees(las, lmf(ws = 5))
x = plot(las)
add_treetops3d(x, ttops)
