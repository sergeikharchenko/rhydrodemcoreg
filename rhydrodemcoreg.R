# Loading of necessary packages
library(raster) ; library(rgdal) ; library(rgeos) ; library(sf) ; library(RSAGA) ; library(foreach) ; library(doParallel)

# set working directory
setwd("C:/DTM_coreg/1_Sat") # EDITABLE PARAMETER. YOU NEED TO TYPE HERE A WORKING FOLDER

time_start <- timestamp()

################# PROCESSING ORIGINAL DEMs #################
# read DEMs
(dem1 <- raster("mejd1_dem.tif")) # EDITABLE PARAMETER. YOU NEED TO TYPE HERE THE FIRST DEM's NAME
# check DEM's view
plot(dem1)
# check DEM's coordinate reference system
if (grep(pattern = "utm", x = crs(dem1)) != 1) {
  center_point <- apply(matrix(extent(dem1), 2, 2), 2, mean)
  center_point <- spTransform(SpatialPoints(coords = matrix(center_point, 1), proj4string = crs(dem1)), "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  utm_zone <- 31 + (coordinates(center_point)[1] %/% 6)
  dem1 <- projectRaster(from = dem1, crs = paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +units=m +no_defs "))
}
# trim "white" edges
dem1 <- trim(dem1)

(dem2 <- raster("mejd2_dem.tif")) # EDITABLE PARAMETER. YOU NEED TO TYPE HERE THE SECOND DEM's NAME
# check DEM's view
plot(dem2)
# check DEM's coordinate reference system
if (!all.equal(crs(dem1), crs(dem2))) {
  dem2 <- projectRaster(dem2, dem1)
}
# trim "white" edges
dem2 <- trim(dem2)

# generate common extent
ext1 <- extent(Reduce(extend,list(dem1, dem2)))

# generate raster template
dem_temp <- raster(res = 0.25, ext = ext1, crs = crs(dem1)) # EDITABLE PARAMETER. Choose resolution of final data

# Unify resolution
system.time(dem1 <- raster::resample(dem1, dem_temp))
system.time(dem1 <- terra::resample(dem1, dem_temp))
dem2 <- resample(dem2, dem_temp)
plot(dem1) # View DEM1
plot(dem2, add = T)  # View DEM2

# generate overlapping area polygon for planar shift/correction
DoD <- dem2 - dem1
DoD_border <- boundaries(x = DoD, type = "outer")
DoD_border[DoD_border == 1] <- NA
DoD_border <- rasterToPolygons(resample(DoD_border, raster(ext = extent(DoD_border), res = 2, crs = crs(DoD_border))), dissolve = T)
DoD_inside_buffer <- gBuffer(DoD_border, width = -100) # EDITABLE PARAMETER. Choose size of inside buffer within overlapping area (default = -100)

# masking two DEM's inside areas for estamation of planar shift
dem1cr <- mask(crop(dem1, DoD_inside_buffer), DoD_inside_buffer)
dem2cr <- resample(mask(crop(dem2, DoD_inside_buffer), DoD_inside_buffer), dem1cr)
par(mfrow = c(1,3))
plot(dem1cr, main = "DEM 1, overlap area")
plot(dem2cr, main = "DEM 2, overlap area")
plot(dem1cr - dem2cr, main = "DEM 1-2, overlap area")
par(mfrow = c(1,1))

# make expand grid of possible shift values (default ±20 cells from the center)
shift_size <- 20*res(dem1cr)[1]
expand1 <- expand.grid(seq(-shift_size,shift_size,res(dem1cr)[1]),seq(-shift_size,shift_size,res(dem1cr)[1]))
# parallel computation of vertical error with different shift along X and Y axes
cl1 <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1)
registerDoParallel(cl1)
err1 <- foreach (m = 1:nrow(expand1), .packages = c("raster"), .combine = "c") %dopar% {
  cellStats(shift(dem2cr, expand1[m,1], expand1[m,2]) - dem1cr, sd)
}
stopCluster(cl1)

# see raster model of vertical error
plot(rasterFromXYZ(cbind(expand1, err1)), main = "Error distribution with \n shifting of dependent DEM")
points(expand1[which.min(err1),], pch = 19)
# write raster of vertical error
writeRaster(rasterFromXYZ(cbind(expand1, err1)), "error.tif", "GTiff", overwrite=TRUE)
dem2_corr <- shift(dem2, expand1[which.min(err1),1], expand1[which.min(err1),2])
dem2_corr <- resample(dem2_corr, dem_temp)
writeRaster(dem2_corr, "dem2_corr.tif", "GTiff", overwrite = T)

# create border lines
bdem1 <- boundaries(x = dem1, type = "inner")
bdem1[bdem1 == 0] <- NA
bdem2 <- boundaries(x = dem2, type = "inner")
bdem2[bdem2 == 0] <- NA

# create cross-points
bdem3 <- bdem1 & bdem2
bdem3 <- rasterToPoints(bdem3)
(ps <- SpatialPoints(coords = bdem3, crs(dem1)))
plot(dem1)
plot(dem2, add = T)
points(ps, pch = 19)

# OR choose most distant pair of points
n1 <- ifelse(which.max(as.matrix(dist(coordinates(ps)))) %% length(ps) == 0, length(ps), which.max(as.matrix(dist(coordinates(ps)))) %% length(ps))
n2 <- 1 + (which.max(as.matrix(dist(coordinates(ps)))) %/% length(ps) )

# OR, if it is necessary to choose not most distant points - use this code
own_points <- click(n = 2) # take points from map
n1 <- which.min(sapply(X = 1:length(ps), FUN = function(x) dist(rbind(own_points[1,], coordinates(ps)[x,1:2]))))
n2 <- which.min(sapply(X = 1:length(ps), FUN = function(x) dist(rbind(own_points[2,], coordinates(ps)[x,1:2]))))

# cut point set to 2-points
ps <- ps[c(n1,n2),]
points(ps, pch = 19, col = "red")
ps$z <- 1:length(ps)
# write Spatial Points layer
writeOGR(ps, "ps.shp", "ps.shp", "ESRI Shapefile")

# extract differences of elevation in endpoints
(zdiff <- extract(dem1, ps) - extract(dem2, ps))
ps$z <- zdiff

# create "shoulder" / basis line and third point
coord <- coordinates(ps)[,1:2]
mean_point <- apply(coord, 2, mean)
angle1 <- pi - ((pi/2) - (atan(apply(coord, 2, diff)[2] / apply(coord, 2, diff)[1])))
shoulder <- 200 # EDITABLE PARAMETER. Length of line (default - 200m)
new_point <- unlist(data.frame(x = mean_point[1] + shoulder * cos(angle1), y = mean_point[2] + shoulder * sin(angle1))) 

# correcting elevation along line
data1 <- data.frame(zdiff = c(zdiff, mean(zdiff), mean(zdiff)), x = c(coordinates(ps)[,1], mean_point[1], new_point[1]), y = c(coordinates(ps)[,2], mean_point[2], new_point[2]))
# view original endpoints
plot(ps, pch = 19)
# view new point for vertical error linear trend construction
points(SpatialPoints(coords = data1[,2:3]), pch = 19, col = "red")

# create 1st linear model of vertical error distribution
lm1 <- lm(zdiff ~ ., data1)
# convert raster to separate points
dem2p <- rasterToPoints(dem2)
# predict vertical error for all points
lm1p <- predict(lm1, data.frame(dem2p))
# create template for corrected DEM2
v_error <- dem2
v_error[!is.na(dem2)] <- lm1p
# view vertical error model
plot(v_error)
# create corrected DEM2
dem2new <- dem2 + v_error
dem2new_backup <- dem2new
# write corrected DEM2
writeRaster(dem2new, "dem2new.tif", "GTiff", overwrite = TRUE)

# extract old and new elevations. Check equality
extract(dem1, ps)
extract(dem2new, ps)

################# TILTING THE DEMs WITH DIFFERENT SLOPES #################
# tilt slave DEM (DEM2)
tilt <- function(i = 0, w = 1) {
  data2 <- data.frame(zdiff = c(0,0,i), x = c(coordinates(ps)[,1], new_point[1]), y = c(coordinates(ps)[,2], new_point[2]))
  
  # create 2nd linear model of vertical error distribution
  lm2 <- lm(zdiff ~ ., data2)
  # predict vertical error for all points
  lm2p <- predict(lm2, data.frame(dem2p))
  # create template for corrected DEM2
  v_error[!is.na(dem2)] <- lm2p
  # create corrected DEM2
  dem2new <- dem2new_backup + v_error
  # resample DEM2
  system.time(dem2new <- resample(dem2new, dem1))
  
  # create DEM of difference
  (dem_diff <- dem2new - dem1)
  # create "zero error" line
  con1 <- rasterToContour(dem_diff, levels = 0)
  con1 <- disaggregate(con1)
  con1 <- con1[which.max(SpatialLinesLengths(con1)),]
  
  # view data. Check "zero error" line are whole inside of overlap polygon
  plot(dem2new, main = paste0("Slope is ",100*i / shoulder," %"))
  plot(dem1, add = T)
  plot(DoD_border, add = T)
  lines(con1)
  points(ps, pch = 19)
  
  # create vector DEMs borders. Border 2
  dem2_border <- dem2new
  dem2_border[!is.na(dem2_border)] <- 1
  writeRaster(dem2_border, "dem2_border.tif", "GTiff", overwrite=TRUE)
  
  # using SAGA GIS for convert raster to vector polygon
  rsaga.geoprocessor("shapes_grid", 6, param = list(
    GRID = "dem2_border.tif",
    POLYGONS = "dem2_border",
    CLASS_ID = 1
  ))
  
  # read polygonal border 2
  polyg1 <- (readOGR("dem2_border.shp"))
  
  # create vector DEMs borders. Border 1
  dem1_border <- dem1
  dem1_border[!is.na(dem1_border)] <- 1
  # plot(dem5res)
  writeRaster(dem1_border, "dem1_border.tif", "GTiff", overwrite=TRUE)
  
  # using SAGA GIS for convert raster to vector polygon
  rsaga.geoprocessor("shapes_grid", 6, param = list(
    GRID = "dem1_border.tif",
    POLYGONS = "dem1_border",
    CLASS_ID = 1
  ))
  
  # read polygonal border 1
  polyg2 <- readOGR("dem1_border.shp")
  
  # view polygons
  plot(polyg1)
  plot(polyg2, add = T)
  
  ##### 
  # combine basic line and polygonal border
  p1_1 <- which.min(sapply(X = 1:nrow(polyg1@polygons[[1]]@Polygons[[1]]@coords), FUN = function(x) dist(rbind(matrix(polyg1@polygons[[1]]@Polygons[[1]]@coords[x,],1), matrix(coordinates(ps)[1,1:2],1)))))
  p1_2 <- which.min(sapply(X = 1:nrow(polyg1@polygons[[1]]@Polygons[[1]]@coords), FUN = function(x) dist(rbind(matrix(polyg1@polygons[[1]]@Polygons[[1]]@coords[x,],1), matrix(coordinates(ps)[2,1:2],1)))))
  
  p2_1 <- which.min(sapply(X = 1:nrow(polyg2@polygons[[1]]@Polygons[[1]]@coords), FUN = function(x) dist(rbind(matrix(polyg2@polygons[[1]]@Polygons[[1]]@coords[x,],1), matrix(coordinates(ps)[1,1:2],1)))))
  p2_2 <- which.min(sapply(X = 1:nrow(polyg2@polygons[[1]]@Polygons[[1]]@coords), FUN = function(x) dist(rbind(matrix(polyg2@polygons[[1]]@Polygons[[1]]@coords[x,],1), matrix(coordinates(ps)[2,1:2],1)))))
  
  m1 <- mean(c(dist(rbind(coordinates(ps)[1,1:2], head(con1@lines[[1]]@Lines[[1]]@coords,1))),dist(rbind(coordinates(ps)[2,1:2], tail(con1@lines[[1]]@Lines[[1]]@coords,1)))))
  m2 <- mean(c(dist(rbind(coordinates(ps)[2,1:2], head(con1@lines[[1]]@Lines[[1]]@coords,1))),dist(rbind(coordinates(ps)[1,1:2], tail(con1@lines[[1]]@Lines[[1]]@coords,1)))))
  
  if (m1 < m2) {
    line1 <- Line(coords = rbind(polyg1@polygons[[1]]@Polygons[[1]]@coords[p1_1,], coordinates(ps)[1,1:2], con1@lines[[1]]@Lines[[1]]@coords, coordinates(ps)[2,1:2], polyg1@polygons[[1]]@Polygons[[1]]@coords[p1_2,]))
    
  } else {
    line1 <- Line(coords = rbind(polyg1@polygons[[1]]@Polygons[[1]]@coords[p1_2,], coordinates(ps)[2,1:2], con1@lines[[1]]@Lines[[1]]@coords, coordinates(ps)[1,1:2], polyg1@polygons[[1]]@Polygons[[1]]@coords[p1_1,] ))
  }
  line1 <- Lines(slinelist = line1, ID = "line1")
  line1 <- SpatialLines(LinesList = list(line1))
  crs(line1) <- crs(dem1)
  
  lines(line1, col ='blue')
  line1$ID <- 1:length(line1)
  unlink(list.files(pattern = "line1"))
  # write basic line 1
  writeOGR(line1, "line1.shp", "line1.shp", "ESRI Shapefile")
  
  if (m1 < m2) {
    line2 <- Line(coords = rbind(polyg2@polygons[[1]]@Polygons[[1]]@coords[p2_1,], coordinates(ps)[1,1:2], con1@lines[[1]]@Lines[[1]]@coords, coordinates(ps)[2,1:2], polyg2@polygons[[1]]@Polygons[[1]]@coords[p2_2,]))
  } else {
    line2 <- Line(coords = rbind(polyg2@polygons[[1]]@Polygons[[1]]@coords[p2_2,], coordinates(ps)[2,1:2], con1@lines[[1]]@Lines[[1]]@coords, coordinates(ps)[1,1:2], polyg2@polygons[[1]]@Polygons[[1]]@coords[p2_1,]))
  }
  line2 <- Lines(slinelist = line2, ID = "line2")
  line2 <- SpatialLines(LinesList = list(line2))
  crs(line2) <- crs(dem1)
  lines(line2, col = "pink")
  line2$ID <- 1:length(line2)
  unlink(list.files(pattern = "line2"))
  # write basic line 2
  writeOGR(line2, "line2.shp", "line2.shp", "ESRI Shapefile")
  
  # cut polygonal border 1 by line 1
  rsaga.geoprocessor("shapes_polygons", 8, param = list(
    POLYGONS = "dem2_border.shp",
    LINES = "line1.shp",
    INTERSECT = "intersect1"
  ))
  
  polyg1 <- (readOGR("intersect1.shp"))
  con1 <- polyg1[which.max(area(polyg1)),]
  
  rsaga.geoprocessor("shapes_polygons", 8, param = list(
    POLYGONS = "dem1_border.shp",
    LINES = "line2.shp",
    INTERSECT = "intersect2"
  ))
  
  polyg2 <- (readOGR("intersect2.shp"))
  con2 <- polyg2[which.max(area(polyg2)),]
  # view cutted area
  plot(con2, col = "green")
  plot(con1, add = T, col = "red")
  
  # create small buffers for overlapping
  con1 <- gBuffer(con1, w = w)
  con2 <- gBuffer(con2, w = w)
  
  res_r <- trim(merge(mask(crop(dem2new, con1), con1), mask(crop(dem1, con2), con2)))
  res_r <- resample(res_r, dem_temp)
  crs(res_r) <- crs(dem_temp)
  plot(res_r, main = i)
  writeRaster(res_r, paste0("partial_DEM_",i,"_m.tif"), "GTiff", overwrite=TRUE)
}

sapply(X = c(-25,25), FUN = function(x) tilt(x)) # X meters. Vertical shift of side (shoulder-end) point

time_finish <- timestamp()

# remove temporary data
rm(list = c("dem1", "dem2", "ext1", "dem_temp", "DoD", "DoD_border", "DoD_inside_buffer", "dem1cr", "dem2cr", "shift_size", "expand1",
            "cl1", "err1", "dem2_corr", "bdem1", "bdem2", "bdem3", "ps", "n1", "n2", "zdiff", "coord", "mean_point", "angle1", "shoulder",
            "new_point", "data1", "lm1", "dem2p", "lm1p", "v_error", "dem2new", "dem2new_backup", "own_points", "tilt",
            "i", "data2", "lm2", "lm2p", "dem_diff", "con1", "con2", "dem1_border", "dem2_border", "polyg1", "polyg2", 
            "p1_1", "p1_2", "p2_1", "p2_2", "m1", "m2", "line1", "line2", "res_r"))
unlink(list.files(pattern = "_border|line|intersect"))

#### staking data
(list_r <- list.files(pattern = "partial"))
stack_1 <- do.call(what = stack, args = list(list_r))
stack_1 <- mean(stack_1 )
plot(stack_1)

# write final DEM
writeRaster(stack_1 , "final_DEM .tif", "GTiff", overwrite=TRUE)