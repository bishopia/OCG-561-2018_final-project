###---Polar map for poster---###


######################################
# THIS SCRIPT SECTION WAS LIFTED FROM:
# https://bit.ly/2Kp2pRE
#######################################

library(sp)
library(maps)
library(rgeos)


# function to slice and dice a map and convert it to an sp() object
maps2sp = function(xlim, ylim, l.out = 100, clip = TRUE) {
  stopifnot(require(maps))
  m = map(xlim = xlim, ylim = ylim, plot = FALSE, fill = TRUE)
  p = rbind(cbind(xlim[1], seq(ylim[1],ylim[2],length.out = l.out)),
            cbind(seq(xlim[1],xlim[2],length.out = l.out),ylim[2]),
            cbind(xlim[2],seq(ylim[2],ylim[1],length.out = l.out)),
            cbind(seq(xlim[2],xlim[1],length.out = l.out),ylim[1]))
  LL = CRS("+init=epsg:4326")
  IDs = sapply(strsplit(m$names, ":"), function(x) x[1])
  stopifnot(require(maptools))
  m = map2SpatialPolygons(m, IDs=IDs, proj4string = LL)
  bb = SpatialPolygons(list(Polygons(list(Polygon(list(p))),"bb")), proj4string = LL)
  
  if (!clip)
    m
  else {
    stopifnot(require(rgeos))
    gIntersection(m, bb)
  }
}

# set colours for map grid
grid.col.light = rgb(0.5,0.5,0.5,0.8)
grid.col.dark = rgb(0.5,0.5,0.5)

# coordinate systems
polar = CRS("+init=epsg:3995")
longlat = CRS("+init=epsg:4326")

# download the blue marble data if it doesn't
# exist
if (!file.exists("blue_marble.tif")) {
  download.file("http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=526312&cs=rgb&format=TIFF&width=5400&height=2700","blue_marble.tif")
}

# read in the raster map and
# set the extent, crop to extent and reproject to polar
r = raster::brick("blue_marble.tif")
e = raster::extent(c(-180,180,55,90))
r_crop = raster::crop(r,e)

# traps NA values and sets them to 1
r_crop[is.na(r_crop)] = 1 
r_polar = raster::projectRaster(r_crop, crs = polar, method = "bilinear")

# some values are not valid after transformation 
# (rgb range = 1 - 255) set these back to 1
# as they seem to be the black areas
r_polar[r_polar < 1 ] = 1

# define the graticule / grid lines by first specifying
# the larger bounding box in which to place them, and
# feeding this into the sp() gridlines function
# finally the grid lines are transformed to
# the EPSG 3995 projection
pts=SpatialPoints(rbind(c(-180,55),c(0,55),c(180,85),c(180,85)), CRS("+init=epsg:4326"))
gl = gridlines(pts, easts = seq(-180,180,30), norths = seq(50,85,10), ndiscr = 100)
gl.polar = spTransform(gl, polar)

# I also create a single line which I use to mark the
# edge of the image (which is rather unclean due to pixelation)
# this line sits at 55 degrees North similar to where I trimmed
# the image
pts=SpatialPoints(rbind(c(-180,55),c(0,55),c(180,80),c(180,80)), CRS("+init=epsg:4326"))
my_line = SpatialLines(list(Lines(Line(cbind(seq(-180,180,0.5),rep(55,721))), ID="outer")), CRS("+init=epsg:4326"))

# crop a map object (make the x component a bit larger not to exclude)
# some of the eastern islands (the centroid defines the bounding box)
# and will artificially cut of these islands
m = maps2sp(c(-180,200),c(55,90),clip = TRUE)

#----- below this point is the plotting routine
# set margins to let the figure "breath" and accommodate labels
par(mar=rep(1,4))

# plot the grid, to initiate the area
# plotRGB() overrides margin settings in default plotting mode
plot(spTransform(gl, polar), lwd=2, lty=2,col="white")

# plot the blue marble raster data
raster::plotRGB(blue_marble, add = TRUE)

# plot grid lines / graticule
lines(spTransform(gl, polar), add = TRUE, lwd=2, lty=2,col=grid.col.light)

# plot outer margin of the greater circle
lines(spTransform(ll, polar), lwd = 3, lty = 1, col=grid.col.dark)

# plot continent outlines, for clarity
plot(spTransform(m, polar), lwd = 1, lty = 1, col = "transparent", border=grid.col.dark, add = TRUE)

# plot longitude labels
l = labels(gl.polar, longlat, side = 1)
l$pos = NULL
text(l, cex = 1, adj = c( 0.5, 2 ),  col = "black")

# plot latitude labels
l = labels(gl.polar, longlat, side = 2)
l$srt = 0
l$pos = NULL
text(l, cex = 1, adj = c(1.2, -1), col = "white")

# After all this you can plot your own site locations etc
# but don't forget to tranform the data from lat / long
# into the arctic polar stereographic projection using
# spTransform()

#-------- antarctic map

library(sp)
library(maps)
library(rgeos)

# function to slice and dice a map and convert it to an sp() object
maps2sp = function(xlim, ylim, l.out = 100, clip = TRUE) {
  stopifnot(require(maps))
  m = map(xlim = xlim, ylim = ylim, plot = FALSE, fill = TRUE)
  p = rbind(cbind(xlim[1], seq(ylim[1],ylim[2],length.out = l.out)),
            cbind(seq(xlim[1],xlim[2],length.out = l.out),ylim[2]),
            cbind(xlim[2],seq(ylim[2],ylim[1],length.out = l.out)),
            cbind(seq(xlim[2],xlim[1],length.out = l.out),ylim[1]))
  LL = CRS("+init=epsg:4326")
  IDs = sapply(strsplit(m$names, ":"), function(x) x[1])
  stopifnot(require(maptools))
  m = map2SpatialPolygons(m, IDs=IDs, proj4string = LL)
  bb = SpatialPolygons(list(Polygons(list(Polygon(list(p))),"bb")), proj4string = LL)
  
  if (!clip)
    m
  else {
    stopifnot(require(rgeos))
    gIntersection(m, bb)
  }
}

# set colours for map grid
grid.col.light = rgb(0.5,0.5,0.5,0.8)
grid.col.dark = rgb(0.5,0.5,0.5)

# coordinate systems
polar = CRS("+init=epsg:3031")
longlat = CRS("+init=epsg:4326")

# download the blue marble data if it doesn't
# exist
if (!file.exists("~/blue_marble.tif")) {
  download.file("http://neo.sci.gsfc.nasa.gov/servlet/RenderData?si=526312&cs=rgb&format=TIFF&width=5400&height=2700","~/blue_marble.tif")
}

# read in the raster map and
# set the extent, crop to extent and reproject to polar
r = raster::brick("~/blue_marble.tif")
e = raster::extent(c(-180,180,-90,-55))
r_crop = raster::crop(r,e)

# traps NA values and sets them to 1
r_crop[is.na(r_crop)] = 1
r_polar = raster::projectRaster(r_crop, crs = polar, method = "bilinear")

# some values are not valid after transformation
# (rgb range = 1 - 255) set these back to 1
# as they seem to be the black areas
r_polar[r_polar < 1 ] = 1

# define the graticule / grid lines by first specifying
# the larger bounding box in which to place them, and
# feeding this into the sp() gridlines function
# finally the grid lines are transformed to
# the EPSG 3995 projection
pts = SpatialPoints(rbind(c(-180,-55),c(0,-55),c(180,-85),c(180,-85)), CRS("+init=epsg:4326"))
gl = gridlines(pts, easts = seq(-180,180,30), norths = seq(-85,-55,10), ndiscr = 100)
gl.polar = spTransform(gl, polar)

# I also create a single line which I use to mark the
# edge of the image (which is rather unclean due to pixelation)
# this line sits at 55 degrees North similar to where I trimmed
# the image
pts = SpatialPoints(rbind(c(-180,-55),c(0,-55),c(180,-80),c(180,-80)), CRS("+init=epsg:4326"))
ll = SpatialLines(list(Lines(Line(cbind(seq(-180,180,0.5),rep(-55,721))), ID="outer")), CRS("+init=epsg:4326"))

# crop a map object (make the x component a bit larger not to exclude)
# some of the eastern islands (the centroid defines the bounding box)
# and will artificially cut of these islands
m = maps2sp(c(-180,200),c(-90,-65),clip = FALSE)

#----- below this point is the plotting routine
# set margins to let the figure "breath" and accommodate labels
par(mar=rep(1,4))

# plot the grid, to initiate the area
# plotRGB() overrides margin settings in default plotting mode
plot(spTransform(gl, polar), lwd=2, lty=2, col="white")

# plot the blue marble raster data
raster::plotRGB(r_polar, add = TRUE)

# plot grid lines / graticule
lines(gl.polar, add = TRUE, lwd=2, lty=2,col=grid.col.light)

# plot outer margin of the greater circle
lines(spTransform(ll, polar), lwd = 3, lty = 1, col=grid.col.dark)

# plot continent outlines, for clarity
plot(spTransform(m, polar), lwd = 1, lty = 1, col = "transparent", border=grid.col.dark, add = TRUE)

# plot longitude labels
l = labels(gl.polar, longlat, side = 3)
l$pos = NULL
text(l, cex = 1, adj = c( 0.5, -1),  col = "black")

# plot latitude labels
l = labels(gl.polar, longlat, side = 2)
l$srt = 0
l$pos = NULL
text(l, cex = 1, adj = c(1.2, -1), col = "white")

# After all this you can plot your own site locations etc
# but don't forget to tranform the data from lat / long
# into the arctic polar stereographic projection using
# spTransform()


############################################
# THIS SCRIPT SECTION IS MY OWN CONTRIBUTION
############################################


#convert station lat long to polar project, add to map
stations <- read.csv("site_project_metadata.csv",header = TRUE)

coordinates(stations) <- c("long","lat") 
proj4string(stations) <- sp::CRS("+init=epsg:4326") 
## now we can transform to the projection in use 
# gps.ant <- sp::spTransform(stations, sp::CRS(polar)) 
gps.ant <- sp::spTransform(stations, polar) 

# try to reorder levels to make things easier
# stations$region <- factor(stations$region, levels = c("Drake Passage","Palmer LTER","Amundsen Sea","Ross Sea", "Ross Shelf"))

#change color schematic
colors <- c('#ffffb3', '#bebada','#fb8072','#80b1d3','#8dd3c7')

#plot points over map
plot(gps.ant, add=TRUE, pch=21, col="red", size=4, bg=colors[stations$region])
