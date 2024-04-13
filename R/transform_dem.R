library(stars)
library(terra)
library(sf)
library(raster)
library(bioRad)

setwd("~/clo/bigbird/wgfd")
load("/Users/at744/clo/bitbucket/birdcast/is-birdcast-observed-images/basemap.RData")
load("/Users/at744/clo/bitbucket/birdcast/is-birdcast-observed-images/vgModel.RData")
lower48_sf <- st_as_sf(lower48)

lat = radarInfo[radarInfo$radar == 'KRIW',]$lat
lon = radarInfo[radarInfo$radar == 'KRIW',]$lon
proj4string_az =CRS(paste("+proj=aeqd +lat_0=",lat," +lon_0=",lon," +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep=""))
proj4string_az

dem_raster = raster::brick("gis/elev/elev48i0100a.tif")

# Create a bounding box around the radar location in WGS84
bbox_latlon <- st_bbox(c(xmin = lon - 1.5, xmax = lon + 1.5, ymin = lat - 1.5, ymax = lat + 1.5), 
                       crs = st_crs(4326))

# Transform the bounding box to match the DEM's CRS
bbox_sf <- st_as_sfc(bbox_latlon)
bbox_transformed <- st_transform(bbox_sf, crs = crs(dem_raster))


# Crop the DEM raster using the transformed bounding box
radar_dem <- raster::crop(dem_raster, as(bbox_transformed, "Spatial"))

writeRaster(radar_dem, filename = "radar_dem.tif", format = "GTiff", overwrite = TRUE)


prj_radar_dem <- raster::projectRaster(from=radar_dem, crs=proj4string_az, res=250)



