#library(bioRad)

#unique_year_week_df <- unique(data.frame(year = closest_times$year, week = closest_times$week))
#unique_year_week_list <- split(unique_year_week_df, seq(nrow(unique_year_week_df)))
#unique_year_week_list[1]
#df_subset = subset(closest_times, year == unique_year_week_list[[1]]$year & week == unique_year_week_list[[1]]$week)
#df_subset$closest_file[1]

#my_pvolfile %>% read_pvolfile() -> my_pvol
#my_vp <- calculate_vp(my_pvolfile, n_layer=50, h_layer=50, sd_vvp_threshold = 1)
#my_ppi_integrated <- integrate_to_ppi(pvol=my_pvol,vp=my_vp,
#                                      xlim = c(-50000, 50000), ylim = c(-50000, 50000),
#                                      res=250)
#plot(my_ppi_integrated, param='VID')

#bm <- "cartolight"
#map(my_ppi_integrated, map=bm, param = "VID", alpha = .8)

closest_times <- read.csv("pvol/KRIW_closest_times.csv")
random_file <- sample(closest_times$closest_file, size = 1, replace = FALSE)
input_pvol_path <- paste0("pvol/", radar, "/", random_file)
apply_mistnet(input_pvol_path) -> my_pvol

get_elevation_angles <- function(pvol) {
  sapply(pvol$scans, function(scan) round(scan$attributes$where$elangle,1))
}

elevation_angles <- get_elevation_angles(my_pvol)


pvol_files = list.files('pvol/biofree/KRIW/raw', full.names = TRUE)
my_pvol = apply_mistnet(pvol_files[1])

my_clean_pvol <- apply_mistnet(pvol_files[1], sep='') %>% clean_pvol()
dbzh = my_clean_pvol$scans[[1]]$params$DBZH
class(dbzh) <- "matrix"
rast(dbzh)




my_scan <- my_pvol$scans[[1]]
my_ppi <- project_as_ppi(my_scan, range_max = 100000)
plot(my_ppi)

occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')

library(terra)

get_elevation_angles <- function(pvol) {
  sapply(pvol$scans, function(scan) scan$attributes$where$elangle)
}

find_proximal_occult_scan_idx <- function(occult, pvol) {
  occult_angles <- get_elevation_angles(occult)
  pvol_angles <- get_elevation_angles(pvol)
  max_occult_angle <- max(occult_angles)
  mapping <- list()
    for (i in seq_along(pvol_angles)) {
    if (pvol_angles[i] <= max_occult_angle) {
      differences = abs(occult_angles - pvol_angles[i])
      closest_idx = which.min(differences)
      mapping[[as.character(i)]] <- closest_idx
    }
  }
  return(mapping)
}

occult_index_mapping <- find_proximal_occult_scan_idx(occult, my_pvol)

for (pvol_idx in seq_along(occult_index_mapping)){

  occult_idx <- occult_index_mapping[[pvol_idx]]

  # Extract information from the pvol scan
  pvol_rscale <- my_pvol$scans[[pvol_idx]]$geo$rscale
  pvol_ascale <- my_pvol$scans[[pvol_idx]]$geo$ascale
  pvol_data <- my_pvol$scans[[pvol_idx]]$params$DBZH
  class(pvol_data) <- "matrix"

  # Extract information from the matched occult scan
  occult_rscale <- occult$scans[[occult_idx]]$geo$rscale
  occult_ascale <- occult$scans[[occult_idx]]$geo$ascale
  occult_data <- occult$scans[[occult_idx]]$params$OCCULT
  class(occult_data) <- "matrix"


  pvol_raster <- rast(pvol_data)
  fact_x <- pvol_rscale / occult_rscale
  fact_y <- pvol_ascale / occult_ascale

  occult_aggregated = aggregate(occult_raster, fact=c(fact_x, fact_y), fun=mean, na.rm=TRUE)
  ext(occult_aggregated) <- ext(pvol_raster)
  occult_resampled <- terra::resample(occult_aggregated, pvol_raster, method="bilinear")

  mask <- occult_resampled > threshold
  pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
  pvol_raster_masked_matrix <- as.matrix(pvol_raster_masked, nrows=pvol_scan_dim[1], ncols=pvol_scan_dim[2])
  class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
  my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix

}


test_pvol_rast = rast(pvol_data)
test_occult_rast = rast(occult_data)



# Create SpatRaster for polar volume scan (pvol)
pvol_extent <- extent(my_pvol$scans[[1]]$geo$lon, my_pvol$scans[[1]]$geo$lon, my_pvol$scans[[1]]$geo$lat, my_pvol$scans[[1]]$geo$lat)
pvol_raster <- rast(nrows=230, ncols=360, extent=pvol_extent, crs="+proj=longlat", resolution=c(my_pvol$scans[[1]]$geo$rscale, my_pvol$scans[[1]]$geo$rscale))

# Create SpatRaster for occultation scan
occult_extent <- extent(occult$scans[[1]]$geo$lon, occult$scans[[1]]$geo$lon, occult$scans[[1]]$geo$lat, occult$scans[[1]]$geo$lat)
occult_raster <- rast(occult$scans[[1]]$params$OCCULT, extent=occult_extent, crs="+proj=longlat", resolution=c(occult$scans[[1]]$geo$rscale, occult$scans[[1]]$geo$rscale))

# Resample the occultation raster to match the polar volume scan resolution
occult_resampled <- resample(occult_raster, pvol_raster, method="bilinear")

# Now you can apply the masking logic based on the threshold as needed
# Assuming you need to mask the polar volume scan where occult exceeds a threshold
blockage_threshold = 0.1
occult_data_resampled <- terra::as.matrix(occult_resampled)
blockage_mask <- occult_data_resampled > blockage_threshold
# Apply the mask to the polar volume scan
# (Further code would be needed here to replace values with NA based on the mask)


write_pvolfile(my_pvol, 'test_pvol.h5')
my_pvol = read_pvolfile('test_pvol.h5')
my_scan <- my_pvol$scans[[1]]
my_ppi <- project_as_ppi(my_scan, range_max = 100000)


range_max = 50000
k = 4 / 3
elev = my_scan$geo$elangle


proj4string <- sp::CRS(paste("+proj=aeqd +lat_0=", my_scan$geo$lat,
                             " +lon_0=", my_scan$geo$lon,
                             " +units=m",
                             sep = ""
))


bboxlatlon <- proj_to_wgs(
  c(-range_max, range_max),
  c(-range_max, range_max),
  proj4string
)@bbox


index <- polar_to_index(
  cartesian_to_polar(sp::coordinates(gridTopo), elev, k = k, lat = attributes(param)$geo$lat, re = re, rp = rp),
  rangebin = attributes(param)$geo$rscale,
  azimbin = attributes(param)$geo$ascale,
  azimstart = ifelse(is.null(attributes(param)$geo$astart), 0, attributes(param)$geo$astart),
  rangestart = ifelse(is.null(attributes(param)$geo$rstart), 0, attributes(param)$geo$rstart)
)


occult= read_pvolfile('occult/h5-150km/KRIW.h5', param='all')
my_occult_scan <- occult$scans[[1]]


occult= read_pvolfile('occult/KRIW.h5', param='all')
my_scan <- occult$scans[[1]]
plot(my_scan, zlim=c(0,1))
my_ppi <- project_as_ppi(my_scan, range_max=150000)
plot(my_ppi,  zlim=c(0,1))




occult_interp_correct <- function(occult_scan, target_nbins) {
  # Original range bins and target range bins
  original_bins <- seq_len(nrow(occult_scan))
  target_bins <- seq_len(target_nbins)

  # Interpolating each ray
  occult_interp <- apply(occult_scan, 2, function(ray) {
    approx(x = original_bins, y = ray, xout = target_bins)$y
  })

  return(occult_interp)
}

occult_scan_data <- my_occult_scan$params$OCCULT  # Extracting OCCULT data
my_occult_scan_interp <- occult_interp_correct(occult_scan_data, 230)


grid_size = 500

class(my_scan) <-  "matrix"


write_pvolfile(occult, 'test_occult.h5', overwrite = TRUE)
test_pvol = read_pvolfile('test_occult.h5', param='all')
my_scan <- test_pvol$scans[[1]]
plot(my_scan)

#fix the plot method for polar volumes plot.scan , plot.ppi
#remake occultation maps



occult= read_pvolfile('KRIW.h5', param='all')
my_scan <- occult$scans[[1]]


my_ppi <- project_as_ppi(occult$scans[[1]], range_max=150000)


x = my_scan
param <- names(x$params)[1]
data <- do.call(function(y) x$params[[y]], list(param))


if(is.null(x$geo$rstart) || length(x$geo$rstart) == 0) {
  rstart <- 0
} else {
  rstart <- x$geo$rstart
}

raster_obj <- raster::raster(t(data), ymn = astart, ymx = astart + 360, xmn = rstart, xmx = rstart + rscale * dimraster[1])


my_ppi <- project_as_ppi(param, range_max = 100000)
plot(my_ppi)

plot(project_as_ppi(occmap$scans[[1]], range_max=150000),zlim=c(0,1))


#Six maps from 2019

indices <- c(length(pvol_files) - 120, length(pvol_files) - 100, length(pvol_files) - 80,
             length(pvol_files) - 60, length(pvol_files) - 40, length(pvol_files) - 20)


for (i in indices) {
  file <- pvol_files[i]
  my_pvol = read_pvolfile(paste('pvol/KRIW/', file, sep=''))
  my_scan <- my_pvol$scans[[1]]
  my_ppi <- project_as_ppi(my_scan,range_max=150000)

  # Define the path for the output plot
  plot_path <- paste("plots/pvol_2007/", file, ".png", sep="")
  # Open a PNG device
  png(filename=plot_path, width=800, height=600)
  print(plot(my_ppi))
  dev.off() # Close the PNG device
  Sys.sleep(1)
}

#How to upsample my occult map to filter my polar volume?

NEXRAD_antenna_height <- readRDS("~/clo/bigbird/wgfd/NEXRAD_antenna_height.rds")
antenna=NEXRAD_antenna_height
radpos = antenna



library(sp)
library(sf)

my_path <- sf::read_sf("umbral_path.geojson")
my_path_xy <- st_zm(my_path, what = "ZM")
my_path_sp <- as(my_path_xy, "Spatial")
lines <- as(my_path_sp, "SpatialLines")

# Extract the Northern and Southern Limits
northern_limit <- lines@lines[[1]]@Lines[[1]]@coords
southern_limit <- lines@lines[[3]]@Lines[[1]]@coords

adjusted_southern_limit <- adjusted_southern_limit[nrow(adjusted_southern_limit):1, ]
coords_polygon <- rbind(northern_limit, adjusted_southern_limit, northern_limit[1, ])

poly <- SpatialPolygons(list(Polygons(list(Polygon(coords_polygon)), "poly1")))
data <- data.frame(id = 1)  # Dummy data frame for SpatialPolygonsDataFrame
row.names(data) <- "poly1"
sp_poly_df <- SpatialPolygonsDataFrame(poly, data)

plot(sp_poly_df, col = 'lightblue', main = "Polygon from Northern and Southern Limits")
saveRDS(sp_poly_df, "umbral_path_poly.rds")


sf::st_read()


sp_poly_df <- readRDS("umbral_path_poly.rds")
#umbral_point <- readRDS("umbral_point.rds")
umbra_trajectory= sf::st_read('2024_eclipse_shapefiles/umbra_lo.shp')

current_utc <- Sys.time()

umbra_trajectory <- umbra_trajectory %>%
  mutate(
    DateTime = as.POSIXct(UTCTime, format = "%H:%M:%S", tz = "UTC"),
    DateTime = DateTime + as.difftime(as.Date(current_utc) - as.Date(DateTime))
  )

closest_feature <- umbra_trajectory %>%
  filter(abs(difftime(DateTime, current_utc, units = "secs")) == min(abs(difftime(DateTime, current_utc, units = "secs")))) %>%
  slice(1)  

closest_feature <- st_transform(closest_feature, CRS("+init=epsg:4326"))
closest_feature_sp <- as(closest_feature, "Spatial")
centroid <- gCentroid(closest_feature_sp)

target_crs <- proj4string(sp_poly_df)
umbral_point <- spTransform(centroid, CRS(target_crs))



library(sf)
library(sp)
library(rgeos)
library(dplyr)

sp_poly_df <- readRDS("umbral_path_poly.rds")
umbra_trajectory <- sf::st_read('2024_eclipse_shapefiles/umbra_lo.shp')

umbra_trajectory <- sf::st_read('2024_eclipse_shapefiles/umbra_lo.shp')
compute_umbral_point <- function(umbra_trajectory, sp_poly_df) {
  tryCatch({
    
    current_utc <- Sys.time() + hours(3)
    print(current_utc)
    
    umbra_trajectory <- umbra_trajectory %>%
      mutate(
        DateTime = as.POSIXct(paste("2024-04-08", UTCTime), format = "%Y-%m-%d %H:%M:%S")
      ) %>%
      filter(abs(difftime(DateTime, current_utc, units = "secs")) == min(abs(difftime(DateTime, current_utc, units = "secs")))) %>%
      slice(1)
    
    closest_feature_transformed <- st_transform(umbra_trajectory, st_crs(sp_poly_df))
    centroid <- st_centroid(closest_feature_transformed)
    
    return(centroid)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    return(NULL)
  })
}

umbral_point <- compute_umbral_point(umbra_trajectory, sp_poly_df)
plot(sp_poly_df)


sp_poly_df <- readRDS("umbral_path_poly.rds")
umbra_trajectory <- sf::st_read('2024_eclipse_shapefiles/umbra_lo.shp')

compute_umbral_point <- function(umbra_trajectory, sp_poly_df) {
  tryCatch({
    
    current_utc <- Sys.time() 
    print(current_utc)
    
    umbra_trajectory <- umbra_trajectory %>%
      mutate(
        DateTime = as.POSIXct(paste("2024-04-08", UTCTime), format = "%Y-%m-%d %H:%M:%S")
      ) %>%
      filter(abs(difftime(DateTime, current_utc, units = "secs")) == min(abs(difftime(DateTime, current_utc, units = "secs")))) %>%
      slice(1)
    
    closest_feature_transformed <- st_transform(umbra_trajectory, st_crs(sp_poly_df))
    centroid <- st_centroid(closest_feature_transformed)
    
    return(centroid)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    return(NULL)
  })
}

umbral_point <- compute_umbral_point(umbra_trajectory, sp_poly_df)
coords <- st_coordinates(umbral_point)
spatial_points <- SpatialPoints(coords, proj4string = CRS(st_crs(umbral_point)$proj4string))
umbral_point_sp <- SpatialPointsDataFrame(spatial_points, data=data.frame(row.names=row.names(spatial_points)))
print(umbral_point_sp)

plot(sp_poly_df)
points(umbral_point_sp)



umbra_trajectory %>%
  mutate(
    DateTime = as.POSIXct(paste("2024-04-08", UTCTime), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  )


library(dplyr)
library(sp)
library(rgdal)

# Ensure 'umbra_trajectory' and 'sp_poly_df' are defined somewhere before this.

compute_umbral_point <- function(umbra_trajectory, sp_poly_df) {
  tryCatch({

    current_utc <- Sys.time()
    
    umbra_trajectory <- umbra_trajectory %>%
      mutate(
        DateTime = as.POSIXct(UTCTime, format = "%H:%M:%S", tz = "UTC"),
        DateTime = DateTime + as.difftime(as.Date(current_utc) - as.Date(DateTime), units = "days")
      ) %>%
      filter(abs(difftime(DateTime, current_utc, units = "secs")) == min(abs(difftime(DateTime, current_utc, units = "secs")))) %>%
      slice(1)
    

    closest_feature_transformed <- st_transform(umbra_trajectory, st_crs(sp_poly_df))
    centroid <- st_centroid(closest_feature_transformed)
    
    return(centroid)
  }, error = function(e) {
    cat("An error occurred: ", e$message, "\n")
    return(NULL)
  })
}

umbral_point <- compute_umbral_point(umbra_trajectory, sp_poly_df)







umbral_point
plot(sp_poly_df)
points(umbral_point)




plot(lower48, bg="black", border="white")
plot(sp_poly_df, add = TRUE, col = adjustcolor("#005e6a", alpha.f=0.5),border=NA)
points(umbral_point, pch = 20,  cex = 1.5)

#plot(umbral_point, add=TRUE, pch = 20, col = "#001315", cex = 1.5)



coords <- lines@lines[[2]]@Lines[[1]]@coords
random_point = coords[600,]
sp_poly_df <- spTransform(sp_poly_df, CRS(proj4string(lower48)))

points_sp <- SpatialPoints(matrix(random_point, nrow = 1, byrow = TRUE),
                           proj4string = CRS("+init=epsg:4326"))

points_sp <- spTransform(points_sp, CRS(proj4string(sp_poly_df)))


png(filename = "your_plot_name.png", width = 800, height = 600, res = 100)
plot(lower48,lwd=2)
plot(sp_poly_df, add = TRUE, col = adjustcolor("#005e6a", alpha.f=0.5),border=NA)
points(umbral_point, pch = 20, col = adjustcolor("#001315", alpha.f=0.7), cex = 4.5)
dev.off()





library(sf)
library(ggplot2)
library(rgeos)

# Define a sequence for x values
x_values <- seq(0, 10, by = 0.1)

# Generate two sinusoidal lines with a phase shift
y_values_line1 <- sin(x_values)
y_values_line2 <- sin(x_values) + 2 # Shifted upwards to make them parallel

# Combine x and y to create coordinates for both lines
coords_line1 <- cbind(x_values, y_values_line1)
coords_line2 <- cbind(x_values, y_values_line2)

# Create LINESTRING objects for both lines
line1 <- st_sfc(st_linestring(coords_line1), crs = 4326)
line2 <- st_sfc(st_linestring(coords_line2), crs = 4326)

# To create a polygon, combine coordinates of line1 with those of line2 in reverse
coords_polygon <- rbind(coords_line1, rev(coords_line2), coords_line1[1, ])

# Create the polygon from combined coordinates
polygon <- st_sfc(st_polygon(list(coords_polygon)), crs = 4326)

# Plot
ggplot() +
  geom_sf(data = st_as_sf(line1), color = "blue") +
  geom_sf(data = st_as_sf(line2), color = "red") +
  geom_sf(data = st_as_sf(polygon), fill = "lightblue", alpha = 0.5) +
  theme_minimal()


# Plot polygons with fill and no borders
plot(lower48, border=NA)
plot(lower48, col=NA, border=adjustcolor("black", alpha.f=0.5), add=TRUE)


lower48_sf <- st_as_sf(lower48)
intersecting_parts <- list()

for (i in seq_len(nrow(lower48_sf))) {
  current_polygon <- lower48_sf[i, ]
  
  for (j in seq_len(nrow(lower48_sf))) {
    if (i != j) {
      intersecting_part <- st_intersection(current_polygon, lower48_sf[j, ])
      if (nrow(intersecting_part) > 0) {
        intersecting_parts[[length(intersecting_parts) + 1]] <- st_geometry(intersecting_part)
      }
    }
  }
}

geom_collection_inner = st_union(combined_geometries)
geometry_list <- list()

for (i in 1:nrow(lower48_sf)) {
  current_geometry <- st_geometry(lower48_sf[i, ])
  geometry_list[[i]] <- current_geometry
}
all_geometries_combined <- do.call(c, geometry_list)
geom_collection_lower48  <- st_union(all_geometries_combined)

outline_geometry <- st_difference(geom_collection_lower48, geom_collection_inner)


sp_poly_df <- readRDS("umbral_path_poly.rds")
sp_poly_sf <- st_as_sf(sp_poly_df)
radar_obscuration = read.csv("radar_obscuration.csv")
coordinates(radar_obscuration) <- ~lon+lat  # Defining the longitude and latitude columns
proj4string(radar_obscuration) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
radar_obscuration_transformed <- spTransform(radar_obscuration, CRS(proj4string(lower48wgs)))
radar_obscuration_sf <- st_as_sf(radar_obscuration_transformed)

bbox_outline <- st_bbox(outline_geometry)
clipped_sp_poly_sf <- st_crop(sp_poly_sf, bbox_outline)
clipped_sp_poly_sf = st_transform(clipped_sp_poly_sf, st_crs(radar_obscuration_sf))

totality_radar <- st_intersection(radar_obscuration_sf, clipped_sp_poly_sf)



ggplot(outline_geometry) +
    geom_sf(fill = NA, color = "black")  +
    theme_void() +
    geom_sf(data = radar_obscuration_sf, inherit.aes = FALSE, color = "red", size = 2)  # Add points

caption_text = 
"  The path of eclipse totality through the network of NEXRAD weather surveillance radars in the continental USA, 8 April 2024. All 143 sites are colored by the maximum amount of 
  obscuration. The path of totality, where obscuration is 100%, is shown in grey, and the 13 sites located within the path of totality are outlined in black. Source: NASA"

upath17 <- geojson_sf('upath17_new.geojson')
upath17_transformed = st_transform(upath17, st_crs(lower48wgs))



my_plot <- ggplot() +
  geom_sf(data = outline_geometry, fill = NA, color = "black") +  # Plot the outline geometry
  geom_sf(data = clipped_sp_poly_sf, fill = adjustcolor("#D3D3D3", alpha = 0.8), color = NA) +  # Polygon layer with transparency
  geom_sf(data = upath17_transformed, fill = adjustcolor("#D3D3D3", alpha = 0.8), color = NA) +  # Polygon layer with transparency
  #geom_sf(data = totality_radar, color = "black", size = 7) +  # Plot points colored by max_percent
  geom_sf(data = radar_obscuration_sf[radar_obscuration_sf$radar == 'KPAH',], aes(color = max_percent), size = 2) + 
  geom_sf_label(data = radar_obscuration_sf[radar_obscuration_sf$radar == 'KPAH',], aes(label = radar),nudge_y = 500000) + # Plot points colored by max_percent
  #scale_color_gradient(low = "#FFFF54", high = "#4D4D4D",  # Color gradient from yellow to dark grey
  #                     name = "obscuration %",  # Label for the color scale in the legend
  #                     limits = c(0, 100),  # Set limits for the color scale
  #                     guide = "colourbar") +
  theme_void() + 
  #theme(plot.caption = element_text(hjust =0, family="sans")) +
  #labs(caption = caption_text) + 
  theme(legend.position = "none") +
  theme(legend.title = element_text(margin = margin(b = 5)))  

ggsave("totality_overlap.png", plot = my_plot, width = 10, height = 8, dpi = 600, bg='#ffffff')





img1 <- readPNG("high_res_plot_aug_21.png")
img2 <- readPNG("high_res_plot_apr_8.png")

grob1 <- rasterGrob(img1, interpolate=TRUE)
grob2 <- rasterGrob(img2, interpolate=TRUE)
grid.arrange(grob1, grob2, ncol=2)

png("combined_high_res.png", width = 5, height = 3, units = "in", res = 600, bg = "#ffffff")
grid.arrange(grob1, grob2, ncol=2)
dev.off()


library(readr)
library(suntools)

x2 <- read_delim("~/Downloads/Example_data_sunriset.csv", delim=";")

x2$sunrise_test <- sunriset(
  as.matrix(data.frame(x2$Longitude, x2$Latitude)),
  as.POSIXct(x2$local_time, format = "%Y-%m-%d %H:%M:%S", tz="UTC"),
  direction = "sunrise",
  POSIXct.out = TRUE
)
x2$sunset_test <- sunriset(
  as.matrix(data.frame(x2$Longitude, x2$Latitude)),
  as.POSIXct(x2$local_time, format = "%Y-%m-%d %H:%M:%S", tz="UTC"),
  direction = "sunset",
  POSIXct.out = TRUE
)

df = data.frame(lat=x2$Latitude, lon=x2$Longitude, time=x2$local_time, sunrise=x2$sunrise_test$time, sunset=x2$sunset_test$time)




seq(as.POSIXct("2023-08-01 00:00:00",
               format = "%Y-%m-%d 00:00:00",
               tz = "Australia/Perth",
               force_tz = TRUE),
    as.POSIXct("2023-08-05 00:00:00",
               format = "%Y-%m-%d 00:00:00",
               tz = "Australia/Perth",
               force_tz = TRUE),
    by = '1 day')






library(tibble)
library(tibbletime)
library(lubridate)

#This doesn't work:
dat <- tibble(Date = seq(as.POSIXct("2023-08-01 00:00:00",
                                    format = "%Y-%m-%d 00:00:00",
                                    tz = "Australia/Perth",
                                    force_tz = TRUE),
                         as.POSIXct("2023-08-05 00:00:00",
                                    format = "%Y-%m-%d 00:00:00",
                                    tz = "Australia/Perth",
                                    force_tz = TRUE),
                         by = '1 day'))
dat
x <- dat[1,1]
x
tz(x)

tz(dat[1,1]$Date)
[1] "Australia/Perth"

Date = dat[1,1]$Date



#Roebuck, AU

dat <- dat %>% mutate(Sunrise = bioRad::sunrise(Date, lon = 122.48, lat = -17.846, tz=tz(Date), force_tz = TRUE),
                      Sunset = bioRad::sunset(Date, lon = 122.48, lat = -17.846, tz=tz(Date), force_tz = TRUE))



dat <- dat %>% mutate(diffSR = Sunrise +as.POSIXct("00-00-01 00:00:00"),
                      diffSS = Sunset + as.POSIXct("00-00-01 00:00:00"))






#bioRad - default tz
dat <- dat %>% mutate(Sunrise = bioRad::sunrise(Date, lon = 122.48, lat = -17.846),
                      Sunset = bioRad::sunset(Date, lon = 122.48, lat = -17.846))
dat

dat <- dat %>% mutate(diffSR = Sunrise +as.POSIXct("00-00-01 00:00:00"),
                      diffSS = Sunset + as.POSIXct("00-00-01 00:00:00"))




#This does
dat <- data.frame(Date = seq(as.POSIXct("2023-08-01 00:00:00",
                                        format = "%Y-%m-%d 00:00:00",
                                        tz = "Australia/Perth",
                                        force_tz = TRUE),
                             as.POSIXct("2023-08-05 00:00:00",
                                        format = "%Y-%m-%d 00:00:00",
                                        tz = "Australia/Perth",
                                        force_tz = TRUE),
                             by = '1 day'))
dat
x <- dat[1,1]
x
tz(x)

dat <- dat %>% mutate(Sunrise = bioRad::sunrise(Date, lon = 122.48, lat = -17.846,
                                                tz = "Australia/Perth", force_tz = TRUE),
                      Sunset = bioRad::sunset(Date, lon = 122.48, lat = -17.846,
                                              tz = "Australia/Perth", force_tz = TRUE))




