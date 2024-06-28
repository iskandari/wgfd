library(terra)
library(raster)
library(bioRad)
library(raster)
library(logger)
library(lubridate)
library(ggplot2)
library(dplyr)

setwd("~/clo/bigbird/wgfd")

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


amplify_dbzh <- function(dbzh, occult, threshold = 0.5, wavelength = 10.71) {
  if (!is.na(dbzh) && !is.nan(dbzh) && !is.na(occult) && !is.nan(occult) && occult < threshold) {
    return(eta_to_dbz(dbz_to_eta(dbzh, wavelength) +
                        dbz_to_eta(dbzh, wavelength) * occult, wavelength))
  } else {
    return(dbzh)
  }
}


clean_weather <- function(my_pvol) {

  #filter out non-mistnet scans
  my_pvol$scans <- Filter(function(scan) "CELL" %in% names(scan$params), my_pvol$scans)
  for (i in seq_along(my_pvol$scans)) {

    # for each scan, extract WEATHER and DBZH parameters
    cell <- my_pvol$scans[[i]]$params$CELL
    dbzh <- my_pvol$scans[[i]]$params$DBZH
    rhohv <- my_pvol$scans[[i]]$params$RHOHV

    param_attrs <- attributes(dbzh)

    class(cell) <- "matrix"
    class(dbzh) <- "matrix"
    class(rhohv) <- "matrix"

    # Set DBZH values to NA where CELL > 1
    dbzh[cell > 1] <- NaN
    dbzh[rhohv >= 0.95] <- NaN

    # Restore original class of dbzh and update with the masked scan
    class(dbzh) <- c("param", "matrix", "array")
    my_pvol$scans[[i]]$params$DBZH <- dbzh
    attributes(my_pvol$scans[[i]]$params$DBZH) <-param_attrs
  }
  return(my_pvol)
}


clean_pvol <- function(file, occult, method, threshold, amplify=TRUE) {

  if (!method %in% c('subtraction', 'occultation')) {
    stop("Invalid method. Please choose either 'subtract' or 'occultation'.")
  }

  apply_mistnet(file) -> my_pvol
  my_pvol = clean_weather(my_pvol)

  # proceed based on the method
  if (method == 'subtraction') {

    log_info("subtracting latent reflectivity")

    for (lyr_idx in 1:terra::nlyr(subtract_mask)){

      pvol_data <- my_pvol$scans[[lyr_idx]]$params$DBZH
      param_attrs <- attributes(pvol_data)
      class(pvol_data) <- "matrix"
      pvol_raster <- terra::rast(pvol_data)

      lyr <- subtract_mask[[lyr_idx]]
      values(lyr)[is.nan(values(lyr))] <- NA

      result <- ifel(!is.na(pvol_raster) & !is.na(lyr),
                     pvol_raster - lyr,  #subtract only if both cells are not NA
                     pvol_raster)

      values <- values(result)
      result_mat <- matrix(values, nrow=nrow(lyr), byrow=TRUE)
      class(result_mat) <- c("param", "matrix", "array")
      my_pvol$scans[[lyr_idx]]$params$DBZH <- result_mat
      attributes(my_pvol$scans[[lyr_idx]]$params$DBZH) <- param_attrs
    }

  } else if (method == 'occultation') {


    #match pvol scan to occult by elevation angle
    occult_index_mapping <- find_proximal_occult_scan_idx(occult, my_pvol)

    for (pvol_idx in seq_along(occult_index_mapping)){

      occult_idx <- occult_index_mapping[[pvol_idx]]

      # compare dimensions of xtract the pvol and occult data
      param_attrs <- attributes(my_pvol$scans[[pvol_idx]]$params$DBZH)
      pvol_scan_dim <- dim(my_pvol$scans[[pvol_idx]]$params$DBZH)
      occult_scan_dim <- dim(occult$scans[[occult_idx]]$params$OCCULT)

      log_info('pvol scan dimensions {pvol_scan_dim}')

      pvol_data <- my_pvol$scans[[pvol_idx]]$params$DBZH
      class(pvol_data) <- "matrix"
      pvol_raster <- rast(pvol_data)

      occult_data <- occult$scans[[occult_idx]]$params$OCCULT
      class(occult_data) <- "matrix"
      occult_raster <- rast(occult_data)

      if (!all(pvol_scan_dim == occult_scan_dim)) {

        log_info('dimensions are not equal, resampling...')

        # extract information from the pvol scan
        pvol_rscale <- my_pvol$scans[[pvol_idx]]$geo$rscale
        pvol_ascale <- my_pvol$scans[[pvol_idx]]$geo$ascale

        # extract information from the matched occult scan
        occult_rscale <- occult$scans[[occult_idx]]$geo$rscale
        occult_ascale <- occult$scans[[occult_idx]]$geo$ascale

        fact_x <- pvol_rscale / occult_rscale
        fact_y <- pvol_ascale / occult_ascale

        occult_aggregated = aggregate(occult_raster, fact=c(fact_x, fact_y), fun=mean, na.rm=TRUE)
        ext(occult_aggregated) <- ext(pvol_raster)
        occult_resampled <- terra::resample(occult_aggregated, pvol_raster, method="bilinear")

        if (amplify==TRUE) {
          dbzh_values <- values(pvol_raster)
          occult_values <- values(occult_resampled)
          amplified_values <- mcmapply(amplify_dbzh, dbzh_values, occult_values,
                                       MoreArgs = list(threshold = 0.5, wavelength = 10.71),
                                       mc.cores =  parallel::detectCores())
          values(pvol_raster) <- amplified_values
        }

        mask <- occult_resampled > threshold
        pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
        values <- values(pvol_raster_masked)

        pvol_raster_masked_matrix <- matrix(values, nrow=nrow(pvol_raster_masked), byrow=TRUE)
        class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
        my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix
        attributes(my_pvol$scans[[pvol_idx]]$params$DBZH) <- param_attrs

      } else {

        log_info('applying mask above {threshold} without resampling')

        if (amplify==TRUE) {

        dbzh_values <- values(pvol_raster)
        occult_values <- values(occult_raster)
        amplified_values <- mcmapply(amplify_dbzh, dbzh_values, occult_values,
                                     MoreArgs = list(threshold = 0.5, wavelength = 10.71),
                                     mc.cores =  parallel::detectCores())
        values(pvol_raster) <- amplified_values

        }

        mask <- occult_raster > threshold
        pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
        values <- values(pvol_raster_masked)
        pvol_raster_masked_matrix <- matrix(values, nrow=nrow(pvol_raster_masked), byrow=TRUE)
        class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
        my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix
        attributes(my_pvol$scans[[pvol_idx]]$params$DBZH) <- param_attrs
      }
    }
  }
  return(my_pvol)
}

process_dem_correction <- function(vp, ppi, dem_raster) {
  # Extract and preprocess the vertical profile data
  data <- vp$data
  data <- na.omit(data[c("dens", "height")])
  data$dens[is.na(data$dens)] <- 0
  data$cum_dens <- cumsum(data$dens)
  max_cum_dens <- max(data$cum_dens, na.rm = TRUE)
  max_height <- max(data$height, na.rm = TRUE)
  data$norm_cum_dens <- data$cum_dens / max_cum_dens
  interp_fun <- approxfun(data$height, data$norm_cum_dens, method = "linear")

  #get bbox from ppi and set the correct order for min and max values
  bbox_values <- c(xmin = ppi$geo$bbox[1], ymin = ppi$geo$bbox[2],
                   xmax = ppi$geo$bbox[3], ymax = ppi$geo$bbox[4])

  #convert bounding box to sfc object and transform CRS
  bbox <- st_as_sfc(st_bbox(bbox_values, crs = 4326))
  crs_albers <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  bbox_transformed <- st_transform(bbox, crs_albers)
  bbox_poly <- as(bbox_transformed, "Spatial")


  #crop the raster using the transformed bounding box
  cropped_raster <- crop(dem_raster, bbox_poly)

  C <- function(h, max_height) {
    interpolated_value <- 1 - interp_fun(h)
    result <- numeric(length(h))
    result[h > max_height] <- 1
    result[is.na(interpolated_value)] <- 1
    result[!is.na(interpolated_value) & h <= max_height] <- interpolated_value[!is.na(interpolated_value) & h <= max_height]
    return(result)
  }

  #apply the function C(h) to the raster values
  values(cropped_raster) <- C(values(cropped_raster), max_height)

  return(rast(cropped_raster))
}

resample_to_ppi <- function(cropped_raster, ppi) {
  target_rows <- dim(ppi)[2]
  target_cols <- dim(ppi)[3]

  #create a template raster that matches ppi dims
  extent_vals <- ext(cropped_raster)
  target_raster <- rast(nrows = target_rows, ncols = target_cols, ext = extent_vals, vals = NA)

  #resample
  resampled_raster <- resample(cropped_raster, target_raster, method = "bilinear")
  return(resampled_raster)
}


dem_raster = raster::brick("gis/elev/elev48i0100a.tif")
subtract_mask= terra::rast('occult/subtraction/KRIW_top_quartile_raster.tif')
my_occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
my_occult_scan <- my_occult$scans[[6]] #0.5 deg
plot(project_as_ppi(my_occult_scan, range_max = 150000), zlim=c(0,1))


pvol_files = list.files('pvol/KCYS/', full.names = TRUE)
my_pvol_file = pvol_files[length(pvol_files)-rand_idx]
my_vp = read_vpfiles(paste0('vp/KRIW/', paste0(basename(pvol_files[length(pvol_files)-rand_idx]),'.h5')))

#clean all scans in a pvol with occultation
my_clean_pvol = clean_pvol(my_pvol_file, my_occult, threshold=0.5, method='occultation', amplify=TRUE)
my_ppi = integrate_to_ppi(my_clean_pvol, my_vp, res = 1000, xlim = c(-150000, 150000), ylim = c(-150000, 150000))
original_plot = plot(my_ppi, param = "VIR") +
  ggtitle(paste("KRIW", as.Date(my_ppi$datetime), sep=" ")) +
  theme(plot.title = element_text(size = 10),
        legend.position = "none")

hc_coeff_raster <- process_dem_correction(my_vp, my_ppi, dem_raster)
hc_coeff <- resample_to_ppi(hc_coeff_raster, my_ppi)

my_corrected_ppi = my_ppi
my_corrected_ppi$data$VIR <-  my_ppi$data$VIR * as.vector(values(hc_coeff))
corrected_plot = plot(my_corrected_ppi, param = "VIR") +
  ggtitle("corrected") +
  theme(plot.title = element_text(size = 10),
        legend.position = "none")

par(mfrow=c(1, 2))
plot(my_ppi, param = "VIR", main = "Original PPI")
plot(my_corrected_ppi, param = "VIR", main = "Corrected PPI")

filename = paste("KRIW", as.Date(my_ppi$datetime), sep="-")

options(scipen = 0)
png(filename, width = 2000, height = 1000, res = 300)
grid.arrange(original_plot, corrected_plot, ncol = 2)
dev.off()



#requires my_occult, dem_raster

process_radar_file <- function(file_path) {
  print(paste("Processing file:", file_path))

  radar = substr(file_path, start = 1, stop = 4)
  my_pvol_file = paste0('pvol/', radar, '/', file_path)
  vp_filename <- paste0('vp/', radar, '/', tools::file_path_sans_ext(basename(file_path)), '.h5')

  if (!file.exists(vp_filename)) {
    print(paste("VP file not found:", vp_filename))
    return(NULL)
  }

  my_vp <- read_vpfiles(vp_filename)

  if (is.null(my_vp)) {
    print("Failed to read VP file.")
    return(NULL)
  }

  if (sum(!is.na(my_vp$data$dens)) == 0) {
    print("No valid data points in VP file.")
    return(NULL)
  }

  print("Starting radar file processing...")

  tryCatch({

    my_orig_pvol = read_pvolfile(my_pvol_file)

    my_clean_pvol = clean_pvol(my_pvol_file, my_occult, threshold=0.5, method='occultation', amplify=TRUE)
    my_ppi = integrate_to_ppi(my_clean_pvol, my_vp, res = 500, xlim = c(-150000, 150000), ylim = c(-150000, 150000))

    hc_coeff_raster <- process_dem_correction(my_vp, my_ppi, dem_raster)
    hc_coeff <- resample_to_ppi(hc_coeff_raster, my_ppi)

    my_corrected_ppi = my_ppi
    my_corrected_ppi$data$VIR <- my_ppi$data$VIR * as.vector(values(hc_coeff))
    my_corrected_ppi$data$VID <- my_ppi$data$VID * as.vector(values(hc_coeff))

    idx = my_corrected_ppi$data$R > 25
    my_corrected_ppi$data$VID[idx] <- NA
    my_corrected_ppi$data$VIR[idx] <- NA

    params <- c("RHOHV", "VRADH", "DBZH", "CELL")
    for (param in params) {
      source_pvol <- if (param == "CELL") my_clean_pvol else my_orig_pvol
      data <- project_parameter(source_pvol, param, 150000, 500)
      my_corrected_ppi <- add_data_to_ppi(my_corrected_ppi, data, param)
    }

  }, error = function(e) {
    print(paste("Error in processing file:", file_path, "Error message:", e$message))
  })

  if (exists("my_corrected_ppi")) {
    return(my_corrected_ppi)
  } else {
    return(NULL)
  }
}



for (j in 21:21) {
    for (i in 1:6) {
    print(paste("j:", j, "i:", i))

      radar = 'KCYS'
      week = j
      hours_after = i

      closest_records = read.csv(paste0('closest_records_', radar, '.csv'))
      closest_records$timestamp <- as.POSIXct(closest_records$timestamp, format = "%Y-%m-%d %H:%M:%S")
      closest_records$week <- week(closest_records$timestamp)
      closest_records = closest_records %>%
        dplyr::select(year, week, closest_record, hours_after_sunset)  %>%
        dplyr::filter(hours_after_sunset == hours_after) %>%
        dplyr::rename(closest_file = closest_record)

      all_files = na.omit(closest_records)
      all_files = all_files[all_files$year >= 2013 & all_files$week == week,]
      files = all_files$closest_file

      numCores <- 12
      cl <- makeCluster(numCores)
      clusterExport(cl, c("read_vpts", "integrate_profile",
                          "read_vpfiles", "resample_to_ppi", "add_data_to_ppi",
                          "amplify_dbzh", "rast", "project_parameter",
                          "clean_pvol", "find_proximal_occult_scan_idx",
                          "get_elevation_angles", "clean_weather",
                          "process_dem_correction", "apply_mistnet",
                          "my_occult", "process_radar_file", "dem_raster"))

      clusterEvalQ(cl, {
        library(dplyr)
        library(bioRad)
        library(vol2birdR)
        library(terra)
        library(raster)
        library(logger)
        library(parallel)
        library(ggplot2)
        library(sf)
        library(gridExtra)
      })

      my_occult= read_pvolfile(paste0('occult/h5-150km/', radar, '/', radar, '.h5'), param='all')
      my_corrected_ppi_list = parLapply(cl, files, process_radar_file)

      clean_list <- Filter(Negate(is.null), my_corrected_ppi_list)
      saveRDS(clean_list, file = paste0("composite_ppi/rds/",radar, "_week", week, "_", hours_after, "_V3_25_corrected_ppi_list.Rds"))

      my_composite_ppi_mean = composite_ppi(clean_list, nx = 300, ny = 300, method="mean")

      saveRDS(my_composite_ppi_mean, file = paste0("composite_ppi/ppi/", radar, "_week", week, "_", hours_after, "_V3_25_mean_composite_ppi.Rds"))

      #png(filename = paste0("composite_ppi/plot/", radar, "_week", week, "_", hours_after, "_V3_25_mean_composite_ppi_plot.png"), width = 1600, height = 1200, res = 300)
      #plot(my_composite_ppi_mean, 'VID') + labs(title=radar, subtitle=paste('week', week, 'mean', sep = ' '))
      #dev.off()

      stopCluster(cl)
      rm(cl)
      gc()
  }
}



add_azimuth_degrees <- function(spatial_grid) {
  if (!inherits(spatial_grid, "SpatialGridDataFrame")) {
    stop("Input must be a SpatialGridDataFrame.")
  }

  coords <- coordinates(spatial_grid)
  azimuth_radians <- atan2(coords[, 2], coords[, 1])
  azimuth_degrees <- (azimuth_radians * 180 / pi)
  azimuth_degrees <- (90 - azimuth_degrees) %% 360
  azimuth_degrees <- ifelse(azimuth_degrees < 0, azimuth_degrees + 360, azimuth_degrees)
  spatial_grid@data$azimuth <- azimuth_degrees
  return(spatial_grid)
}


add_azimuth_degrees <- function(spatial_grid, lon, lat, reference_angle = 0) {
  if (!inherits(spatial_grid, "SpatialGridDataFrame")) {
    stop("Input must be a SpatialGridDataFrame.")
  }

  # Convert input degrees to radians
  lon_center = lon * pi / 180
  lat_center = lat * pi / 180

  coords <- coordinates(spatial_grid)
  lon_coords = coords[, 1] * pi / 180
  lat_coords = coords[, 2] * pi / 180

  azimuth_radians <- atan2(sin(lon_coords - lon_center),
                           cos(lat_center) * tan(lat_coords) - sin(lat_center) * cos(lon_coords - lon_center))

  azimuth_degrees <- azimuth_radians * 180 / pi
  azimuth_degrees <- ifelse(azimuth_degrees < 0, azimuth_degrees + 360, azimuth_degrees)
  azimuth_degrees <- (azimuth_degrees - reference_angle + 360) %% 360
  spatial_grid@data$azimuth <- azimuth_degrees

  return(spatial_grid)
}




# -------------------------------------------------------------------------

vp_stat_df <- na.omit(vp_stat_df)

vp_stat_df$radians <- vp_stat_df$dd * pi / 180

weighted_cos <- sum(cos(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)
weighted_sin <- sum(sin(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)

mean_angle_radians <- atan2(weighted_sin, weighted_cos)

mean_angle_degrees <- mean_angle_radians * 180 / pi
mean_angle_degrees <- ifelse(mean_angle_degrees < 0, mean_angle_degrees + 360, mean_angle_degrees)

mean_angle_degrees


vp_stat_df <- data.frame(
  vid = unlist(vp_stat_list$vid),
  dd = unlist(vp_stat_list$dd)
)

vp_stat_df$dd_adjusted <- ifelse(vp_stat_df$dd < 180, vp_stat_df$dd + 360, vp_stat_df$dd)

# Remove rows with NA in either column before computing the weighted mean
vp_stat_df <- na.omit(vp_stat_df)

# Compute the weighted mean of dd, weighted by vid
weighted_mean_dd <- weighted.mean(vp_stat_df$dd_adjusted, w = vp_stat_df$vid, na.rm = TRUE)

# Display the computed weighted mean
print(weighted_mean_dd)





center_lon <- mean_composite_ppi$geo$lon[1]
center_lat < -mean_composite_ppi$geo$lat[1]

angle_degrees <- weighted_mean_dd
angle_degrees <- (90 - angle_degrees) %% 360
angle_degrees <- ifelse(angle_degrees < 0, angle_degrees + 360, angle_degrees)

distance_km <- 150  # Distance in kilometers

# Prepare data frame for ggplot
center_data <- data.frame(
  x = center_lon,
  y = center_lat,
  angle = angle_degrees * pi / 180,
  radius = distance_km / 111.32  # Approx conversion factor from degrees to kilometers
)

# ----------------------------------------------------------------------


mean_composite_ppi$data <- add_azimuth_degrees(mean_composite_ppi$data,
                                               lon=mean_composite_ppi$geo$lon[1],
                                               lat=mean_composite_ppi$geo$lat[1],
                                               reference_angle = weighted_mean_dd)

data_for_plot <- as.data.frame(mean_composite_ppi$data)
data_for_gam <- na.omit(data_for_plot[, c("VID", "azimuth")])

gam_fit <- gam(VID ~ s(azimuth, bs = "cc"), data = data_for_gam)
summary(gam_fit)
data_for_plot$fitted_values <- predict(gam_fit, newdata = data_for_plot, type = "response")


# Create the plot
plot <- ggplot(data_for_plot, aes(x = azimuth, y = VID)) +
  geom_point(aes(color = "Observed"), alpha = 0.01, size = 0.5) +
  geom_line(aes(y = fitted_values, color = "Fitted"), size = 1) +
  labs(x = "Azimuth (degrees)", y = "VID", title = "GAM fit VID ~ Azimuth") +
  scale_color_manual(name = "", values = c("Observed" = "blue", "Fitted" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 50)



plot(mean_composite_ppi) +
  geom_spoke(data=center_data, aes(x=x, y=y, angle = angle, radius = radius),
             color = "black", size = 0.5, arrow = arrow(length = unit(0.2, "cm"))) +
  ggtitle('KYCS week 19 hour 3')



library(ggplot2)

# Create a data frame with x values ranging from 0 to 360
data <- data.frame(x = seq(0, 360, by = 1))

# Calculate y using the sine function to oscillate between 0.5 and 1.5
data$y <- 1 + 0.5 * sin(2 * data$x * pi / 180)  # Amplitude 0.5, midline 1

# Plotting the data
ggplot(data, aes(x = x, y = y)) +
  geom_line(col="red") +  # Plot the curve
  labs(x = "azimuth", y = "coeff") +
  ggtitle("VID ") +
  theme_minimal()






library(png)
library(grid)
library(gridExtra)
library(magrittr)


img1 <- readPNG("KCYS_week18_3_V1_mean_composite_ppi_plot.png")
img2 <- readPNG("KCYS_week18_3_V2_mean_composite_ppi_plot.png")
img3 <- readPNG("KCYS_week18_3_V3_mean_composite_ppi_plot.png")
img4 <- readPNG("KCYS_week18_3_V3_50_mean_composite_ppi_plot.png")
img5 <- readPNG("KCYS_week18_3_V3_25_mean_composite_ppi_plot.png")
img6 <- readPNG("KCYS_week18_3_V3_10_mean_composite_ppi_plot.png")

# Create grobs and annotate with text in the top right corner
annotate_img <- function(img, label) {
  grob <- rasterGrob(img, interpolate = TRUE)
  txt <- textGrob(label, x = 0.95, y = 0.95, hjust = 1, vjust = 1, gp = gpar(col = "white", fontface = "bold"))
  grobTree(grob, txt)
}

grob1 <- annotate_img(img1, "V1")
grob2 <- annotate_img(img2, "V2")
grob3 <- annotate_img(img3, "V2 R > 100")
grob4 <- annotate_img(img4, "V2 R > 50")
grob5 <- annotate_img(img5, "V2 R > 25")
grob6 <- annotate_img(img6, "V2 R > 10")

final_grob <- arrangeGrob(grob1, grob2, grob3, grob4, grob5, grob6,  ncol = 2, nrow = 3)
png("combined_image.png", width = 1600, height = 1600, res = 300)  # Adjust size as needed
grid.draw(final_grob)
dev.off()


# Get a list of all .Rds files
files <- list.files(pattern = "\\.Rds$", recursive = TRUE, full.names = TRUE)

process_file <- function(file_path) {
  my_composite_ppi_mean <- readRDS(file_path)
  r <- rast(my_composite_ppi_mean$data)
  r_rounded <- app(r, round, digits = 2)
  polys <- as.polygons(r_rounded, round = FALSE, aggregate = FALSE) |> st_as_sf()
  output_file <- gsub("Rds$", "geojson", file_path)
  st_write(polys, output_file, driver = "GeoJSON", quiet = TRUE)
  return(paste("Processed and saved:", output_file))
}

results <- mclapply(files, process_file, mc.cores = 10)
print(results)


files <- list.files(pattern = "\\.Rds$", recursive = TRUE, full.names = TRUE)


process_file <- function(file_path) {
  my_composite_ppi_mean <- readRDS(file_path)

  sgdf <- my_composite_ppi_mean$data
  vid_layer <- rast(sgdf["VID"])

  name_parts <- strsplit(file_path, '_')
  week <- gsub('week', '' , name_parts[[1]][2])
  hour <- name_parts[[1]][3]
  print(paste0(week,hour))

  layer_name <- paste("week", gsub("w", "", week), "h", hour, sep="_")
  names(vid_layer) <- layer_name
  return(vid_layer)
}

all_layers <- lapply(files, process_file)
combined_raster <- rast(all_layers)
r_rounded <- app(combined_raster, round, digits = 2)
polys <- as.polygons(r_rounded, round = FALSE, aggregate = FALSE) |> st_as_sf()
st_write(polys, 'vid.geojson', driver = "GeoJSON", quiet = TRUE)








i = 0
for (file_path in files) {
  i <- i + 1
  data <- readRDS(file_path)
  print(file_path)

  current_min <- min(data[[1]]$VID, na.rm = TRUE)
  current_max <- max(data[[1]]$VID, na.rm = TRUE)

  global_min[[i]] <- current_min
  global_max[[i]] <- current_max

}




# Get the list of RDS files
file_names <- list.files(pattern = "\\.Rds$", full.names = TRUE)

# Directory where the plots will be saved
output_dir <- "raster/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through each RDS file
for (file_name in file_names) {
  # Read the RDS file
  my_composite_ppi_mean <- readRDS(file_name)

  vid_layer <- raster(my_composite_ppi_mean$data, layer="VID")

  vid_vals <- values(vid_layer)
  vid_vals[vid_vals == 0] <- NA
  values(vid_layer) <- vid_vals

  parts <- strsplit(basename(file_name), "_")[[1]]
  week <- parts[2]
  hour <- sub("V3_25_mean_composite_ppi.Rds", "", parts[3])

  # Prepare the filename for the output image
  tif_filename <- paste0(output_dir, gsub("\\.Rds$", ".tif", basename(file_name)))
  writeRaster(vid_layer, filename=tif_filename, format="GTiff", overwrite=TRUE)

    # Plot the data
  #p = plot(my_composite_ppi_mean, 'VID') + labs(title = paste("KCYS", week, "-", "hour", hour))
  #ggsave(png_filename, plot = p, width = 16, height = 12, units = "in", dpi = 300)

}

my_composite_ppi_mean <- readRDS(file_names[9])

ggplot() +
  annotation_map_tile(type = "osm") +
  layer_spatial(my_composite_ppi_mean$data, aes(alpha=0.5, color = VID)) +  # Ensure this is color
  scale_color_viridis_c()

rlayer = raster(my_composite_ppi_mean$data)
writeRaster(rlayer, filename="output_raster.tif", format="GTiff", overwrite=TRUE)







#KRIW
#20130430_051424

my_pvol = read_pvolfile("/Users/at744/clo/bigbird/wgfd/pvol/KRIW/KRIW20130430_051424_V06.gz")

sum(is.na(my_pvol$scans[[1]]$params$DBZH) & !is.nan(my_pvol$scans[[1]]$params$DBZH))
sum(is.nan(my_pvol$scans[[1]]$params$DBZH))














#my_composite_ppi_max = composite_ppi(clean_list, nx = 300, ny = 300, method="max")
#png(filename = paste0("KCYS_week", week, "_", hours_after, "_max_composite_ppi_plot.png"), width = 1600, height = 1200, res = 300)
#plot(my_composite_ppi_max, 'VID') + labs(title='KCYS', subtitle=paste('week', week, 'max', sep = ' '))
#dev.off()



#MISC
-------------
#my_corrected_spat_raster <- rast(my_corrected_ppi$data)
#writeRaster(my_corrected_spat_raster, filename = paste0('tif/', radar, '/', basefile, '.tif'), format = "GTiff", overwrite=TRUE)

#original_plot = plot(my_ppi, param = "VIR") +
#  ggtitle(paste(radar, as.Date(my_ppi$datetime), sep=" ")) +
#  theme(plot.title = element_text(size = 10), legend.position = "none")

#corrected_plot = plot(my_corrected_ppi, param = "VIR") +
#  ggtitle("corrected") +
#  theme(plot.title = element_text(size = 10), legend.position = "none")

# Save the plots
#filename <- paste0("compare/KRIW/KRIW-", as.Date(my_ppi$datetime), ".png")
#options(scipen = 0)
#png(filename, width = 2000, height = 1000, res = 300)
#grid.arrange(original_plot, corrected_plot, ncol = 2)
#dev.off()

#my_corrected_spat_raster <- rast(my_corrected_ppi$data)
#out_filename <- paste0('tif/', radar, '/', tools::file_path_sans_ext(basename(file_path)), '.tif')

my_corrected_ppi

print(paste("Writing raster to file:", out_filename))
terra::writeRaster(my_corrected_spat_raster, filename = out_filename, overwrite=TRUE)
print("File written successfully.")


process_radar_file <- function(file_path) {
  my_pvol_file = file_path
  vp_filename <- paste0('vp/', radar, '/', tools::file_path_sans_ext(basename(file_path)), '.h5')
  my_vp <- read_vpfiles(vp_filename)
  basefile = tools::file_path_sans_ext(basename(file_path))

  if (sum(!is.na(my_vp$data$dens))){
    tryCatch({
      # Perform radar file processing
      my_clean_pvol = clean_pvol(my_pvol_file, my_occult, threshold=0.5, method='occultation', amplify=TRUE)
      my_ppi = integrate_to_ppi(my_clean_pvol, my_vp, res = 1000, xlim = c(-150000, 150000), ylim = c(-150000, 150000))

      hc_coeff_raster <- process_dem_correction(my_vp, my_ppi, dem_raster)
      hc_coeff <- resample_to_ppi(hc_coeff_raster, my_ppi)

      my_corrected_ppi = my_ppi
      my_corrected_ppi$data$VIR <- my_ppi$data$VIR * as.vector(values(hc_coeff))


    }, error = function(e) {
      cat("Error in processing file:", file_path, "\nError message:", e$message, "\n")
    }, finally = {
      cat("Finished processing, whether success or error, for file:", file_path, "\n")
    })
  } else {
    cat("No valid data points in VP file for:", file_path, "\n")
  }
}


project_parameter <- function(pvol_data, param, range_max, grid_size) {
  param_data = pvol_data$scans[[1]]$params[[param]]
  return(project_as_ppi(param_data, range_max = range_max, grid_size = grid_size)$data[[1]])
}

add_data_to_ppi <- function(ppi, data, param_name) {
  ppi$data[[param_name]] <- data
  return(ppi)
}


ppi_raster <- rast(my_composite_ppi_mean$data)
hill_transformed <- project(hill, crs(ppi_raster))

transformed_hill <- rast("transformed_hillshade.tif")

plot(transformed_hill, col = grey(0:100/100), main = "Hillshade with PPI Overlay")
plot(ppi_raster, add = TRUE, col = rainbow(100, alpha = 0.5))  # Adjust alpha for transparency


hill_df <- as.data.frame(transformed_hill, xy = TRUE, na.rm = TRUE)
colnames(hill_df) <- c("lon", "lat", "hillshade")

ppi_df <- as.data.frame(ppi_raster, xy = TRUE, na.rm = TRUE)
colnames(ppi_df) <- c("lon", "lat", "ppi")


plot(ppi) + ggtitle('hello')

hill_df <- as.data.frame(transformed_hill, xy = TRUE, na.rm = TRUE)
colnames(hill_df) <- c("lon", "lat", "hillshade")



vp_stat_list <- list()

for (i in 1:length(files)) {
  tryCatch({
    # Attempt to process each file
    print(i)
    filename <- paste0('vp/KCYS/', tools::file_path_sans_ext(basename(files[i])), '.h5')
    my_vp <- read_vpfiles(filename)

    vp_stat_list$vid[[i]] <- integrate_profile(my_vp)$vid
    vp_stat_list$dd[[i]] <- integrate_profile(my_vp)$dd

  }, error = function(e) {
    # Error handling
    cat("Error in processing file", i, ": ", e$message, "\n")
  })
}

as.data.frame(vp_stat_list)

#KABR
#KAKQ
#KABX
#KBGM


files <- list.files(pattern = "KLIX.*\\.csv$", recursive = TRUE)

for (i in 1:length(files)){
print(files[i])
bioRad::read_vpts(files[i])
}

process_radar_files <- function(radar) {
  cat("Processing files for radar:", radar, "\n")

  # Define the pattern and list files
  pattern <- paste0(radar, ".*\\.csv$")
  files <- list.files(pattern = pattern, recursive = TRUE)

  # Initialize a vector to collect filenames that cause errors
  error_files <- c()

  # Process each file
  for (file_path in files) {
    result <- tryCatch({
      # Attempt to read the file
      vpts_data <- bioRad::read_vpts(file_path)
      NULL  # Return NULL if no error occurs
    }, error = function(e) {
      return(file_path)  # Return the problematic file path if an error occurs
    })

    if (!is.null(result)) {
      error_files <- c(error_files, result)
    }
  }

  cat("Completed processing for radar:", radar, "\n")
  return(list(radar = radar, error_files = error_files))
}


radars <- radarInfo$radar

# Set up a cluster
num_cores <- 10  # Reserve one core for system operations
cl <- makeCluster(num_cores)
clusterExport(cl, "read_vpts")  # Export necessary variables
clusterEvalQ(cl, {
  library(bioRad)  # Load necessary libraries in each cluster node
})

# Use parLapply to process files in parallel
results <- parLapply(cl, radars, process_radar_files)










process_radar_data <- function(radar) {
       files <- list.files(pattern = radar, recursive = TRUE)

         if (length(files) > 0) {
               vpts_data <- bioRad::read_vpts(files)
               saveRDS(vpts_data, file = paste0(radar, ".Rds"))
           } else {
                 cat("No files found for radar: ", radar, "\n")
             }
   }

num_cores <- 10
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("list.files", "read_vpts"))

parLapply(cl, radarInfo$radar, process_radar_data)


bioRad::read_vpts(files[1:3])

library(bioRad)

for (j in 1:length(radarInfo$radar)){

  radar = radarInfo$radar[j]
  files <- list.files(pattern = radar, recursive = TRUE)

    continue_reading <- TRUE
    for (i in 245:250) {
      if (continue_reading) {
        tryCatch({
          vpts_data <- bioRad::read_vpts(files[244:i])
          cat("Successfully read ", i, "files.\n")
        }, error = function(e) {
          # Handle errors: print message and update flag to stop the loop
          cat("Error encountered while reading ", i, "files:\n")
          cat("Problem file:", files[i], "\n")
          cat("Error message:", e$message, "\n")
          continue_reading <- FALSE  # Set flag to FALSE to break the loop
        })
      } else {
        break  # Exit the loop if an error has been encountered
      }
    }
}

bioRad::read_vpts(files[240])
