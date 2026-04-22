library(terra)
library(raster)
library(bioRad)
library(raster)
library(logger)
library(lubridate)
library(ggplot2)
library(dplyr)
library(mgcv)
library(units)
library(gridExtra)
library(pracma)

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


project_parameter <- function(pvol_data, param, range_max, grid_size) {
  param_data = pvol_data$scans[[1]]$params[[param]]
  return(project_as_ppi(param_data, range_max = range_max, grid_size = grid_size)$data[[1]])
}

add_data_to_ppi <- function(ppi, data, param_name) {
  ppi$data[[param_name]] <- data
  return(ppi)
}



dem_raster = raster::brick("gis/elev/elev48i0100a.tif")
my_occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
my_occult_scan <- my_occult$scans[[6]] #0.5 deg
plot(project_as_ppi(my_occult_scan, range_max = 150000), zlim=c(0,1))



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

closest_records = read.csv(paste0('closest_records_all.csv'))


  for (j in c(18:22, 34:40)) {
    for (i in c(3:3)) {
    print(paste("j:", j, "i:", i))

      radar = 'KCBX'
      week = j
      hours_after = i

      closest_records <- closest_records[closest_records$radar == radar,]
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

      my_occult= read_pvolfile(paste0('occult/h5-150km/', radar, '/', radar, '.h5'), param='all')

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

      my_corrected_ppi_list = parLapply(cl, files, process_radar_file)

      clean_list <- Filter(Negate(is.null), my_corrected_ppi_list)
      #saveRDS(clean_list, file = paste0("composite_ppi/rds/",radar, "_week", week, "_", hours_after, "_V3_25_corrected_ppi_list.Rds"))

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

  #convert input degrees to radians
  lon_center = lon * pi / 180
  lat_center = lat * pi / 180

  coords <- coordinates(spatial_grid)
  lon_coords = coords[, 1] * pi / 180
  lat_coords = coords[, 2] * pi / 180

  azimuth_radians <- atan2(sin(lon_coords - lon_center),
                           cos(lat_center) * tan(lat_coords) - sin(lat_center) * cos(lon_coords - lon_center))

  azimuth_degrees <- azimuth_radians * 180 / pi
  spatial_grid@data$azimuth <- ifelse(azimuth_degrees < 0, azimuth_degrees + 360, azimuth_degrees)

  azimuth_degrees <- ifelse(azimuth_degrees < 0, azimuth_degrees + 360, azimuth_degrees)
  azimuth_degrees <- (azimuth_degrees - reference_angle + 360) %% 360
  spatial_grid@data$relative_azimuth <- azimuth_degrees

  return(spatial_grid)
}


# -------------------------------------------------------------------------
# Add relative azimuth degrees to composite ppi data

mean_composite_ppi <- readRDS(paste0('composite_ppi/ppi/', names(dataframe_list)))


mean_composite_ppi + theme(legend.position = "none")

mean_composite_ppi$data <- add_azimuth_degrees(mean_composite_ppi$data
                                               ,lon=my_vp$attributes$where$lon
                                               ,lat=my_vp$attributes$where$lat
                                               ,reference_angle=vp_df_list[[1]]$stats$mean_angle_degrees)
library(gridExtra)
library(ggplot2)

modified_composite_ppis <- list()

for (composite_name in names(dataframe_list)) {

  mean_composite_ppi <- readRDS(composite_name)

  vp_stat_df <- dataframe_list[[composite_name]]$data
  vp_stat_df$radians <- vp_stat_df$dd * pi / 180
  weighted_cos <- sum(cos(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)
  weighted_sin <- sum(sin(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)
  mean_angle_radians <- atan2(weighted_sin, weighted_cos)
  angle_degrees <- mean_angle_radians * 180 / pi
  angle_degrees <- ifelse(angle_degrees < 0, angle_degrees + 360, angle_degrees)

  print(angle_degrees)

  grid_data <- mean_composite_ppi$data
  lon <- mean_composite_ppi$geo$lon[1]
  lat <- mean_composite_ppi$geo$lat[1]

  #apply the transformation to the grid_data
  modified_grid_data <- add_azimuth_degrees(
    grid_data,
    lon = lon,
    lat = lat,
    reference_angle = angle_degrees
  )

  mean_composite_ppi$data <- modified_grid_data
  modified_composite_ppis[[composite_name]] <- mean_composite_ppi
}



resample_dem_raster <- function(cropped_raster, ppi) {
  target_rows <- dim(ppi)[2]
  target_cols <- dim(ppi)[3]

  #create a template raster
  extent_vals <- ext(cropped_raster)
  target_raster <- rast(nrows = target_rows, ncols = target_cols, ext = extent_vals, vals = NA)

  #resample
  resampled_raster <- resample(cropped_raster, target_raster, method = "bilinear")
  return(resampled_raster)
}



deg2num<-function(lat_deg, lon_deg, zoom){
  lat_rad <- lat_deg * pi /180
  n <- 2.0 ^ zoom
  xtile <- floor((lon_deg + 180.0) / 360.0 * n)
  ytile = floor((1.0 - log(tan(lat_rad) + (1 / cos(lat_rad))) / pi) / 2.0 * n)
  return( c(xtile, ytile))
  #  return(paste(paste("https://tile.openstreetmap.org", zoom, xtile, ytile, sep="/"),".png",sep=""))
}



#### Experimental 
process_pvol <- function(file_path) {
  tryCatch({

    #output_dir_raster <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/raster/"
    #output_dir_image <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/image/"
    #filename = basename(tools::file_path_sans_ext(file_path))
    my_pvol <- apply_mistnet(file_path) %>% clean_pvol()
    #return(my_pvol)


    my_clean_pvol = apply_mistnet("/Users/at744/clo/bigbird/wgfd/pvol/KRIW/KRIW20190919_041656_V06") %>% clean_pvol()
    my_param <- get_param(my_clean_pvol$scans[[1]], 'DBZH')
    my_ppi <- project_as_ppi(my_param, range_max = 150000)
    my_new_ppi = my_ppi

    subtract_latent_dbz <- function(ppi_dbz, upsampled_raster) {

      eta_values <- sapply(as.matrix(ppi_dbz), dbz_to_eta_safe)
      upsampled_matrix <- matrix(upsampled_raster)
      upsampled_matrix[is.na(upsampled_matrix) | is.nan(upsampled_matrix)] <- 0

      result_matrix <- ifelse(is.na(eta_values), NA, eta_values - upsampled_matrix)
      result_matrix[result_matrix < 0] <- 0
      result_dbz <- sapply(result_matrix, eta_to_dbz_safe)
      return(result_dbz)
    }

    my_new_ppi$data$DBZH <- subtract_latent_dbz(my_new_ppi$data, upsampled_raster)


    raster_df = as.data.frame(my_raster, xy=TRUE) %>% na.omit()
    p <- (
      ggplot(data = raster_df) +
        geom_raster(aes(x = x, y = y, fill=DBZH)) +
        scale_fill_viridis_c(limit = c(-35, 55)) +
        theme_void() +
        theme(
          panel.background = element_rect(fill = 'black'),
          plot.background = element_rect(fill = "black"),
          legend.position = "bottom",
          legend.text = element_text(color = "#E4EBF7"),  # Set legend text color
          legend.title = element_text(color = "#E4EBF7"),  # Set legend title color
        )
    )
    
   ggsave(filename = paste0(output_dir_image, filename, ".png"), plot = p, width = 10, height = 8, dpi = 300)

  }, error = function(e) {
    message("Error processing file: ", file_path, "\nError message: ", e$message)
    return(NULL)  # Return NULL for this file to indicate the error
  })
}


numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("process_pvol", "clean_pvol", "apply_mistnet",
                              "get_param", "project_as_ppi"))
clusterEvalQ(cl, {
  library(dplyr)
  library(bioRad)
  library(raster)
  library(magrittr)
  library(ggplot2)
})

file_path <- files[1]

parLapply(cl, files, process_pvol)
lapply(files[1], process_pvol)

######
#write corrected vid to geojson


radar = 'KRIW'
setwd("~/clo/bigbird/wgfd/composite_ppi/composite/KRIW/corrected")
pattern <- paste0(".*", radar, ".*\\.Rds$")
files <- list.files(pattern = pattern, recursive = FALSE, full.names = TRUE)


process_file <- function(file_path) {
  my_composite_ppi_mean <- readRDS(file_path)

  sgdf <- my_composite_ppi_mean$data
  vid_layer <- rast(sgdf["adjusted_VID"])

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
st_write(polys, paste0('~/clo/bigbird/wgfd/composite_ppi/composite/', radar, '/corrected/vid.geojson'), driver = "GeoJSON", quiet = TRUE)




# Get a list of all .Rds files
files <- list.files(pattern = "\\.Rds$", recursive = TRUE, full.names = TRUE)

desired_bands <- c(
  "VIR",
  "VID",
  "R",
  "overlap",
  "eta_sum",
  "eta_sum_expected",
  "RHOHV",
  "VRADH",
  "DBZH",
  "CELL",
  "adjusted_VID"
)


##write out polygonized rasters/geojsons
process_file <- function(file_path) {

  tryCatch({

  my_composite_ppi_mean <- readRDS(file_path)

  r <- rast(my_composite_ppi_mean$data)
  r_subset <- r[[desired_bands]]
  r_rounded <- app(r_subset, round, digits = 2)
  output_file <- sub("\\.Rds$", ".geojson", file_path)

  #writeRaster(r_rounded, filename = paste0('raster/', output_file), overwrite = TRUE, filetype = "GTiff")

  polys <- as.polygons(r_rounded, round = FALSE, aggregate = FALSE) |> st_as_sf()
  output_file <- gsub("Rds$", "geojson", file_path)
  st_write(polys, output_file, driver = "GeoJSON", quiet = TRUE)
  return(paste("Processed and saved:", output_file))
  }, error = function(e) {
    # Return error information
    return(list(status = "error", file = file_path, message = e$message))
  })
}

mclapply(files, process_file, mc.cores = 10)


# -------------------------------------------------------------------------
# Add relative azimuth degrees to composite ppi data

mean_composite_ppi <- readRDS(paste0('composite_ppi/ppi/', names(dataframe_list)))


mean_composite_ppi + theme(legend.position = "none")

mean_composite_ppi$data <- add_azimuth_degrees(mean_composite_ppi$data
                                               ,lon=my_vp$attributes$where$lon
                                               ,lat=my_vp$attributes$where$lat
                                               ,reference_angle=vp_df_list[[1]]$stats$mean_angle_degrees)
library(gridExtra)
library(ggplot2)


modified_composite_ppis <- list()

for (composite_name in names(dataframe_list)) {

  mean_composite_ppi <- readRDS(composite_name)

  vp_stat_df <- dataframe_list[[composite_name]]$data
  vp_stat_df$radians <- vp_stat_df$dd * pi / 180
  weighted_cos <- sum(cos(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)
  weighted_sin <- sum(sin(vp_stat_df$radians) * vp_stat_df$vid) / sum(vp_stat_df$vid)
  mean_angle_radians <- atan2(weighted_sin, weighted_cos)
  angle_degrees <- mean_angle_radians * 180 / pi
  angle_degrees <- ifelse(angle_degrees < 0, angle_degrees + 360, angle_degrees)

  print(angle_degrees)

  grid_data <- mean_composite_ppi$data
  lon <- mean_composite_ppi$geo$lon[1]
  lat <- mean_composite_ppi$geo$lat[1]

  #apply the transformation to the grid_data
  modified_grid_data <- add_azimuth_degrees(
    grid_data,
    lon = lon,
    lat = lat,
    reference_angle = angle_degrees
  )

  mean_composite_ppi$data <- modified_grid_data
  modified_composite_ppis[[composite_name]] <- mean_composite_ppi
}




# -------------------------------------------------------------------------
# Get direction dd and vid from each vp and store

# Initialize an empty list to store results
vp_stat_list <- list()

# Iterate over each set of matched VP files
for (composite_name in names(matching_vp_files)) {
  vp_files <- matching_vp_files[[composite_name]]

  # Initialize a list to store the results for this composite file
  vp_stats_for_composite <- list(vid = list(), dd = list())

  # Iterate through each file in the current set of matched files
  for (i in seq_along(vp_files)) {
    tryCatch({
      # Construct the full path for the VP file
      filename <- file.path("vp/KCYS", vp_files[i])
      print(paste("Processing file:", filename))

      # Process the VP file
      my_vp <- read_vpfiles(filename)

      # Calculate and store stats
      vp_stats <- integrate_profile(my_vp)
      vp_stats_for_composite$vid[[i]] <- vp_stats$vid
      vp_stats_for_composite$dd[[i]] <- vp_stats$dd

    }, error = function(e) {
      # Error handling
      cat("Error in processing file", filename, ": ", e$message, "\n")
    })
  }

  # Store the stats for this composite file in the main list
  vp_stat_list[[composite_name]] <- vp_stats_for_composite
}

# Now vp_stat_list is structured with each composite ppi having its own list of vid and dd
print(vp_stat_list)


----------------------------------------------------------------------


vp_stat_list$KCYS_week18_3_V3_25_mean_composite_ppi.Rds


vp_stat_list <- list()

for (i in 1:length(matching_vp_files[1])) {
  tryCatch({
    # Attempt to process each file
    print(i)
    filename <- paste0('vp/KCYS/', tools::file_path_sans_ext(basename(files[i])), '.h5')
    my_vp <- read_vpfiles(filename)

    print(integrate_profile(my_vp)$vid)
    print(integrate_profile(my_vp)$dd)
    #vp_stat_list$vid[[i]] <- integrate_profile(my_vp)$vid
    #vp_stat_list$dd[[i]] <- integrate_profile(my_vp)$dd

  }, error = function(e) {
    # Error handling
    cat("Error in processing file", i, ": ", e$message, "\n")
  })
}

vp_stat_list = as.data.frame(vp_stat_list)

# -------------------------------------------------------------------------
#List to dataframes

dataframe_list <- list()

# Loop over each composite file entry in vp_stat_list
for (composite_name in names(vp_stat_list)) {
  stats <- vp_stat_list[[composite_name]]

  if (length(stats$vid) > 0 && length(stats$dd) > 0) {
    df <- data.frame(
      vid = unlist(stats$vid),  # Unlist to convert list of numbers into a vector
      dd = unlist(stats$dd),    # Same for dd
      stringsAsFactors = FALSE  # Avoid factor conversion
    )

    dataframe_list[[composite_name]] <- df
  } else {
    dataframe_list[[composite_name]] <- data.frame(vid = numeric(0), dd = numeric(0))
  }
}

# Optional: Print one of the dataframes to check structure
if (length(dataframe_list) > 0) {
  print(head(dataframe_list[[1]]))
}



# -------------------------------------------------------------------------
# Weighted means


# Loop over each dataframe in the list
for (name in names(dataframe_list)) {
  # Access the dataframe
  df <- dataframe_list[[name]]

  # Remove rows with NA values
  df <- na.omit(df)

  df$dd_adjusted <- ifelse(df$dd < 180, df$dd + 360, df$dd)

  df$radians <- df$dd_adjusted * pi / 180

  weighted_cos <- sum(cos(df$radians) * df$vid) / sum(df$vid)
  weighted_sin <- sum(sin(df$radians) * df$vid) / sum(df$vid)

  mean_angle_radians <- atan2(weighted_sin, weighted_cos)

  mean_angle_degrees <- mean_angle_radians * 180 / pi
  mean_angle_degrees <- ifelse(mean_angle_degrees < 0, mean_angle_degrees + 360, mean_angle_degrees)

  weighted_mean_dd <- weighted.mean(df$dd_adjusted, w = df$vid, na.rm = TRUE)

  results <- list(
    weighted_mean_dd = weighted_mean_dd,
    mean_angle_degrees = mean_angle_degrees
  )

  dataframe_list[[name]] <- list(data = df, stats = results)
}

# Print results for verification
print(dataframe_list[[1]]$stats)




#Extract discrete threshold boundaries using some eta_mean


library(bioRad)
library(raster) 
library(sp)
library(sf)


station <- "KCYS"
weeks   <- c(18, 19, 20, 21, 36, 37, 38, 39)

my_ppis <- lapply(weeks, function(w) {
  obj_name <- paste0(station, "_week", w, "_3_V3_25_mean_composite_ppi_adjusted")
  get(obj_name, envir = .GlobalEnv) 
})


eta_list <- lapply(my_ppis, function(x) x$data$eta_sum_expected)
eta_mean <- Reduce("+", eta_list) / length(eta_list)

mean_ppi <- my_ppis[[1]]
mean_ppi$data$eta_sum_expected <- eta_mean

thresholds <- c(50) 

lines_list <- list()

for(th in thresholds) {
  z_masked <- ifelse(eta_mean >= th, eta_mean, NA)
  
  r <- raster(mean_ppi$data) 
  values(r) <- z_masked  #
  
  if(all(is.na(z_masked))) {
    next
  }
  poly_sp <- rasterToPolygons(r, 
                              fun = function(x) !is.na(x),
                              dissolve = TRUE)
  poly_sf <- st_as_sf(poly_sp)
  
  poly_union <- st_union(poly_sf)
  boundary_sf <- st_boundary(poly_union)
  
  boundary_sf <- st_sf(
    threshold = th,        
    geometry = boundary_sf
  )
  
  lines_list[[length(lines_list) + 1]] <- boundary_sf
}


all_boundaries_sf <- do.call(rbind, lines_list)
#st_write(all_boundaries_sf, "KCYS_threshold_boundary.geojson", driver = "GeoJSON")


