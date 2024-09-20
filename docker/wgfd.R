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
library(dplyr)
library(parallel)
library(microbenchmark)


# Retrieve environment variables
radar <- Sys.getenv("RADAR")
week <- as.numeric(Sys.getenv("WEEK"))
hours_after <- as.numeric(Sys.getenv("HOURS_AFTER"))

cat("Radar:", radar, "\n")
cat("Week:", week, "\n")
cat("Hours After:", hours_after, "\n")


dem_raster = raster::brick("gis/elev/elev48i0100a.tif")

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

# numCores <- 10
# cl <- makeCluster(numCores)

my_occult= read_pvolfile(paste0('occult/h5-150km/', radar, '/', radar, '.h5'), param='all')


#parallel processing______________________________________________________
# clusterExport(cl, c("read_vpts", "integrate_profile",
#                     "read_vpfiles", "resample_to_ppi", "add_data_to_ppi",
#                     "amplify_dbzh", "rast", "project_parameter",
#                     "clean_pvol", "find_proximal_occult_scan_idx",
#                     "get_elevation_angles", "clean_weather",
#                     "process_dem_correction", "apply_mistnet",
#                     "my_occult", "process_radar_file", "dem_raster"))

# clusterEvalQ(cl, {
#   library(dplyr)
#   library(bioRad)
#   library(vol2birdR)
#   library(terra)
#   library(raster)
#   library(logger)
#   library(parallel)
#   library(ggplot2)
#   library(sf)
#   library(gridExtra)
# })

#my_corrected_ppi_list = parLapply(cl, files, process_radar_file)


#sequential processing______________________________________________________
benchmark_results <- microbenchmark(
my_corrected_ppi_list <- lapply(files, process_radar_file)
)

print(benchmark_results)

clean_list <- Filter(Negate(is.null), my_corrected_ppi_list)
#saveRDS(clean_list, file = paste0("composite_ppi/rds/",radar, "_week", week, "_", hours_after, "_V3_25_corrected_ppi_list.Rds"))

my_composite_ppi_mean = composite_ppi(clean_list, nx = 300, ny = 300, method="mean")

saveRDS(my_composite_ppi_mean, file = paste0('output/', radar, "_week", week, "_", hours_after, "_V3_25_mean_composite_ppi.Rds"))

#png(filename = paste0("composite_ppi/plot/", radar, "_week", week, "_", hours_after, "_V3_25_mean_composite_ppi_plot.png"), width = 1600, height = 1200, res = 300)
#plot(my_composite_ppi_mean, 'VID') + labs(title=radar, subtitle=paste('week', week, 'mean', sep = ' '))
#dev.off()

stopCluster(cl)
rm(cl)
gc()


