library(terra)
library(bioRad)
library(logger)

#setwd("~/clo/bigbird/wgfd")

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

clean_weather <- function(my_pvol) {
  
  #filter out non-mistnet scans
  my_pvol$scans <- Filter(function(scan) "WEATHER" %in% names(scan$params), my_pvol$scans)
  for (i in seq_along(my_pvol$scans)) {
    
    # for each scan, extract WEATHER and DBZH parameters
    weather <- my_pvol$scans[[i]]$params$WEATHER
    dbzh <- my_pvol$scans[[i]]$params$DBZH
    param_attrs <- attributes(dbzh)
    
    class(weather) <- "matrix"
    class(dbzh) <- "matrix"
    
    # Set DBZH values to NA where WEATHER > 0.45
    dbzh[weather > 0.45] <- NA
    
    # Restore original class of dbzh and update with the masked scan
    class(dbzh) <- c("param", "matrix", "array")
    my_pvol$scans[[i]]$params$DBZH <- dbzh
    attributes(my_pvol$scans[[i]]$params$DBZH) <-param_attrs
  }
  return(my_pvol)
}


clean_pvol <- function(file, occult, method, threshold) {
  
  if (!method %in% c('subtract', 'occultation')) {
    stop("Invalid method. Please choose either 'subtract' or 'occultation'.")
  }
  
  apply_mistnet(file) -> my_pvol
  my_pvol = clean_weather(my_pvol)
  
  # proceed based on the method
  if (method == 'subtract') {
    
    log_info("subtracting latent reflectivity")
    
    for (lyr_idx in 1:terra::nlyr(subtract_mask)){
      
      pvol_data <- my_pvol$scans[[lyr_idx]]$params$DBZH
      param_attrs <- attributes(pvol_data)
      class(pvol_data) <- "matrix"
      pvol_raster <- rast(pvol_data)
      
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
        
        mask <- occult_resampled > threshold
        pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
        values <- values(pvol_raster_masked)
        pvol_raster_masked_matrix <- matrix(values, nrow=nrow(pvol_raster_masked), byrow=TRUE)
        class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
        my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix
        attributes(my_pvol$scans[[pvol_idx]]$params$DBZH) <- param_attrs
        
      } else {
        
        log_info('applying mask above {threshold} without resampling')
        
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


#Example

my_occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
my_occult_scan <- my_occult$scans[[6]] #0.5 deg
plot(project_as_ppi(my_occult_scan, range_max = 150000), zlim=c(0,1))

pvol_files = list.files('pvol/KRIW/', full.names = TRUE)
my_pvol_file = pvol_files[length(pvol_files)-1]

#clean all scans in a pvol with > 0.2 occultation
my_clean_pvol = clean_pvol(my_pvol_file, my_occult, method='occultation', threshold = 0.2)


