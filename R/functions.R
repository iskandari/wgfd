library(terra)
library(raster)
library(bioRad)
library(raster)
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

eta_to_dbz(dbz_to_eta(dbzh, wavelength=10.71) + 
              + dbz_to_eta(dbzh, wavelength=10.71) * occult)


            
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
    param_attrs <- attributes(dbzh)
    
    class(cell) <- "matrix"
    class(dbzh) <- "matrix"
    
    # Set DBZH values to NA where CELL > 1
    dbzh[cell > 1] <- NA
    
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

#Get height correction coefficient
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


subtract_mask= terra::rast('occult/subtraction/KRIW_top_quartile_raster.tif')
my_occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
my_occult_scan <- my_occult$scans[[6]] #0.5 deg
plot(project_as_ppi(my_occult_scan, range_max = 150000), zlim=c(0,1))


rand_idx = sample(70:100, 1)

pvol_files = list.files('pvol/KRIW/', full.names = TRUE)
my_pvol_file = pvol_files[length(pvol_files)-rand_idx]
my_vp = read_vpfiles(paste0('vp/KRIW/', paste0(basename(pvol_files[length(pvol_files)-rand_idx]),'.h5')))

#clean all scans in a pvol with occultation
my_clean_pvol = clean_pvol(my_pvol_file, my_occult,threshold=0.5, method='occultation', amplify=TRUE)
my_ppi = integrate_to_ppi(my_clean_pvol, my_vp, res = 250, xlim = c(-150000, 150000), ylim = c(-150000, 150000))
plot(my_ppi, param = "VIR")

hc_coeff_raster <- process_dem_correction(my_vp, my_ppi, dem_raster)
hc_coeff <- resample_to_ppi(hc_coeff_raster, my_ppi)

my_corrected_ppi = my_ppi
my_corrected_ppi$data$VIR <-  my_ppi$data$VIR * as.vector(values(hc_coeff))
plot(my_corrected_ppi, param = "VIR")


par(mfrow=c(1, 2))
plot(my_ppi, param = "VIR", main = "Original PPI")
plot(my_corrected_ppi, param = "VIR", main = "Corrected PPI")



