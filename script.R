library(bioRad)
library(terra)
library(stringr)
library(vol2birdR)
library(logger)
setwd("~/clo/bigbird/wgfd")

---

#helper functions

clean_pvol <- function(my_pvol) {

  #filter out non-mistnet scans
  my_pvol$scans <- Filter(function(scan) "BIOLOGY" %in% names(scan$params), my_pvol$scans)
  for (i in seq_along(my_pvol$scans)) {
    # For each scan, extract CELL and DBZH parameters
    cell <- my_pvol$scans[[i]]$params$CELL
    dbzh <- my_pvol$scans[[i]]$params$DBZH

    class(cell) <- "matrix"
    class(dbzh) <- "matrix"

    # Set DBZH values to NA where CELL >= 1
    dbzh[cell >= 1] <- NA

    # Restore original class of dbzh and update the pvol object
    class(dbzh) <- c("param", "matrix", "array")
    my_pvol$scans[[i]]$params$DBZH <- dbzh
  }
  return(my_pvol)
}

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

---
#test single week

radar = "KRIW"
year = 2019
wn = 23
threshold <- 0.1 #occultation threshold

closest_times <- read.csv(paste0("pvol/", radar, "_closest_times.csv"))
df_filtered <- subset(closest_times, closest_file != "")

unique_year_week_df <- unique(data.frame(year = closest_times$year, week = closest_times$week))
unique_year_week_list <- split(unique_year_week_df, seq(nrow(unique_year_week_df)))
filtered_list <- Filter(function(x) x$year == year, unique_year_week_list)

df_subset = subset(closest_times, year == filtered_list[[wn]]$year & week == filtered_list[[wn]]$week & closest_file != '')
df_subset$closest_file <- gsub("_MDM", "", df_subset$closest_file)

----
#time full script

occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
ppi_list = list()

generate_ppis <- function(closest_file) {

  input_pvol_path <- paste0("pvol/", radar, "/", closest_file)
  input_vp_path <- paste0("vp/", radar, "/", tools::file_path_sans_ext(closest_file), ".h5")

  my_vp = read_vpfiles(input_vp_path)

  apply_mistnet(input_pvol_path) -> my_pvol
  my_pvol = clean_pvol(my_pvol)

  #match pvol scan to occult by elevation angle
  occult_index_mapping <- find_proximal_occult_scan_idx(occult, my_pvol)

  for (pvol_idx in seq_along(occult_index_mapping)){

    occult_idx <- occult_index_mapping[[pvol_idx]]

    # compare dimensions of xtract the pvol and occult data
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
      pvol_raster_masked_matrix <- matrix(values(pvol_raster_masked), nrow = nrow(pvol_raster_masked), ncol = ncol(pvol_raster_masked))
      class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
      my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix

    } else {

      log_info('applying mask without resampling')

      mask <- occult_raster > threshold
      pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
      values <- values(pvol_raster_masked)
      pvol_raster_masked_matrix <- matrix(values, nrow=nrow(pvol_raster_masked), byrow=TRUE)
      class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
      my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix
    }
  }

  my_ppi = integrate_to_ppi(my_pvol, my_vp, res = 250, xlim = c(-100000, 100000), ylim = c(-100000, 100000))
  output_png_file <- paste0("plots/vid/",radar,"/corrected/",
                            radar, "_", closest_file,
                            ".png")
  png(filename=output_png_file, width=800, height=600)
  print(plot(my_ppi, param = "VID"))
  dev.off()

  #ppi_list <<- c(ppi_list, list(my_ppi))
}


lapply(df_subset$closest_file, generate_ppis)
my_composite = composite_ppi(ppi_list, param = "VID", method = "mean", res=250)


output_png_file <- paste0("plots/vid/",radar,"/corrected/",
                          radar, "_",
                          paste(filtered_list[[wn]]$year,
                                filtered_list[[wn]]$week, sep="_"),
                          ".png")

# Open a PNG device
png(filename=output_png_file, width=800, height=600)
print(plot(my_composite))
dev.off()

r <- rast(my_composite$data)

output_tif_file <- paste0("tif/", radar, "/",
                          radar, "_",
                          paste(filtered_list[[wn]]$year,
                                filtered_list[[wn]]$week, sep="_"),
                          ".tif")

terra::writeRaster(r, output_tif_file, overwrite=TRUE)
end_time <- Sys.time()

---
# process vp
# use multi-year average for the same week number across years
# 2013 - fall first week of september, spring pick a week around May 10

# WEEKLY AGGREGATES
# for the full aggregate use the peak week numbers
# december during the day - look at the vid


year = 2019
radar <- "KCYS"
file_dir <- paste("/Users/at744/clo/bigbird/wgfd/pvol/", radar, sep="")
output_dir <- paste("/Users/at744/clo/bigbird/wgfd/vp/", radar, sep="")
files <- list.files(file_dir, pattern = paste0("KCYS", year), full.names = TRUE)

for (file in files) {
  message("Processing file: ", file)

  output_file <- paste0(output_dir, "/", gsub(".*(KCYS2019.+)$", "\\1.h5",  basename(tools::file_path_sans_ext(file))))

  config <- vol2bird_config()
  config$clutterMap <- paste('/Users/at744/clo/bitbucket/birdcast/is-birdcast-observed-profile/occult/', radar, ".h5", sep = "")
  config$useClutterMap <- TRUE
  config$clutterValueMin <- 0.05
  config$useMistNet <- TRUE
  config$nLayers <- 50
  config$layerThickness <- 100
  config$maxNyquistDealias <- 50
  config$stdDevMinBird <- 1

  vol2bird(file = file, output_file, config = config)

  message("Output saved to: ", output_file)
}

---

pvol_dir <- 'pvol/biofree/KRIW/raw'
pvol_files <- list.files(pvol_dir, full.names = TRUE)

process_pvol <- function(file_path) {
  tryCatch({
    apply_mistnet(file_path) %>%
      clean_pvol()
  }, error = function(e) {
    message("Error processing file: ", file_path, "\nError message: ", e$message)
    return(NULL)  # Return NULL for this file to indicate the error
  })
}

cleaned_pvols <- lapply(pvol_files, process_pvol)
cleaned_pvols <- cleaned_pvols[!sapply(cleaned_pvols, is.null)]

mat_list <- list()

for (i in seq_along(cleaned_pvols)) {
  dbzh = cleaned_pvols[[i]]$scans[[1]]$params$DBZH
  class(dbzh) <- "matrix"
  mat_list[[i]] <- dbzh
}

#all_matrices <- array(unlist(mat_list), dim = c(dim(mat_list[[1]]), length(mat_list)))
#median_matrix <- apply(all_matrices, c(1, 2), median)
#median_matrix[median_matrix == 0] <- NA
#class(median_matrix) <- c("param", "matrix", "array")

min_matrix <- Reduce(function(x, y) pmin(x, y, na.rm = TRUE), mat_list)
class(min_matrix) <- c("param", "matrix", "array")
my_pvol$scans[[1]]$params$DBZH <- min_matrix

plot(my_pvol$scans[[1]])
my_pvol = cleaned_pvols[1]
my_pvol$scans[[1]]$params$DBZH <- min_matrix

cleaned_pvols[[1]]$scans[[1]]$params$DBZH

----------

#match pvol scan to occult by elevation angle
occult_index_mapping <- find_proximal_occult_scan_idx(occult, my_pvol)

for (pvol_idx in seq_along(occult_index_mapping)){

  occult_idx <- occult_index_mapping[[pvol_idx]]

  # compare dimensions of xtract the pvol and occult data
  pvol_scan_dim <- dim(my_pvol$scans[[pvol_idx]]$params$DBZH)
  occult_scan_dim <- dim(occult$scans[[occult_idx]]$params$OCCULT)

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
    pvol_raster_masked_matrix <- as.matrix(pvol_raster_masked, nrows=pvol_scan_dim[1], ncols=pvol_scan_dim[2])
    class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
    my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix

  } else {

    log_info('applying mask without resampling')

    mask <- occult_raster > threshold
    pvol_raster_masked <- terra:::ifel(mask, NA, pvol_raster)
    pvol_raster_masked_matrix <- as.matrix(pvol_raster_masked, nrows=pvol_scan_dim[1], ncols=pvol_scan_dim[2])
    class(pvol_raster_masked_matrix) <- c("param", "matrix", "array")
    my_pvol$scans[[pvol_idx]]$params$DBZH <- pvol_raster_masked_matrix

  }
}

