library(bioRad)
library(terra)
library(stringr)
library(vol2birdR)
library(logger)
library(parallel)
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
    param_attrs <- attributes(dbzh)

    class(cell) <- "matrix"
    class(dbzh) <- "matrix"

    # Set DBZH values to NA where CELL >= 1
    dbzh[cell >= 1] <- NA

    # Restore original class of dbzh and update the pvol object
    class(dbzh) <- c("param", "matrix", "array")
    my_pvol$scans[[i]]$params$DBZH <- dbzh
    attributes(my_pvol$scans[[i]]$params$DBZH) <-param_attrs
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
subtract_mask = terra::rast(paste0("occult/subtraction/",radar,"_top_quartile_raster.tif"))
filter_year = 2013
wn = 19

#threshold <- 0.1 #occultation threshold

closest_times <- read.csv(paste0("pvol/", radar, "_closest_times.csv"))
df_filtered <- subset(closest_times, closest_file != "")
df_subset <- subset(df_filtered, week == wn & year >= filter_year)
df_subset$closest_file <- gsub("_MDM", "", df_subset$closest_file)


#unique_year_week_df <- unique(data.frame(year >= closest_times$year, week = closest_times$week))
#unique_year_week_list <- split(unique_year_week_df, seq(nrow(unique_year_week_df)))
#filtered_list <- Filter(function(x) x$year == year, unique_year_week_list)
#df_subset = subset(closest_times, year == filtered_list[[wn]]$year & week == filtered_list[[wn]]$week & closest_file != '')
#df_subset$closest_file <- gsub("_MDM", "", df_subset$closest_file)

----
#time full script

occult= read_pvolfile('occult/h5-150km/KRIW/KRIW.h5', param='all')
ppi_list = list()

generate_ppis <- function(closest_file, method) {

  if (!method %in% c("subtract", "occult")) {
    stop("Invalid method. Please choose either 'subtract' or 'occut'.")
  }

  input_pvol_path <- paste0("pvol/", radar, "/", closest_file)
  input_vp_path <- paste0("vp/", radar, "/", tools::file_path_sans_ext(closest_file), ".h5")

  if (!file.exists(input_vp_path)) {
    log_info(paste("VP file does not exist:", input_vp_path, "- Skipping file."))
    return(NULL)
  }

  my_vp = read_vpfiles(input_vp_path)

  apply_mistnet(input_pvol_path) -> my_pvol
  my_pvol = clean_pvol(my_pvol)

  # Proceed based on the method
  if (method == "subtract") {

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

  } else if (method == "occult") {

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
        values <- values(pvol_raster_masked)
        pvol_raster_masked_matrix <- matrix(values, nrow=nrow(pvol_raster_masked), byrow=TRUE)
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
  }

  log_info("integrating to ppi")

  my_ppi = integrate_to_ppi(my_pvol, my_vp, res = 250, xlim = c(-100000, 100000), ylim = c(-100000, 100000))

  #output_png_file <- paste0("plots/vid/",radar,"/corrected/",
  #                          radar, "_", closest_file,
  #                          ".png")
  #png(filename=output_png_file, width=800, height=600)
  #print(plot(my_ppi, param = "VID"))
  #dev.off()

  ppi_list <<- c(ppi_list, list(my_ppi))
}

lapply(df_subset[-(1:22),]$closest_file, function(x) generate_ppis(x, method='subtract'))
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


radar <- "KCYS"
year_start = 2013
year_end = 2013
year_pattern <- paste0(radar, "(", paste(year_start:year_end, collapse="|"), ")")
file_dir <- paste("/Users/at744/clo/bigbird/wgfd/pvol/", radar, sep="")
output_dir <- paste("/Users/at744/clo/bigbird/wgfd/vp/", radar, sep="")
files <- list.files(file_dir, pattern = year_pattern, full.names = TRUE)


for (file in files) {
  message("Processing file: ", file)
  pattern <- paste0(".*(", radar, year, ".+)$")
  output_file <- paste0(output_dir, "/", gsub(pattern, "\\1.h5", basename(tools::file_path_sans_ext(file))))

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

process_file <- function(file, radar, output_dir) {

  filename <- basename(tools::file_path_sans_ext(file))
  message("Processing file: ", file)
  output_file <- paste0(output_dir, "/", filename, ".h5")

  config <- vol2bird_config()
  config$clutterMap <- paste0('/Users/at744/clo/bitbucket/birdcast/is-birdcast-observed-profile/occult/', radar, ".h5")
  config$useClutterMap <- TRUE
  config$clutterValueMin <- 0.05
  config$useMistNet <- TRUE
  config$nLayers <- 50
  config$layerThickness <- 100
  config$maxNyquistDealias <- 50
  config$stdDevMinBird <- 1

  # Execute vol2bird
  vol2bird(file = file, output_file, config = config)
  message("Output saved to: ", output_file)
}

numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("vol2bird", "vol2bird_config", "process_file", "output_dir", "radar"))
clusterEvalQ(cl, {
  library(dplyr)
  library(vol2birdR)
})

parLapply(cl, sample(files), function(file) process_file(file, radar, output_dir))

---

#construct subtraction rasters

cleaned_pvols <- lapply(pvol_files, process_pvol)
cleaned_pvols <- cleaned_pvols[!sapply(cleaned_pvols, is.null)]

elevation_angles_vector <- sapply(cleaned_pvols, function(pvol) {
  pvol$scans[[1]]$attributes$where$elangle
})


pvol_dir <- 'pvol/biofree/KRIW'
pvol_files <- list.files(pvol_dir, full.names = TRUE)[1:100]

process_pvol <- function(file_path) {
  tryCatch({
    apply_mistnet(file_path) %>%
      clean_pvol()
  }, error = function(e) {
    message("Error processing file: ", file_path, "\nError message: ", e$message)
    return(NULL)  # Return NULL for this file to indicate the error
  })
}

numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("process_pvol", "apply_mistnet", "clean_pvol"))
clusterEvalQ(cl, {
  library(dplyr)
  library(bioRad)
})


raster_list = list()

pvol_dir <- 'pvol/biofree/KRIW'
all_pvol_files <- setdiff(list.files(pvol_dir, full.names = TRUE), list.dirs(pvol_dir, recursive = FALSE, full.names = TRUE))
num_files <- length(all_pvol_files)
target_param = 'VRADH'

# Loop through files in chunks of 100, starting from 101 to the end
for (start_idx in seq(0, num_files, by = 100)) {
  end_idx <- min(start_idx + 99, num_files) # Ensure we don't go beyond the array length
  pvol_files <- all_pvol_files[start_idx:end_idx]

  # Process each chunk of files
  cleaned_pvols <- parLapply(cl, pvol_files, process_pvol)

  for (i in seq_along(cleaned_pvols)) {
    param <- cleaned_pvols[[i]]$scans[[1]]$params[[target_param]]
    if (!is.null(param) && all(dim(param) == c(1201, 720))) {
      class(param) <- "matrix"
      raster_list[[length(raster_list) + 1]] <- rast(param) # Append new raster to the list
    }
  }

  # Optionally clear cleaned_pvols to free memory
  rm(cleaned_pvols)
  gc() # Garbage collection
}


cleaned_pvols <- parLapply(cl, pvol_files, process_pvol)
raster_list <- list()

for (i in seq_along(cleaned_pvols)) {
  dbzh = cleaned_pvols[[i]]$scans[[1]]$params$DBZH
    if (!is.null(dbzh) && all(dim(dbzh) == c(1201, 720))) {
      class(dbzh) <- "matrix"
      raster_list[[i]] <- rast(dbzh)
    }
}

rm(cleaned_pvols)
gc()


stacked_raster <- rast(raster_list)
median_raster <- app(stacked_raster, fun = median, na.rm = TRUE)
count_non_na <- function(x) {
  sum(!is.na(x))
}

count_raster <- app(stacked_raster, fun = count_non_na)

count_values <- values(count_raster)
median_values <- values(median_raster)
threshold = 250
valid_indices <- which(count_values >= threshold & !is.na(count_values))
valid_median_values <- median_values[valid_indices]


#plot distribution of values form randomly selected cells
valid_cells <- which(values(count_raster) >= threshold)
selected_cell <- sample(valid_cells, 1)
coords <- xyFromCell(count_raster, cell=selected_cell)
points <- vect(coords, crs=crs(count_raster))
selected_cell_values <- terra::extract(stacked_raster, points)

long_format <- pivot_longer(selected_cell_values, cols = starts_with("lyr"),
                            names_to = "Layer", values_to = "Value")
values_vector <- long_format$Value
data_for_plot <- data.frame(Values = values_vector)


plot_file_path = "pvol/biofree/KRIW/cell_sample/"
plot_file_name <- paste("violin_plot_x", coords[1], "y", coords[2], ".png", sep="")
plot_file_path_full <- paste0(plot_file_path, plot_file_name)
png(filename = plot_file_path_full)
vioplot(values_vector, main=paste("x:", coords[1], "y:", coords[2],
                                  "\n", "n=", length(na.omit(values_vector)), sep=" "), ylab='DBZH')
dev.off()



adjusted_raster <- ifel(median_raster > 10, median_raster, NA)


count_values <- as.vector(values(count_raster))
median_values <- as.vector(values(median_raster))

data_for_clustering <- data.frame(Count = count_values, Median = median_values)
data_for_clustering <- na.omit(data_for_clustering)
data_for_clustering <- data_for_clustering[!is.infinite(rowSums(data_for_clustering)), ]


wss <- numeric(10)
for (k in 1:10) {
  set.seed(123)
  clusters <- kmeans(data_for_clustering, centers = k, nstart = 20)
  wss[k] <- clusters$tot.withinss
}


plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing k")



set.seed(123)
clusters <- kmeans(data_for_clustering, centers = 2)

cluster_assignment <- clusters$cluster

# Optionally, plot the results
plot(data_for_clustering$Count, data_for_clustering$Median, col = cluster_assignment)
legend("topright", legend = unique(cluster_assignment), col = 1:3, pch = 1)




top_quartile_threshold <- quantile(valid_median_values, 0.75)
top_quartile_raster <- ifel(median_raster > top_quartile_threshold, median_raster, NA)
values = values(top_quartile_raster)

values = values(median_raster[[1]])
values_with_na = values[is.nan(values)] <- NA
raster_matrix <- matrix(values_with_na, nrow=nrow(median_raster[[1]]), byrow=TRUE)

class(raster_matrix) <- c("param", "matrix", "array")
my_pvol$scans[[1]]$params$DBZH  <- raster_matrix
attributes(my_pvol$scans[[1]]$params$DBZH) <- attributes(my_proxy_pvol$scans[[1]]$params$DBZH)



top_quartile_raster_masked_matrix <- matrix(values, nrow=nrow(top_quartile_raster), byrow=TRUE)
top_quartile_raster_masked_matrix[is.nan(top_quartile_raster_masked_matrix)] <- NA
class(top_quartile_raster_masked_matrix) <- c("param", "matrix", "array")


my_proxy_pvol = cleaned_pvols[[100]]
my_pvol = my_proxy_pvol
my_pvol$scans[[1]]$params$DBZH <- top_quartile_raster_masked_matrix
attributes(my_pvol$scans[[1]]$params$DBZH) <- attributes(my_proxy_pvol$scans[[1]]$params$DBZH)
plot(project_as_ppi(my_pvol$scans[[1]], range_max=100000), main='Median DBZH')

raster_output = list()
raster_output[[2]] <- top_quartile_raster
stacked_raster <- rast(raster_output)
terra::writeRaster(stacked_raster, filename = "occult/subtraction/KRIW_top_quartile_raster.tif", overwrite = TRUE)


output_tif_file <- paste0("tif/", radar, "/",
                          radar, "_",
                          paste(filtered_list[[wn]]$year,
                                filtered_list[[wn]]$week, sep="_"),
                          ".tif")

terra::writeRaster(r, output_tif_file, overwrite=TRUE)


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








