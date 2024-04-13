library(bioRad)
library(terra)
library(stringr)
library(vol2birdR)
library(logger)
library(parallel)
library(tidyr)
library(vioplot)
setwd("~/clo/bigbird/wgfd")

---

#helper functions

clean_pvol <- function(my_pvol) {

  #filter out non-mistnet scans
  my_pvol$scans <- Filter(function(scan) "WEATHER" %in% names(scan$params), my_pvol$scans)
  for (i in seq_along(my_pvol$scans)) {
    # For each scan, extract CELL and DBZH parameters
    weather <- my_pvol$scans[[i]]$params$WEATHER
    dbzh <- my_pvol$scans[[i]]$params$DBZH
    param_attrs <- attributes(dbzh)

    class(weather) <- "matrix"
    class(dbzh) <- "matrix"

    # Set DBZH values to NA where CELL >= 1
    dbzh[weather > 0.45] <- NA

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
file_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/raw"
output_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/mistnet_clean"
files <- list.files(file_dir, full.names = TRUE)


for (file in files) {
  message("Processing file: ", file)
  
  filename = basename(tools::file_path_sans_ext(file))
  output_file <- paste0(output_dir, filename, '.csv')
  radar = substr(filename, 1, 4)
  
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
  radar = substr(filename, 1, 4)
  output_file <- paste0(output_dir, "/", filename, ".csv")

  config <- vol2bird_config()
  config$clutterMap <- paste0('/Users/at744/clo/bitbucket/birdcast/is-birdcast-observed-profile/occult/', radar, ".h5")
  config$useClutterMap <- TRUE
  config$clutterValueMin <- 0.05
  config$useMistNet <- TRUE
  config$nLayers <- 50
  config$layerThickness <- 100
  config$maxNyquistDealias <- 50
  config$stdDevMinBird <- 1

  # Execute vol2bird inside tryCatch to handle potential errors
  tryCatch({
    vol2bird(file = file, output_file, config = config)
    message("Output saved to: ", output_file)
  }, error = function(e) {
    message("Error processing file ", file, ": ", e$message)
  })
}

numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("vol2bird", "vol2bird_config", "process_file", "output_dir", "radar"))
clusterEvalQ(cl, {
  library(dplyr)
  library(vol2birdR)
})

parLapply(cl, sample(files), function(file) process_file(file, radar, output_dir))

file_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/additional4"

output_files <- basename(tools::file_path_sans_ext(list.files(output_dir, full.names = TRUE)))
input_files <- basename(tools::file_path_sans_ext(list.files(file_dir, full.names = TRUE)))
files <- list.files(file_dir, full.names = TRUE)


#Undefined error: 0
"/Users/at744/clo/bigbird/wgfd/pvol/eclipse/KRGX20240407_183101_V06"
"/Users/at744/clo/bigbird/wgfd/pvol/eclipse/KHDX20240407_191228_V06"
"/Users/at744/clo/bigbird/wgfd/pvol/eclipse/KPUX20240406_184953_V06"

full_paths <- paste0("/Users/at744/clo/bigbird/wgfd/pvol/eclipse/", missing_files)
---
  
radar_obscuration = read.csv("new_radar_obscuration2.csv")
radar_obscuration = radar_obscuration[radar_obscuration$closest_file != '',]
radar_obscuration$closest_file_path =  paste0("/Users/at744/clo/bigbird/wgfd/vp/eclipse/", radar_obscuration$closest_file, '.csv')

file_paths <- radar_obscuration$closest_file_path
files_exist <- file.exists(file_paths)
missing_files <- file_paths[!files_exist]

# Function to transform the file path into the desired S3 path format
transform_to_s3_path <- function(file_path) {
  # Extract the file name without the directory path and extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Extract the date and radar from the file name
  radar <- substr(file_name, 1, 4)
  year <- substr(file_name, 5, 8)
  month <- substr(file_name, 9, 10)
  day <- substr(file_name, 11, 12)
  
  # Construct the S3 path
  s3_path <- paste(year, month, day, radar, file_name, sep="/")
  
  return(s3_path)
}

# Apply the transformation to each file path
s3_paths <- sapply(missing_files, transform_to_s3_path)

# Create a dataframe with the resulting S3 paths
s3_paths_df <- data.frame(s3_path = s3_paths)

# Print the dataframe
print(s3_paths_df)

radar_obscuration = read.csv("new_radar_obscuration.csv")
radar_obscuration = radar_obscuration[radar_obscuration$closest_file != '',]
radar_obscuration$closest_file_path =  paste0("/Users/at744/clo/bigbird/wgfd/vp/eclipse/", radar_obscuration$closest_file, '.csv')
radar_obscuration$timestamp <- as.POSIXct(radar_obscuration$timestamp, tz = "UTC")
radar_obscuration$date <- as.Date(radar_obscuration$timestamp, tz = "UTC")

my_vpts = read_vpts("/Users/at744/clo/bigbird/wgfd/vp/eclipse/KABR20240407_165720_V06.csv")

required_features = colnames(example_vp$data)
vpts_to_vp <- function(my_vpts, required_features) {
  
  # Check if my_vpts is of class "vpts"
  if (!is.vpts(my_vpts)) {
    stop("The provided object is not of class 'vpts'.")
  }
  
  # Check required features
  required_features <- colnames(example_vp$data)
  if (length(setdiff(required_features, colnames(as.data.frame(my_vpts)))) != 0) {
    stop("The provided my_vpts doesn't have all required features.")
  }
  
  # Ensure non-null components
  if (is.null(my_vpts$radar) || is.null(my_vpts$datetime) || is.null(my_vpts$attributes)) {
    stop("radar, datetime, or attributes in my_vpts are NULL.")
  }
  
  # Construct the object
  my_vp <- list(
    data = as.data.frame(my_vpts),
    radar = my_vpts$radar,
    datetime = my_vpts$datetime,
    attributes = my_vpts$attributes
  )
  
  # Assign class vp
  class(my_vp) <- "vp"
  
  return(my_vp)
}

  

my_vpts = read_vpts(closest_file_path)
my_vp = vpts_to_vp(my_vpts)
my_integrated_vp = integrate_profile(my_vp)
my_integrated_vp$vir

# Define the function to process each file
process_file <- function(row) {
  if (file.exists(row$closest_file_path)) {
    my_vpts <- read_vpts(row$closest_file_path)
    my_vp <- vpts_to_vp(my_vpts)
    my_integrated_vp <- integrate_profile(my_vp)
    
    # Return a list or data frame containing both vir and minutes_from_eclipse
    return(data.frame(vir = my_integrated_vp$vir, minutes_from_eclipse = row$minutes_from_eclipse))
  }
  # Return NULL or an empty data frame if the condition is not met
  return(data.frame(vir = NA, minutes_from_eclipse = row$minutes_from_eclipse))
}

numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, c("read_vpts", "vpts_to_vp", "integrate_profile"))
clusterEvalQ(cl, {
  library(dplyr)
  library(bioRad)
})


radar_obscuration = read.csv("eclipse_day.csv")

radar_obscuration = radar_obscuration[radar_obscuration$closest_file != '',]
radar_obscuration = radar_obscuration[radar_obscuration$radar == 'KPAH',]
radar_obscuration$closest_file_path =  paste0("/Users/at744/clo/bigbird/wgfd/vp/eclipse/", radar_obscuration$closest_file, '.csv')
radar_obscuration$timestamp <- as.POSIXct(radar_obscuration$timestamp, tz = "UTC")
radar_obscuration$date <- as.Date(radar_obscuration$timestamp, tz = "UTC")

radar_obscuration_subset = radar_obscuration
#radar_obscuration_subset = radar_obscuration[radar_obscuration$max_percent < 80,]
#radar_obscuration_subset = radar_obscuration[radar_obscuration$max_percent >= 80 & radar_obscuration$max_percent < 95,]
#radar_obscuration_subset = radar_obscuration[radar_obscuration$max_percent >= 95,]
#radar_obscuration_subset = radar_obscuration[radar_obscuration$max_percent == 100,]

radar_obscuration_list <- split(radar_obscuration_subset, seq(nrow(radar_obscuration_subset)))
results <- parLapply(cl, radar_obscuration_list, process_file)
vir_data_frame <- as.data.frame(do.call(rbind, results))
vir_data_frame$minutes_from_eclipse <- as.numeric(as.integer(as.character(vir_data_frame$minutes_from_eclipse)))

#vir_data_frame_80 <- vir_data_frame
#vir_data_frame_80_95 <- vir_data_frame
#vir_data_frame_95 <- vir_data_frame
#vir_data_frame_100 <- vir_data_frame

vir_data_frame_80$vir <- as.numeric(unlist(vir_data_frame_80$vir))
vir_data_frame_80_95$vir <- as.numeric(unlist(vir_data_frame_80_95$vir))
vir_data_frame_95$vir <- as.numeric(unlist(vir_data_frame_95$vir))
vir_data_frame_100$vir <- as.numeric(unlist(vir_data_frame_100$vir))


library(mgcv)


colors <- c("<80%" = "#E5E455", "80−95%" = "#848151", "95−100%" = "#4D4D4D")

vir_data_frame_2017 
vir_data_frame_2024 

p <- ggplot() +
  #geom_smooth(data = vir_data_frame_80, aes(x = minutes_from_eclipse, y = vir, color = "<80%"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  #geom_smooth(data = vir_data_frame_80_95, aes(x = minutes_from_eclipse, y = vir, color = "80−95%"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  #geom_smooth(data = vir_data_frame_95, aes(x = minutes_from_eclipse, y = vir, color = "95−100%"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  geom_smooth(data = vir_data_frame, aes(x = minutes_from_eclipse, y = vir, color = "black"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 0.25) +
  scale_color_manual(values = colors, name = "obscuration") +
  labs(
    x = "minutes from sunset",
    y = expression(VIR ~ cm^2 ~ km^-2),
    title = "Aug 21 2017 KPAH"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5),  # Half as thick axis lines
    axis.ticks = element_line(color = "black", size = 0.5),  # Half as thick tick marks
    axis.ticks.length = unit(0.125, "cm")  # Adjust tick mark length if necessary
  ) + coord_cartesian(xlim=c(-120, 120)) + 
  coord_cartesian(ylim = c(0, 600))


ggsave("high_res_plot_aug_21_KPAH.png", plot = p, width = 10, height = 6, dpi = 300, bg='white')

ggplot(vir_data_frame_80, aes(x = minutes_from_eclipse, y = vir)) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 0.5) + 
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"), color = "black", size = 1) +
  labs(
    x = "Minutes from 'eclipse' maximum",
    y = expression(VIR ~ cm^2 ~ km^-2),
    title = "Apr 8 2024"
  ) +
  theme_minimal() +
  theme(panel.border = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,axis.line = element_line()
  )



base_plot <- ggplot() +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 0.25) +
  scale_color_manual(values = c("black"), name = "obscuration") +
  labs(
    x = "minutes from sunset",
    y = expression(VIR ~ cm^2 ~ km^-2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.125, "cm")
  ) +
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(-120, 120), ylim = c(0, 600))

plot_2017 <- base_plot +
  geom_smooth(data = vir_data_frame_2017, aes(x = minutes_from_eclipse, y = vir, color = "Year 2017"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  labs(title = "Aug 21 2017 KPAH")

# Plot for 2024
plot_2024 <- base_plot +
  geom_smooth(data = vir_data_frame_2024, aes(x = minutes_from_eclipse, y = vir, color = "Year 2024"), method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
  labs(title = "Apr 8 2024 KPAH")

# Arrange plots side by side
grid.arrange(plot_2017, plot_2024, ncol = 2)

arranged_plots <- arrangeGrob(plot_2017, plot_2024, ncol = 2)

# Now use ggsave to save the arranged plot to a PNG file
ggsave("KPAH.png", plot = arranged_plots, width = 16, height = 8, dpi = 300)



#construct subtraction rasters

cleaned_pvols <- lapply(pvol_files, process_pvol)
cleaned_pvols <- cleaned_pvols[!sapply(cleaned_pvols, is.null)]

elevation_angles_vector <- sapply(cleaned_pvols, function(pvol) {
  pvol$scans[[1]]$attributes$where$elangle
})


file_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/raw"
output_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/mistnet_clean/"
files <- list.files(file_dir, full.names = TRUE)


process_pvol <- function(file_path) {
  tryCatch({
    
    #output_dir_raster <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/raster/"
    #output_dir_image <- "/Users/at744/clo/bigbird/wgfd/pvol/biofree/KRIW/image/"
    filename = basename(tools::file_path_sans_ext(file_path))
    my_pvol <- apply_mistnet(file_path) %>% clean_pvol()
    
    #my_param <- get_param(my_pvol$scans[[1]], 'DBZH')
    #my_ppi <- project_as_ppi(my_param, range_max = 100000)
    #my_raster <- rast(my_ppi$data)
    #writeRaster(my_raster, filename = paste0(output_dir_raster, filename,".tif"), overwrite=TRUE)
    
    #raster_df = as.data.frame(my_raster, xy=TRUE) %>% na.omit()
    # p <- (
    #   ggplot(data = raster_df) +
    #     geom_raster(aes(x = x, y = y, fill=DBZH)) +
    #     scale_fill_viridis_c(limit = c(-35, 55)) +
    #     theme_void() +
    #     theme(
    #       panel.background = element_rect(fill = 'black'),
    #       plot.background = element_rect(fill = "black"),
    #       legend.position = "bottom",
    #       legend.text = element_text(color = "#E4EBF7"),  # Set legend text color
    #       legend.title = element_text(color = "#E4EBF7"),  # Set legend title color
    #     )
    # )
    # 
   # ggsave(filename = paste0(output_dir_image, filename, ".png"), plot = p, width = 10, height = 8, dpi = 300)
    
  }, error = function(e) {
    message("Error processing file: ", file_path, "\nError message: ", e$message)
    return(NULL)  # Return NULL for this file to indicate the error
  })
  return(my_pvol)
}


numCores <- 10
cl <- makeCluster(numCores)
clusterExport(cl, varlist = c("process_pvol", "clean_pvol", "apply_mistnet",
                              "get_param", "project_as_ppi", "rast", "writeRaster",
                              "ggplot", "ggsave", "scale_fill_viridis_c", "theme",
                              "element_text","element_rect", "geom_raster",
                              "theme_void", "aes" ,"ggsave"))
clusterEvalQ(cl, {
  library(dplyr)
  library(bioRad)
  library(terra)
  library(magrittr)
  library(ggplot2)
})

file_path <- files[1]

parLapply(cl, files, process_pvol)
lapply(files[1], process_pvol)

raster_df = as.data.frame(my_raster, xy=TRUE) %>% na.omit()

  (
    ggplot(data = raster_df) +
      geom_raster(aes(x = x, y = y, fill=DBZH)) +
      scale_fill_viridis_c() +
      theme_void() +
      theme(
        panel.background = element_rect(fill = 'black'),
        plot.background = element_rect(fill = "black"),
        legend.position = "bottom",
        legend.text = element_text(color = "#E4EBF7"),  # Set legend text color
        legend.title = element_text(color = "#E4EBF7"),  # Set legend title color
      )
  )


raster_list = list()

pvol_dir <- 'pvol/biofree/KRIW/raw'
all_pvol_files <- setdiff(list.files(pvol_dir, full.names = TRUE), list.dirs(pvol_dir, recursive = FALSE, full.names = TRUE))
num_files <- length(all_pvol_files)
target_params = c('DBZH', 'VRADH', 'RHOHV','PHIDP')

# Loop through files in chunks of 100, starting from 101 to the end
for (start_idx in seq(0, num_files, by = 100)) {
  end_idx <- min(start_idx + 99, num_files) # Ensure we don't go beyond the array length
  pvol_files <- all_pvol_files[start_idx:end_idx]
  
  cleaned_pvols <- parLapply(cl, pvol_files, process_pvol)
  
  for (i in seq_along(cleaned_pvols)) {
    layers = list()
    for (target_param in target_params) {
      param <- cleaned_pvols[[i]]$scans[[1]]$params[[target_param]]
      if (!is.null(param) && all(dim(param) == c(1201, 720))) {
        class(param) <- "matrix"
        layers[[length(layers) + 1]] <- rast(param) 
      }
    }
    
    if (length(layers) == length(target_params)) {
      combined_raster <- terra::rast(layers)
      names(combined_raster) <- target_params
      raster_list[[length(raster_list) + 1]] <- combined_raster
    }
  }
  
  rm(cleaned_pvols)
  gc() 
}



#cleaned_pvols <- parLapply(cl, pvol_files, process_pvol)
#raster_list <- list()

for (i in seq_along(cleaned_pvols)) {
  dbzh = cleaned_pvols[[i]]$scans[[1]]$params$DBZH
    if (!is.null(dbzh) && all(dim(dbzh) == c(1201, 720))) {
      class(dbzh) <- "matrix"
      raster_list[[i]] <- rast(dbzh)
    }
}

rm(cleaned_pvols)
gc()


--------

lat = radarInfo[radarInfo$radar == 'KRIW',]$lat
lon = radarInfo[radarInfo$radar == 'KRIW',]$lon
proj4string_az =paste("+proj=aeqd +lat_0=",lat," +lon_0=",lon," +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep="")

dem_raster = brick("gis/elev/elev48i0100a.tif")
lower48_sf_transformed = sf::st_transform(lower48_sf, sf::st_crs(crs(dem_raster)))
wyoming_sf= lower48_sf_transformed[lower48_sf_transformed$name == 'Wyoming',]

# Crop the DEM using the extent of Wyoming
wyoming_extent <- sf::st_bbox(wyoming_sf)
wyoming_dem <- crop(dem_raster, wyoming_extent)

#Trasnform Wyoming DEM
dem_stars <- st_as_stars(wyoming_dem)
dem_transformed <- st_transform(dem_stars, crs = proj4string_az)


ggplot() + geom_stars(data = dem_transformed)





-------

dbzh_layers = list()
vradh_layers = list()# Assuming you have a raster_list with DBZH, VRADH, and RHOHV layers
rhohv_layers = list()

# Extract individual layers and add them to respective lists
for (raster in raster_list) {
  dbzh_layers[[length(dbzh_layers) + 1]] <- raster[["DBZH"]]
  vradh_layers[[length(vradh_layers) + 1]] <- raster[["VRADH"]]
  rhohv_layers[[length(rhohv_layers) + 1]] <- raster[["RHOHV"]]
}

dbzh_stack <- terra::rast(dbzh_layers)
vradh_stack <- terra::rast(vradh_layers)
rhohv_stack <- terra::rast(rhohv_layers)

median_with_threshold <- function(x) {
  if (sum(!is.na(x)) >= 100) {
    return(median(x, na.rm = TRUE))
  } else {
    return(NA)
  }
}

dbzh_stack_agg_mean = aggregate(dbzh_stack, fact = 4, fun = "mean")
#vradh_stack_agg_mean = aggregate(vradh_stack, fact = 4, fun = "mean")
#rhohv_stack_agg_mean = aggregate(rhohv_stack, fact = 3, fun = "mean")

dbzh_median_raster <- app(dbzh_stack_agg_mean, fun = median_with_threshold)
#vradh_median_raster <- app(vradh_stack_agg_mean, fun = median_with_threshold)
#rhohv_median_raster <- app(rhohv_stack_agg_mean, fun = median_with_threshold)

probability_echo_not_na <- function(x) {
  prob_not_na = sum(!is.na(x)) / length(x)
  return(prob_not_na)
}

dbzh_prob_echo_raster <- app(dbzh_stack_agg_mean, fun = probability_echo_not_na)
#vradh_prob_echo_raster <- app(vradh_stack_agg_mean, fun = probability_echo_not_na)
#rhohv_prob_echo_raster <- app(rhohv_stack_agg_mean, fun = probability_echo_not_na)

--------
combined_median_raster <- c(dbzh_median_raster, vradh_median_raster, rhohv_median_raster)
m <- any(is.na(combined_median_raster))
x <- mask(combined_median_raster, m, maskvalue=1)

dbzh_median_raster_masked = mask(dbzh_median_raster, m, maskvalue=1)
vradh_median_raster_masked = mask(vradh_median_raster, m, maskvalue=1)
rhohv_median_raster_masked = mask(rhohv_median_raster, m, maskvalue=1)
--------

combined_median_raster <- c(dbzh_median_raster, vradh_median_raster)
m <- any(is.na(dbzh_median_raster))

dbzh_median_raster_masked = mask(dbzh_median_raster, m, maskvalue=1)
vradh_median_raster_masked = mask(vradh_median_raster, m, maskvalue=1)

dbzh_median_values <- values(dbzh_median_raster_masked)
vradh_median_values <- values(vradh_median_raster_masked)
valid_values_df <- data.frame(dbzh = dbzh_median_values, vradh=vradh_median_values)
valid_values_df <- na.omit(valid_values_df)
colnames(valid_values_df) <- c('dbzh', 'vradh')


m <- any(is.na(dbzh_median_raster))
dbzh_prob_echo_raster_masked <- mask(dbzh_prob_echo_raster, m, maskvalue=1)
dbzh_median_raster_masked = mask(dbzh_median_raster, m, maskvalue=1)
dbzh_median_values <- values(dbzh_median_raster_masked)
dbzh_prob_echo_values <- values(dbzh_prob_echo_raster_masked)

valid_values_df <- data.frame(median = dbzh_median_values, prob_echo=dbzh_prob_echo_values)
valid_values_df <- na.omit(valid_values_df)
colnames(valid_values_df) <- c('median', 'prob')


eps = 0.5
minPts = 5

# Performing DBSCAN
dbscan_result = dbscan(valid_values_df, eps = eps, minPts = minPts)

# Adding cluster assignments to the dataframe
valid_values_df$cluster = dbscan_result$cluster


ggplot(valid_values_df, aes(x = median, y = prob)) +
  geom_point(alpha = 0.05) +  # Use alpha to make points semi-transparent if there's overplotting
  labs(x = "Median DBZH", y = "Probability of Echo", title = "Median DBZH vs. Probability of Echo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    plot.background = element_rect(fill = "white", color = NA),  # White background
    panel.background = element_rect(fill = "white", color = NA)  # White panel background
  )


for (k in 1:10) {
  set.seed(123)  # Set seed for reproducibility
  clusters <- kmeans(valid_values_df, centers = k, nstart = 20)
  wss[k] <- clusters$tot.withinss
}

# Plot the WSS for each k to visualize the elbow method
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing k")

k=3
clusters_kmeans <- kmeans(valid_values_df, centers = k, nstart = 20)
valid_values_df$cluster <- clusters_kmeans$cluster

ggplot(valid_values_df, aes(x = median, y = prob, color = factor(cluster))) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_d(name = "Cluster") + # Corrected for discrete data
  labs(x = "Median DBZH", y = "Probability of Echo", title = "K-means Clustering of Median DBZH vs. Probability of Echo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )

# Plotting
ggplot(valid_values_df, aes(x = median, y = prob, color = factor(cluster))) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_d(name = "Cluster") + # For discrete data
  labs(x = "Median DBZH", y = "Probability of Echo", title = "DBSCAN Clustering of Median DBZH vs. Probability of Echo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  )









raster_df = as.data.frame(rhohv_median_raster, xy=TRUE) %>% na.omit()
p <- ggplot(data = raster_df, aes(x = x, y = y, fill = lyr.1)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_void() + 
  theme(
    axis.text.x = element_text(size = 8), # Smaller font size for x-axis
    axis.text.y = element_text(size = 8), # Smaller font size for y-axis
    panel.background = element_rect(fill = 'black'),
    plot.background = element_rect(fill = "black"),
    legend.text = element_text(color = "#E4EBF7"),  # Set legend text color
    legend.title = element_text(color = "#E4EBF7"),  
    
  ) +
  labs(fill = "RHOHV")

ggsave("rhohv_median_raster.png", plot = p, width = 10, height = 8, dpi = 300)




raster_stack = np.stack([dbzh_median_raster_masked, vradh_median_raster_masked, rhohv_median_raster_masked], axis=0)



# Extract values from the median rasters where the cells are valid
dbzh_values <- values(dbzh_median_raster_masked)
vradh_values <- values(vradh_median_raster_masked)
rhohv_values <- values(rhohv_median_raster_masked)

# Create a dataframe with the valid values from DBZH, VRADH, and RHOHV
valid_values_df <- data.frame(DBZH = dbzh_values, VRADH = vradh_values, RHOHV = rhohv_values)
valid_values_df = na.omit(valid_values_df)
colnames(valid_values_df) <- c('DBZH', 'VRADH', 'RHOHV')

valid_indices <- which(count_values >= threshold & !is.na(count_values))
valid_median_values <- median_values[valid_indices]

#plot distribution of values form randomly selected cells
valid_cells <- which(values(count_raster) >= threshold)
selected_cell <- sample(valid_cells, 1)
coords <- xyFromCell(count_raster, cell=selected_cell)
points <- vect(coords, crs=crs(count_raster))
selected_cell_values <- terra::extract(dbzh_raster, points)



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



png(filename = "wss.png")
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing k")
dev.off()

k <- 3 # Number of clusters chosen based on the elbow method
clusters_kmeans <- kmeans(valid_values_df, centers = 2, nstart = 20)
valid_values_df$cluster <- clusters_kmeans$cluster
pairs(valid_values_df[, c("DBZH", "VRADH", "RHOHV")], col = valid_values_df$cluster, pch = 20,
      main = "Pairwise scatter plots with two clusters")



data_scaled <- scale(valid_values_df)
dist_matrix <- dist(data_scaled, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc, hang = -1) 

wss <- numeric(10)

# Perform k-means clustering for different numbers of clusters (from 1 to 10)
for (k in 1:10) {
  set.seed(123)  # Set seed for reproducibility
  clusters <- kmeans(valid_values_df, centers = k, nstart = 20)
  wss[k] <- clusters$tot.withinss
}

# Plot the WSS for each k to visualize the elbow method
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares",
     main = "Elbow Method for Choosing k")




library(factoextra)



wss <- numeric(10)
silhouette_scores <- numeric(10)

# Perform k-means clustering and calculate silhouette for different numbers of clusters
for (k in 2:10) {  # silhouette analysis requires at least 2 clusters to be meaningful
  set.seed(123)  # Set seed for reproducibility
  clusters <- kmeans(valid_values_df, centers = k, nstart = 20)
  ss <- cluster::silhouette(clusters$cluster, dist(valid_values_df))
  silhouette_scores[k] <- mean(ss[, 3])  # average silhouette width
}

# Determine the optimal number of clusters based on the highest average silhouette score
optimal_clusters <- which.max(silhouette_scores)

optimal_clustering <- kmeans(valid_values_df, centers = optimal_clusters, nstart = 20)
silhouette_plot <- silhouette(optimal_clustering$cluster, dist(valid_values_df))
fviz_silhouette(silhouette_plot)


set.seed(123)  # Ensure reproducibility
optimal_clusters <- which.max(silhouette_scores)  # Assuming you have computed this from silhouette analysis
clusters <- kmeans(valid_values_df, centers = optimal_clusters, nstart = 20)

# Plot the data points, color them by cluster assignment
plot(valid_values_df[,1], valid_values_df[,2], col=clusters$cluster, pch=20, main="K-Means Clustering", xlab="DBZH", ylab="VRADH")
legend("topright", legend=paste("Cluster", unique(clusters$cluster)), fill=unique(as.numeric(clusters$cluster)), title="Clusters")


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








