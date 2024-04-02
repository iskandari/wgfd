library(bioRad)
library(terra)
library(stringr)
library(vol2birdR)
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


---
#test single week


radar = "KRIW"
year = 2019
wn = 23

closest_times <- read.csv(paste0("pvol/", radar, "_closest_times.csv"))
df_filtered <- subset(closest_times, closest_file != "")

unique_year_week_df <- unique(data.frame(year = closest_times$year, week = closest_times$week))
unique_year_week_list <- split(unique_year_week_df, seq(nrow(unique_year_week_df)))
filtered_list <- Filter(function(x) x$year == year, unique_year_week_list)


df_subset = subset(closest_times, year == filtered_list[[wn]]$year & week == filtered_list[[wn]]$week & closest_file != '')
df_subset$closest_file <- gsub("_MDM", "", df_subset$closest_file)


----
#time full script

start_time <- Sys.time()
occult= read_pvolfile('occult/h5-150km/KRIW/12/KRIW.h5', param='all')

ppi_list = list()
generate_ppis <- function(closest_file) {

  input_pvol_path <- paste0("pvol/", radar, "/", closest_file)
  input_vp_path <- paste0("vp/", radar, "/", tools::file_path_sans_ext(closest_file), ".h5")

  my_vp = read_vpfiles(input_vp_path)

  apply_mistnet(input_pvol_path) -> my_pvol
  my_pvol = clean_pvol(my_pvol)

  blockage_threshold <- 0.1
  num_scans_to_process <- min(length(occult$scans), length(my_pvol$scans))

  for(scan_idx in 1:num_scans_to_process) {

    occult_data <- occult$scans[[scan_idx]]$params$OCCULT
    class(occult_data) <- "matrix"
    blockage_mask <- occult_data > blockage_threshold

    for(param_name in names(my_pvol$scans[[scan_idx]]$params)) {
      param_data <- my_pvol$scans[[scan_idx]]$params[[param_name]]
      class(param_data) <- "matrix"
      param_data[blockage_mask] <- NA
      class(param_data) <- c("param", "matrix", "array")
      my_pvol$scans[[scan_idx]]$params[[param_name]] <- param_data
    }
  }

  my_ppi = integrate_to_ppi(my_pvol, my_vp, res = 250, xlim = c(-100000, 100000), ylim = c(-100000, 100000))
  ppi_list <<- c(ppi_list, list(my_ppi))
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
#process vp

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

