library(bioRad)
library(terra)
library(stringr)

---
#helper functions

clean_pvol <- function(my_pvol) {
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

closest_times <- read.csv("pvol/KRIW_closest_times.csv")
df_filtered <- subset(closest_times, closest_file != "")

unique_year_week_df <- unique(data.frame(year = closest_times$year, week = closest_times$week))
unique_year_week_list <- split(unique_year_week_df, seq(nrow(unique_year_week_df)))
filtered_list <- Filter(function(x) x$year == 2019, unique_year_week_list)

df_subset = subset(closest_times, year == filtered_list[[23]]$year & week == filtered_list[[23]]$week)
df_subset$closest_file <- gsub("_MDM", "", df_subset$closest_file)

----
#time full script

start_time <- Sys.time()

ppi_list = list()
generate_ppis <- function(closest_file) {

  input_pvol_path <- paste0("pvol/", radar, "/", closest_file)
  input_vp_path <- paste0("vp/", radar, "/", closest_file, ".h5")

  my_vp = read_vpfiles(input_vp_path)

  apply_mistnet(input_pvol_path) -> my_pvol
  my_pvol = clean_pvol(my_pvol)

  my_ppi = integrate_to_ppi(my_pvol, my_vp, res = 250, xlim = c(-100000, 100000), ylim = c(-100000, 100000))
  ppi_list <<- c(ppi_list, list(my_ppi))
}

lapply(df_subset$closest_file, generate_ppis)
my_composite = composite_ppi(ppi_list, param = "VID", method = "mean", res=250)

r <- rast(my_composite$data)

output_tif_file <- paste0("tif/KRIW/",
                          radar, "_",
                          paste(filtered_list[[23]]$year,
                                filtered_list[[23]]$week, sep="_"),
                          ".tif")

terra::writeRaster(r, output_tif_file, overwrite=TRUE)

end_time <- Sys.time()

---
#process vp

file_dir <- "/Users/at744/clo/bigbird/wgfd/pvol/KRIW"
output_dir <- "/Users/at744/clo/bigbird/wgfd/vp/KRIW"
radar <- "KRIW"

files <- list.files(file_dir, pattern = "KRIW2019", full.names = TRUE)

for (file in files) {
  message("Processing file: ", file)

  output_file <- paste0(output_dir, "/", gsub(".*(KRIW2019.+)$", "\\1.h5", basename(file)))

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





