library(sf)
library(suntools)
library(lubridate)

locations <- list(
  KCYS = list(
    code = 'KCYS',
    lon = -104.8061,
    lat = 41.15194
  ),
  KRIW = list(
    code = 'KRIW',
    lon = -108.4772,
    lat = 43.06611
  )
)

for (location in locations) {

  all_dates = vector("list", 2)
  datelist = c()
  years <- 2007:2019
  spring_period <- list(name="spring", start = "03-01", end = "06-15")
  fall_period <- list(name="fall", start = "08-15", end = "12-01")
  periods = list(spring_period, fall_period)

  generate_dates <- function(year, start, end) {
    seq(as.Date(paste(as.character(year), start, sep="-")), as.Date(paste(as.character(year), end, sep="-")), by="day")
  }

names(all_dates) <- sapply(periods, function(x) x$name)

  for (period in periods){
    datelist <- c()
    sunset_list <- c()
    for (year in years){
      dates <- as.character(generate_dates(year, period$start, period$end))
      datelist <- append(datelist, dates)
      for (date in dates){
        sunset_info <- sunriset(
                       matrix(c(location$lon, location$lat), nrow = 1),
                       as.POSIXct(date, tz = "UTC"),
                       direction='sunset',
                       POSIXct.out=TRUE)
        sunset_date = sunset_info$time + (3600 * 3)
        sunset_list = append(sunset_list, sunset_date)
      }
    }
    all_dates[[period$name]]$date <-   datelist
    all_dates[[period$name]]$sunset <- format(sunset_list, format="%Y-%m-%d %H:%M:%S")
  }

process_season <- function(season_data, season_name) {
  df <- data.frame(
    date = as.Date(season_data$date),
    sunset = as.POSIXct(season_data$sunset, format="%Y-%m-%d %H:%M:%S"),
    season = season_name,
    stringsAsFactors = FALSE
  )
  # Extract year from the date
  df$year <- format(df$date, "%Y")
  df$week <- week(df$date)
  df$radar <- location$code
  return(df)
}

all_seasons_df <- do.call(rbind, lapply(names(all_dates), function(season_name) {
  process_season(all_dates[[season_name]], season_name)
}))

write.csv(all_seasons_df, paste0(location$code, "_sunset_times.csv"), row.names = FALSE)

}



