######################################
#      LOAD REQUIRED LIBRARIES       #
######################################

library(dataRetrieval)
library(waterData)
library(tidyverse)
library(sf)
library(ncdf4)
library(lubridate)
library(nlme)
library(mgcv)
library(gt)

######################################
#      SOURCE HELPER FUNCTIONS       #
######################################

source('~/bighole-hydro/R/functions.R')

######################################
#    GET USGS STREMFLOW AND TEMP     #
######################################

#define the states of interest 
state_codes = c("MT")

# Initialize a list to store results
site_list = list()

# Loop through states and get sites with stream temp data
for (st in state_codes) {
  message(paste("Fetching sites in", st))
  sites = tryCatch(
    whatNWISdata(stateCd = st, parameterCd = "00010", service = "dv"),
    error = function(e) NULL
  )
  
  if (!is.null(sites) && nrow(sites) > 0) {
    site_list[[st]] = sites
  }
  
  Sys.sleep(0.5)  # Be polite to the API
}

# Combine all results into one data frame
# and then convert to SF object
spatial_meta =  bind_rows(site_list) |>
  as_tibble()  |> 
  distinct(site_no, .keep_all = TRUE) |>
  st_as_sf(coords = c('dec_long_va', 'dec_lat_va')) |>
  st_set_crs('EPSG:4269') |>
  st_transform('EPSG:4326')

# read in HUC8 dataset and pull out the big hole
big_hole = read_sf('https://data.climate.umt.edu/drought-indicators/fgb/watersheds.fgb') |>
  filter(NAME == 'Big Hole')

# compute spatial intersection with the bighole watershed
big_hole_sites = spatial_meta |>
  st_intersection(big_hole) |>
  filter(!site_no == '453136112420301') # this is a well, not a stream gage

# get the Hourly data:
# Initialize list for storage
temp_usgs_list = list()

# loop through sites to extract the data
for (site in  big_hole_sites$site_no) {
  message(paste("Getting hourly data for site", site))
  
  tryCatch({
    # Read unit (instantaneous) values
    temp_data = readNWISuv(siteNumbers = site,
                            parameterCd = c("00010", '00060'),
                            startDate = "1900-01-01",
                            endDate = Sys.Date()) # Running on 8/21/2025 - Updateable
    
    if (!is.null(temp_data) && nrow(temp_data) > 0) {
      temp_data$site_no = site
      temp_usgs_list[[site]] = temp_data
    }
  }, error = function(e) {
    message(paste("Failed for site", site, ":", e$message))
  })
  
  Sys.sleep(0.3)
}

######################################
#    GET StAGE STREMFLOW AND TEMP    #
######################################

# import StAGE meta data for sites in the Big Hole from:
# https://gis.dnrc.mt.gov/apps/stage/gage-report/
StAGE_meta = read_csv('~/bighole-hydro/StAGE-data/site-meta.csv') |>
  st_as_sf(coords = c('lon', 'lat')) |>
  st_set_crs('EPSG:4326') |>
  mutate(network = 'StAGE') |>
  select(network, site_no = site_id,
         station_nm = name)

# import the data, clean and bind
all_StAGE_data = list.files("~/bighole-hydro/StAGE-data/raw", full.names = TRUE) |>
  purrr::map(readr::read_csv, skip = 1) |>
  purrr::map(clean_stage_tbl) |>
  dplyr::bind_rows() |>
  rename(dateTime = Timestamp) |>
  select(site_no, dateTime, temp, discharge) |>
  drop_na(temp, discharge)

######################################
#      MERGE StAGE AND USGS DATA     #
######################################

# merge site metadata
all_site_meta = big_hole_sites |>
  select(network = agency_cd, site_no, station_nm) |>
  bind_rows(StAGE_meta)

# plot sites with HUC8 boundry and flowlines
plot_big_hole_sites_nhd(
  sites   = all_site_meta,
  huc8    = "10020004",
  label   = TRUE,
  label_col = "site_no",
  save_path = "~/bighole-hydro/figs/big_hole_sites_nhd.png",
  flowlines_url = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer/3"
)

# Combine all into one data frame
all_hourly_temp = bind_rows(temp_usgs_list) |>
  as_tibble() |>
  rename(temp = X_00010_00000,
         discharge = X_00060_00000) |>
  select(site_no, dateTime, temp, discharge) |>
  drop_na(temp, discharge) |>
  bind_rows(all_StAGE_data) |>
  mutate(hour = lubridate::hour(dateTime),
         date = as.Date(dateTime)) |> 
  group_by(site_no, date, hour) |>
  summarise(max_hourly_temp = max(temp),
            median_hourly_discharge = median(discharge)) |>
  ungroup()

# Identify "complete" May–Oct seasons
# at least 80% complete record
completeness_lambda = 0.8

complete_seasons =
  all_hourly_temp |>
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::filter(date >= as.Date(paste0(year, "-05-01")) &
                  date <= as.Date(paste0(year, "-10-31"))) |>
  dplyr::group_by(site_no, year) |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_non_na = sum(!is.na(max_hourly_temp)),
    expected_hours = 24L * as.integer(as.Date(paste0(year, "-10-31")) -
                                        as.Date(paste0(year, "-05-01")) + 1L),
    .groups = "drop"
  ) |>
  dplyr::filter(n_rows > expected_hours*completeness_lambda, n_non_na > expected_hours*completeness_lambda) |>
  dplyr::select(site_no, year) |>
  mutate(valid_site_year = paste0(site_no, '_', year))

# compute exceedance 
exceedance_hours = all_hourly_temp |>
  mutate(year = lubridate::year(date)) |>
  group_by(site_no, year) |>
  summarise(cumulative_hours_21 = sum(max_hourly_temp >= 21),
            cumulative_hours_25 = sum(max_hourly_temp >= 25)) |>
  mutate(valid_site_year = paste0(site_no, '_', year)) |>
  filter(valid_site_year %in% complete_seasons$valid_site_year) |>
  select(-valid_site_year) |>
  left_join(all_site_meta |>
              select(site_no, station_nm) |>
              st_drop_geometry(), by = 'site_no')|>
  add_count(site_no, name = "n_obs") |>
  filter(n_obs > 3)

# compute GAM for each exceedance facet
smooth_21 = fit_per_station_gam(exceedance_hours |>
                                  rename(cumulative_hours = cumulative_hours_21), by = 0.1)

smooth_25 = fit_per_station_gam(exceedance_hours |>
                                  rename(cumulative_hours = cumulative_hours_25) , by = 0.1)

# Using precomputed smoothers:
plot_21 = plot_exceedance(
  df = exceedance_hours,
  y_col = "cumulative_hours_21",
  smooth_df = smooth_21,
  facet_var = "station_nm",
  save_path = "~/bighole-hydro/figs/big_hole_temp_21C.png"
)

plot_25 = plot_exceedance(
  df = exceedance_hours,
  y_col = "cumulative_hours_25",
  smooth_df = smooth_25,
  facet_var = "station_nm",
  save_path = "~/bighole-hydro/figs/big_hole_temp_25C.png"
)

#####################################
#   EXCEEDANCE PROBABILITY MODEL    #
#####################################
break_point = as.Date("2017-01-01")

# get daily min and max air temp for the sites
# takes a few min to retrieve from GRIDMET server
daily_temps =
  get_gridmet_daily_tmmx_tmmn_sf(all_site_meta, out_units = "C")

# compute the daily max temperature for each site
# filter for days with complete data to ensure correct
# max temperature computation
# 399,500 compete hourly records as of 8/20/2025
daily_max_complete =
  all_hourly_temp |>
  dplyr::group_by(site_no, date) |>
  dplyr::summarise(
    n_temp = length(max_hourly_temp),
    n_flow = length(median_hourly_discharge),
    daily_max_stream_temp = max(max_hourly_temp),
    daily_median_discharge = median(median_hourly_discharge),
    .groups = "drop"
  ) |>
  filter(n_temp == 24 & n_flow == 24) |>
  select(-c(n_temp, n_flow)) 

#compute the site_years that exist both before and after 2017
symetric_site_years = 
  all_hourly_temp |>
  mutate(site_year = paste0(site_no, '_', year(date))) |>
  filter(site_year %in% unique(complete_seasons$valid_site_year)) |>
  group_by(site_no) |> 
  summarise(
    min_year = min(year(date), na.rm = TRUE),
    max_year = max(year(date), na.rm = TRUE)
  ) |>
  filter(min_year < year(break_point) & max_year > year(break_point))

# merge the datasets for analysis
merged = daily_max_complete |>
  left_join(daily_temps |>
  select(site_no, date, tmmx, tmmn)) |>
  mutate(site_year = paste0(site_no,'_',year(date))) 
  
# compute the GAMM for both all data and for the subset of gages that meet
# the symetric_site_years critearia above
gamm_all_data = fit_exceedance_bam(merged, threshold = 21, break_date = break_point)
gamm_filtered_data = fit_exceedance_bam(
  merged  |>
    #filter for the season of interest
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::filter(date >= as.Date(paste0(year, "-05-01")) &
                  date <= as.Date(paste0(year, "-10-31")) &
                  #filter for compltete site_years as defined above using the 80% threshold
                  site_no %in% symetric_site_years$site_no), 
  threshold = 21, 
  break_date = break_point
  )

# convert the table to an image to describe results
make_exceedance_fixed_effects_table(
  fit       = gamm_all_data$fit,
  threshold = 21,
  html_path = "~/bighole-hydro/figs/exceedance_fixed_effects_all_data.html",
  png_path  = "~/bighole-hydro/figs/exceedance_fixed_effects_all_data.png",
  interp_text = interpret_exceedance_gamm(gamm_all_data$fit, break_date = break_point)
)

make_exceedance_fixed_effects_table(
  fit       = gamm_filtered_data$fit,
  threshold = 21,
  html_path = "~/bighole-hydro/figs/exceedance_fixed_effects_symetric.html",
  png_path  = "~/bighole-hydro/figs/exceedance_fixed_effects_symetric.png",
  interp_text = interpret_exceedance_gamm(gamm_filtered_data$fit, break_date = break_point)
)

#####################################
#   PLOT EXCEEDANCE MODEL RESULTS   #
#####################################

# By-site (faceted) predictions at median flow & tmmn, period = most common:
plot_exceedance_vs_tmmx(gamm_all_data$fit, gamm_all_data$dataset, mode = "by_site", 
                        save_path = '~/bighole-hydro/figs/exceedance_model_results_by_site.png')

# Pre vs Post by-site (colored), at medians:
plot_exceedance_vs_tmmx(gamm_all_data$fit, gamm_all_data$dataset, mode = "prepost_by_site", 
                        save_path = '~/bighole-hydro/figs/exceedance_model_results_prepost_by_site.png')

# Single fixed-effects curve (random effects excluded):
plot_exceedance_vs_tmmx(gamm_all_data$fit, gamm_all_data$dataset, mode = "fixed", 
                        save_path = '~/bighole-hydro/figs/exceedance_model_results_fixed.png')

# Single fixed-effects pre/post overlay:
plot_exceedance_vs_tmmx(gamm_all_data$fit, gamm_all_data$dataset, mode = "fixed_prepost", 
                        save_path = '~/bighole-hydro/figs/exceedance_model_results_fixed_prepost.png')

#####################################
#  SUMMERISE RESULTS WITH EXAMPLES  #
#####################################

summary_all = summarize_exceedance_gamm(gamm_all_data$fit, gamm_all_data$dataset)
summary_all

#####################################
#           COMPUTE TRENDS          #
#####################################

plot_airtemp_trends(daily_temps, save_path = "~/bighole-hydro/figs/airtemp_trends.png")

#for time periods with temp data
plot_streamflow_trends(
  merged,
  months = 5:10,
  end_year = 2024,
  save_path = "~/bighole-hydro/figs/flow_trends.png"
)

#For USGS gages with long records
usgs_df = dplyr::bind_rows(temp_usgs_list) |> 
  tibble::as_tibble()

p_flow_usgs = plot_streamflow_trends_usgs(
  usgs_raw   = usgs_df,
  months     =  5:10,
  end_year   = 2024,
  save_path  = "~/bighole-hydro/figs/flow_trends_usgs.png"
)

plot_exceedance_vs_airtemp_from_bases(
  daily_temps, exceedance_hours,
  response = "cumulative_hours_21",
  save_path = "~/bighole-hydro/figs/exceedance_vs_airtemp_21.png"
)

res = analyze_seasonal_exceedance(
  daily_temps         = daily_temps,
  daily_max_complete  = daily_max_complete,
  exceedance_hours    = exceedance_hours,
  response            = "cumulative_hours_21",
  hold_at             = "overall_median",
  save_prefix         = "~/bighole-hydro/figs/seasonal_exceedance",
  height = 7, caption_width = 90
)

# Start from the table with site/year + metrics

plot_z_anomalies_by_site(
  res$data,
  save_path = "~/bighole-hydro/figs/z_anomalies_by_site.png"
)
