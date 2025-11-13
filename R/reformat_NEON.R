#' Clean NEON metabolism parameter data and prep for two-station stream metabolism modeling
#'
#' @importFrom streamMetabolizer convert_UTC_to_solartime
#' @importFrom fluxr o2sat
#'
#' @param data dataset
#' @param startdate start date for filtering
#' @param enddate end date for filtering
#' @param calculate_light logical, should the function calculate light from geographic position
#' @param single_station_vars character string, which variables were measured at only one location
#'
#' @return dataframe of formatted parameter data
#'
#' @export
reformat_NEON <- function(data, startdate, enddate, 
                          calculate_light = TRUE,
                          single_station_vars){
  ## Grab sensor locations #####################################################
  loc <- data$locations
  
  # If there's no position end date, put in todays date
  loc$positionEndDateTime[is.na(loc$positionEndDateTime)] <- 
    with_tz(Sys.time(), "UTC")
  
  # Grab sensor locations that overlap with the date interval
  sensorloc <-
    loc %>%
    mutate(timeInterval = 
             interval(positionStartDateTime, positionEndDateTime)) %>%
    filter(interval(startdate, enddate) %within% timeInterval) %>%
    group_by(HOR.VER) %>%
    summarise(latitude = mean(latitude), 
              longitude = mean(longitude),
              elevation_m = mean(elevation_m)) %>%
    rename(location = HOR.VER)
  
  # Change location names
  sensorloc$location <- 
    str_replace_all(sensorloc$location, c("101.100" = "upstream",
                                          "102.100" = "downstream"))
  
  # Join lat longs to the metabolism parameters df
  cleaned <- 
    full_join(data$metabolism.parameters, sensorloc, by = "location") 
  
  ## Solar time ################################################################
  # Using longitude column, calculate solar time
  cleaned$solar.time <- 
    streamMetabolizer::convert_UTC_to_solartime(cleaned$datetime_UTC,
                                                cleaned$longitude)
  
  # Make sure solar.time doesn't have any hanging seconds (i.e. 00:00:01)
  cleaned$solar.time <- lubridate::round_date(cleaned$solar.time, "minute")
  
  ## Filter data to interval of interest #######################################
  filtint <- interval(startdate, enddate)
  
  cleaned <- filter(cleaned, solar.time %within% filtint)
  
  ## Calculated light ##########################################################
  if(calculate_light){
    cleaned$light.calc <- 
      streamMetabolizer::calc_light(cleaned$solar.time, 
                                    cleaned$latitude, 
                                    cleaned$longitude)
  }
  
  ## Remove variables that were not measured at a sensor location ##############
  nosensor <- function(df, var){
    # Subset variable from big df, and remove na observations
    obs <- na.omit(select(df, c(solar.time, all_of(var))))
    
    # Remove the sensor variable from the big df
    df <- select(df, -all_of(var))
    
    # Merge new column for sensor variable back into big df
    df <- full_join(df, obs, relationship = "many-to-many", by = "solar.time")
    return(df)
  }
  
  # Loop through single_station_vars
  for(i in seq_along(single_station_vars)){
    currentvar <- single_station_vars[i]
    cleaned <- nosensor(cleaned, currentvar)
  }
  
  
  ## Calculate DO saturation ###################################################
  cleaned$DO.sat_mgL <- 
    fluxr::o2sat(cleaned$watertemp_C, cleaned$pressure_kPa, 
                 pressUnits = "kPa", 
                 outUnits = "mg/L")
  
  ## Change near-0 PAR values to actually 0 ####################################
  cleaned$par_umolm2s <-
    case_when(cleaned$par_umolm2s <= 0 & cleaned$par_umolm2s > -1 ~ 0,
              cleaned$par_umolm2s > 0 ~ cleaned$par_umolm2s)
  
  ## Rename columns to be fluxr-friendly #######################################
  cleaned <- 
    rename(cleaned, DO.obs = DO.obs_mgL, DO.sat = DO.sat_mgL,
           temp.water = watertemp_C, discharge = discharge_m3s, 
           pressure = pressure_kPa, light.data = par_umolm2s) %>%
    arrange(solar.time) %>%
    select(siteID, solar.time, location, DO.obs, DO.sat, temp.water, 
           pressure, discharge, light.data, any_of("light.calc"))
  
  ## Return data to user #######################################################
  return(cleaned)
}