#' Request data products from NEON API for two-station metabolism modeling
#'
#' \code{request_NEON} communicates with NEON API to request DO,
#'   water temperature, PAR, discharge, depth, and water quality data for the
#'   time period and monitoring station(s) of interest (using NEONScience's
#'   \code{neonUtilities::loadByProduct}). DO
#'   percent saturation is corrected for local elevation
#'   using NEONScience's \code{localPressureDO::calcBK_eq}.
#'
#'   \code{request_NEON} additionally requests salt and conductivity slug data
#'   necessary for reaeration rate (K) calculations, and displays an
#'   interactive GUI to the user to identify peaks in slug transport using
#'   NEONScience's \code{reaRate::def.calc.reaeration}. Afterwards,
#'   K600 values are calculated based on the user's input.
#'
#'   Returned is a list containing:
#'   (1) \code{data}: a formatted dataframe with all raw data necessary for
#'   the user to model two-station stream metabolism,
#'   (2) \code{k600_clean}:describe here,
#'   (3) \code{k600_fit}: describe here,
#'   (4) \code{k600_expanded}: describe here.
#'
#' @importFrom magrittr %>%
#'
#' @param sitecode character string specifying 4-letter NEON site code to
#'   request data from (ex. `"HOPB"`). Can be more than one site
#'   (ex.`c("HOPB", "BLDE")` but be warned data pull will take longer)
#' @param startdate YYYY-MM character string defining start year and month for
#'    data request
#' @param enddate YYYY-MM character string defining end year and month for
#'    data request
#'
#' @seealso \url{https://github.com/NEONScience/NEON-utilities/neonUtilities/}
#'    for details on \code{neonUtilities} package
#'
#' @export
request_NEON <- function(sitecode, startdate, enddate, interval = 15, expanded = FALSE){
  # Input parameters for metabolism modeling
  params <- c("DP1.20288.001", "DP1.20053.001", "DP1.00024.001", "DP4.00130.001",
              "DP1.00004.001")
  names <- c("waterqual", "temp", "par", "discharge", "barometer")
  names(params) <- names

  if(expanded){
    params <- c(params, "DP1.20033.001", "DP1.20093.001", "DP1.20072.001")
    names(params) <- c(names, "no3", "waterchem", "aquatictransects")
  }

  # Pull NEON data from api
  for (i in seq_along(params)){
    dpID <- params[i]

    # Pull NEON data for parameter, saving the data to the variable name from
    # names(params)
    # Pull basic only to save time
    assign(names(dpID),
           tryCatch({
             neonUtilities::loadByProduct(dpID = dpID, site = sitecode,
                                          startdate = startdate,
                                          enddate = enddate,
                                          package = "basic",
                                          check.size = F)
             },
             error = function(err){
               NA
               }))
  }

  # To get the distance between S1 and S2, need to pull geomorphology data
  # Geomorphology surveys are only done once every 5 years, so pull all the
  # data that's available, rather than just for the period of interest
  geomorph <-
    neonUtilities::loadByProduct(dpID = "DP4.00131.001",
                                 site = sitecode,
                                 package = "basic",
                                 check.size = F)
  # Format geomorphology data
  geo_data <-
    geomorph$geo_surveySummary %>%
    filter(surveyBoutTypeID == "geomorphology") %>%
    select(siteID, eventID, startDate, meanBankfullWidth, S1habitatID,
           S2habitatID, sensorSetLength) %>%
    rename(datetime_UTC = startDate,
           habitatID_upstream = S1habitatID,
           habitatID_downstream = S2habitatID,
           distance_between_sensors_m = sensorSetLength)

  # Format water quality data
  wq_data <-
    waterqual$waq_instantaneous %>%
    select(siteID, horizontalPosition, startDateTime,
                  dissolvedOxygen, localDissolvedOxygenSat, pH, chlorophyll,
                  turbidity, specificConductance) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           startDateTime = round_date(startDateTime, unit = "minute")) %>%
    rename(datetime_UTC = startDateTime,
           location = horizontalPosition, DO.obs_mgL = dissolvedOxygen,
           DO.sat_mgL = localDissolvedOxygenSat, chlorophyll_ugL = chlorophyll,
           turbidity_NTU = turbidity,
           specificConductance_uScm = specificConductance)
  wq_data$location <-
    str_replace_all(wq_data$location,
                    c("100|101|110|S1" = "upstream",
                      "102|200|S2" = "downstream"))

  if(!expanded){
    wq_data <-
      wq_data %>%
      select(-chlorophyll_ugL, -turbidity_NTU, -specificConductance_uScm)
  }

  # Format temperature data
  temp_data <-
    temp$TSW_5min %>%
    select(siteID, horizontalPosition, startDateTime, surfWaterTempMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           startDateTime = round_date(startDateTime, unit = "minute")) %>%
    rename(datetime_UTC = startDateTime, watertemp_C = surfWaterTempMean,
           location = horizontalPosition)
  temp_data$location <-
    str_replace_all(temp_data$location,
                    c("100|101|110|S1" = "upstream",
                      "102|200|S2" = "downstream"))


  # Format air pressure
  press_data <-
    barometer$BP_1min %>%
    select(siteID, horizontalPosition, startDateTime, staPresMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           startDateTime = round_date(startDateTime, unit = "minute")) %>%
    rename(datetime_UTC = startDateTime, pressure_kPa = staPresMean,
           location = horizontalPosition)
  press_data$location <-
    str_replace_all(press_data$location,
                    c("100|101|110|S1" = "upstream",
                      "102|200|S2" = "downstream"))

  # Format discharge and water depth
  if(!all(is.na(discharge))){
    discharge_data <-
      discharge$csd_continuousDischarge %>%
      select(siteID, endDate, stationHorizontalID, #equivalentStage,
             maxpostDischarge) %>%
      mutate(endDate = lubridate::with_tz(endDate, tz = "UTC"),
             endDate = round_date(endDate, unit = "minute"),
             # convert discharge, measured as L/s, to m3/s
             maxpostDischarge = maxpostDischarge / 1000) %>%
      rename(datetime_UTC = endDate, #depth_m = equivalentStage,
             discharge_m3s = maxpostDischarge,
             location = stationHorizontalID)
    discharge_data$location <-
      str_replace_all(discharge_data$location,
                      c("100|101|110|S1" = "upstream",
                        "102|131|132|200|S2" = "downstream"))
  }

  # Format light
  par_data <-
    par$PARPAR_1min %>%
    select(siteID, horizontalPosition, startDateTime, PARMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           startDateTime = round_date(startDateTime, unit = "minute")) %>%
    rename(datetime_UTC = startDateTime, par_umolm2s = PARMean,
           location = horizontalPosition)
  par_data$location <-
    str_replace_all(par_data$location,
                    c("100|101|110|S1" = "upstream",
                      "102|200|S2" = "downstream"))

  # Cut to just 15-min measurements
  if(interval == 15){
    wq_data <- filter(wq_data, minute(datetime_UTC) %in% c(0, 15, 30, 45))
    temp_data <- filter(temp_data, minute(datetime_UTC) %in% c(0, 15, 30, 45))
    press_data <- filter(press_data, minute(datetime_UTC) %in% c(0, 15, 30, 45))
    discharge_data <- filter(discharge_data, minute(datetime_UTC) %in% c(0, 15, 30, 45))
    par_data <- filter(par_data, minute(datetime_UTC) %in% c(0, 15, 30, 45))
  }

  if(expanded){
    if(!all(is.na(no3))){
      # Format nitrate data
      no3_data <-
        no3$NSW_15_minute %>%
        select(siteID, horizontalPosition, startDateTime, surfWaterNitrateMean) %>%
        mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
               startDateTime = round_date(startDateTime, unit = "minute")) %>%
        rename(datetime_UTC = startDateTime, location = horizontalPosition,
               nitrate_umolL = surfWaterNitrateMean)
      no3_data$location <-
        str_replace_all(no3_data$location,
                        c("100|101|110|S1" = "upstream",
                          "102|200|S2" = "downstream"))
    }
    if(!all(is.na(waterchem))){
      wc_data <-
        waterchem$swc_externalLabDataByAnalyte %>%
        select(siteID, namedLocation, startDate, analyte, analyteConcentration,
               analyteUnits, sampleCondition) %>%
        mutate(namedLocation = stringr::str_extract(string = namedLocation,
                                                    pattern = "\\w\\d$")) %>%
        rename(location = namedLocation)
      wc_data$location <-
        str_replace_all(wc_data$location,
                        c("100|101|110|S1" = "upstream",
                          "102|200|S2" = "downstream"))
    }
    if(!all(is.na(aquatictransects))){
      transect_data <-
        aquatictransects$apc_pointTransect %>%
        select(siteID, namedLocation, pointNumber, transectDistance,
               collectDate, habitatType, substrate, remarks) %>%
        group_by(siteID, namedLocation, pointNumber) %>%
        arrange(transectDistance, .by_group = TRUE)

      ##NOTE this is for calculating percent cover - build this in to another function
      #percCover_full <-
      # coverData %>%
      #dplyr::count(substrate) %>%
      #dplyr::rename(transect = pointNumber, observationCount = n) %>%
      #dplyr::mutate(sumObservations = sum(observationCount),
      #             percentCover = observationCount / sumObservations * 100)

      #percCover_summary <-
      # percCover_full %>%
      #dplyr::select(transect, substrate, percentCover) %>%
      #tidyr::pivot_wider(names_from = substrate, values_from = percentCover) %>%
      #dplyr::rename(coarseWoodyDebris = CWD, leafLitter = `leaf litter`) %>%
      #replace(is.na(.), 0) %>%
      #dplyr::ungroup(transect) %>%
      #dplyr::select(-transect) %>%
      #dplyr::summarise_all(mean)
    }
  }

  # Join metabolism parameters into params dataframe
  params <-
    full_join(wq_data, temp_data, by = c("siteID","location", "datetime_UTC")) %>%
    full_join(., press_data, by = c("siteID","location", "datetime_UTC")) %>%
    full_join(., discharge_data, by = c("siteID","location", "datetime_UTC")) %>%
    full_join(., par_data, by = c("siteID","location", "datetime_UTC"))

  # Grab sensor positions during study period
  loc_data <-
    waterqual$sensor_positions_20288 %>%
    mutate(positionStartDateTime = ymd_hms(positionStartDateTime),
           positionEndDateTime = ymd_hms(positionEndDateTime)) %>%
   # filter(positionStartDateTime <= min(params$datetime_UTC)) %>%
    #filter(is.na(positionEndDateTime) |
     #        positionEndDateTime >= max(params$datetime_UTC)) %>%
    select(siteID, HOR.VER, sensorLocationID, sensorLocationDescription,
           positionStartDateTime, positionEndDateTime,
           locationReferenceLatitude, locationReferenceLongitude,
           locationReferenceElevation) %>%
    rename(latitude = locationReferenceLatitude,
           longitude = locationReferenceLongitude,
           elevation_m = locationReferenceElevation)

  # Output
  list(
    metabolism.parameters = params,
    no3 = if(exists("no3_data")) no3_data else NA,
    waterchemistry = if(exists("wc_data")) wc_data else NA,
    transects = if(exists("transect_data")) transect_data else NA,
    locations = loc_data,
    geomorphic = geo_data
  )

}
