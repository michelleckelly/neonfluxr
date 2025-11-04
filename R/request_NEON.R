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

  # Format water quality data
  wq_data <-
    waterqual$waq_instantaneous %>%
    select(siteID, horizontalPosition, startDateTime,
                  dissolvedOxygen, localDissolvedOxygenSat, pH, chlorophyll,
                  turbidity, specificConductance) %>%
    mutate(horizontalPosition = dplyr::recode(horizontalPosition,
                                              "101" = "upstream",
                                              "110" = "upstream",
                                              "100" = "upstream",
                                              "102" = "downstream",
                                              "200" = "downstream",
                                              "101" = "upstream",
                                              "S1" = "upstream",
                                              "S2" = "downstream"),
                  startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    rename(datetime_UTC = startDateTime,
           location = horizontalPosition, DO.obs_mgL = dissolvedOxygen,
           DO.sat_mgL = localDissolvedOxygenSat, chlorophyll_ugL = chlorophyll,
           turbidity_NTU = turbidity, specificConductance_uScm = specificConductance)

  if(!expanded){
    wq_data <-
      wq_data %>%
      select(-chlorophyll_ugL, -turbidity_NTU, -specificConductance_uScm)
  }

  # Format temperature data
  temp_data <-
    temp$TSW_5min %>%
    select(siteID, horizontalPosition, startDateTime, surfWaterTempMean) %>%
    mutate(horizontalPosition = recode(horizontalPosition,
                                       "101" = "upstream",
                                       "110" = "upstream",
                                       "100" = "upstream",
                                       "102" = "downstream",
                                       "200" = "downstream",
                                       "101" = "upstream",
                                       "S1" = "upstream",
                                       "S2" = "downstream"),
           startDateTime = lubridate::with_tz(startDateTime, tz = "UTC")) %>%
    rename(datetime_UTC = startDateTime, watertemp_C = surfWaterTempMean,
           location = horizontalPosition)

  # Format air pressure
  press_data <-
    barometer$BP_1min %>%
    select(siteID, horizontalPosition, startDateTime, staPresMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           horizontalPosition = recode(horizontalPosition,
                                       "101" = "upstream",
                                       "110" = "upstream",
                                       "100" = "upstream",
                                       "102" = "downstream",
                                       "200" = "downstream",
                                       "101" = "upstream",
                                       "S1" = "upstream",
                                       "S2" = "downstream")) %>%
    rename(datetime_UTC = startDateTime, pressure_kPa = staPresMean,
           location = horizontalPosition)

  # Format discharge and water depth
  if(!all(is.na(discharge))){
    discharge_data <-
      discharge$csd_continuousDischarge %>%
      select(siteID, endDate, stationHorizontalID, equivalentStage,
             maxpostDischarge) %>%
      mutate(endDate = lubridate::with_tz(endDate, tz = "UTC"),
             stationHorizontalID = recode(stationHorizontalID,
                                          "101" = "upstream",
                                          "110" = "upstream",
                                          "100" = "upstream",
                                          "102" = "downstream",
                                          "200" = "downstream",
                                          "132" = "downstream",
                                          "131" = "downstream",
                                          "101" = "upstream",
                                          "S1" = "upstream",
                                          "S2" = "downstream"),
             # convert discharge, measured as L/s, to m3/s
             maxpostDischarge = maxpostDischarge / 1000) %>%
      rename(datetime_UTC = endDate, depth_m = equivalentStage,
             discharge_m3s = maxpostDischarge,
             location = stationHorizontalID)
  }

  # Format light
  par_data <-
    par$PARPAR_1min %>%
    select(siteID, horizontalPosition, startDateTime, PARMean) %>%
    mutate(startDateTime = lubridate::with_tz(startDateTime, tz = "UTC"),
           horizontalPosition = recode(horizontalPosition,
                                       "101" = "upstream",
                                       "110" = "upstream",
                                       "100" = "upstream",
                                       "102" = "downstream",
                                       "200" = "downstream",
                                       "101" = "upstream",
                                       "S1" = "upstream",
                                       "S2" = "downstream")) %>%
    rename(datetime_UTC = startDateTime, par_umolm2s = PARMean,
           location = horizontalPosition)

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
               horizontalPosition = recode(horizontalPosition,
                                           "101" = "upstream",
                                           "110" = "upstream",
                                           "100" = "upstream",
                                           "102" = "downstream",
                                           "200" = "downstream",
                                           "101" = "upstream",
                                           "S1" = "upstream",
                                           "S2" = "downstream")) %>%
        rename(datetime_UTC = startDateTime, location = horizontalPosition,
               nitrate_umolL = surfWaterNitrateMean)
    }
    if(!all(is.na(waterchem))){
      wc_data <-
        waterchem$swc_externalLabDataByAnalyte %>%
        select(siteID, namedLocation, startDate, analyte, analyteConcentration,
               analyteUnits, sampleCondition) %>%
        mutate(namedLocation = stringr::str_extract(string = namedLocation,
                                                    pattern = "\\w\\d$"),
               namedLocation = recode(namedLocation,
                                      "101" = "upstream",
                                      "110" = "upstream",
                                      "100" = "upstream",
                                      "102" = "downstream",
                                      "200" = "downstream",
                                      "101" = "upstream",
                                      "S1" = "upstream",
                                      "S2" = "downstream")) %>%
        rename(location = namedLocation)
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

  # Grab average sensor longitude for upstream and downstream locations
  lat_upstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceLatitude[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "101")]
  )
  long_upstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceLongitude[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "101")]
    )
  elev_upstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceElevation[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "101")]
  )

  lat_downstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceLatitude[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "102")]
  )
  long_downstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceLongitude[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "102")]
  )
  elev_downstream <- mean(
    waterqual$sensor_positions_20288$locationReferenceElevation[
      str_detect(waterqual$sensor_positions_20288$HOR.VER, "102")]
  )

  # Output
  list(
    metabolism.parameters = params,
    no3 = if(exists("no3_data")) no3_data else NA,
    waterchemistry = if(exists("wc_data")) wc_data else NA,
    transects = if(exists("transect_data")) transect_data else NA,
    upstream = c(latitude = lat_upstream, longitude = long_upstream,
                 elevation = elev_upstream),
    downstream = c(latitude = lat_downstream, longitude = long_downstream,
                   elevation = elev_downstream)
  )

}
