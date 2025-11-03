# `neonfluxr`
`neonfluxr` is an R package for pulling & reformatting NEON (https://data.neonscience.org/) river data for two-station metabolism modeling. 

`neonfluxr` is intended for use in conjunction with the `fluxr` package (https://github.com/michelleckelly/fluxr).

## What does it do?

1. `request_NEON()`  
    - Pulls and formats water quality, temperature, PAR, nitrate, discharge, water chemistry, aquatic plant counts, and barometric pressure data products during the time period of choice into a single time-series dataframe
      
2. `clean_NEON()` 
    - Takes products returned from `request_NEON()` and checks for obviously erronious sensor data (i.e. things like discharge < 0), equal 15-minute breaks throughout the time series, and the amount of missing data 
    - Any missing data is filled using an ARIMA model (check out the documentation for [`auto.arima()`](https://www.rdocumentation.org/packages/forecast/versions/8.16/topics/auto.arima) for more info)
    - K600 and travel time between sensor stations are extrapolated for each time point using the relationship between field measurements and discharge 
    - Column names and formatting are changed to jive with `streamMetabolizer`'s requirements for metabolism modeling

