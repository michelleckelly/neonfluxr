# `neonfluxr`
`neonfluxr` is a wrapper package for [`neonUtilities`](https://cran.r-project.org/web/packages/neonUtilities/index.html) that: 
1. Pulls [NEON](https://data.neonscience.org/) data for two-station stream metabolism modeling
2. Reformats that data for use with [`streamMetabolizer`](https://github.com/DOI-USGS/streamMetabolizer) or [`fluxr`](https://github.com/michelleckelly/fluxr)

## `request_NEON()`
`request_NEON()` uses [`neonUtilities::loadByProduct`](https://cran.r-project.org/web/packages/neonUtilities/refman/neonUtilities.html#loadByProduct) to request [water quality](https://data.neonscience.org/data-products/DP1.20288.001), [surface water temperature](https://data.neonscience.org/data-products/DP1.20053.001), [photosynthetically active radiation](https://data.neonscience.org/data-products/DP1.00024.001), [continuous discharge](https://data.neonscience.org/data-products/DP4.00130.001), and [barometric pressure](https://data.neonscience.org/data-products/DP1.00004.001) data products from the NEON API for a specified site and dates. 

`request_NEON()` then formats the output of `loadByProduct`, creating a single dataframe of dissolved oxygen, water temperature, light, water depth, discharge, and air pressure data at the upstream and downstream sensor stations.



## `clean_NEON()` 

Takes products returned from `request_NEON()` and checks for obviously erronious sensor data (i.e. things like discharge < 0), equal 15-minute breaks throughout the time series, and the amount of missing data 

Any missing data is filled using an ARIMA model (check out the documentation for [`auto.arima()`](https://www.rdocumentation.org/packages/forecast/versions/8.16/topics/auto.arima) for more info)

K600 and travel time between sensor stations are extrapolated for each time point using the relationship between field measurements and discharge 

Column names and formatting are changed to jive with `streamMetabolizer` and `fluxr`'s requirements for metabolism modeling

