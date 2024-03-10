# MendotaCyanobacteriaForecast
Code repository for Lake Mendota seasonal cyanobacteria forecasting models.

March 2024 Update

App is available (app.R) to build season-ahead (JJA) cyanobacteria biomass and beach closing forecasts. The app uses USGS guage data API, GHCN data for Dane County Regional Airport, and NOAA ERSSTv5 to build forecasts for a chosen year. 


To call the app from GitHub, open R and run:
shiny::runGitHub('MendotaCyanobacteriaForecast', 'mrwbeal')
