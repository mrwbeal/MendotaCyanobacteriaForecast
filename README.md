# MendotaCyanobacteriaForecast
Code repository for Lake Mendota seasonal cyanobacteria forecasting models.

August 2023 Update

App is available (app.R) to build season-ahead (JJA) cyanobacteria biomass and beach closing forecasts. The app uses USGS guage data API and the ERA5 API to build forecasts for a chosen year. 

To call the app from GitHub, open R and run:
shiny::runGitHub('MendotaCyanobacteriaForecast', 'mrwbeal')
