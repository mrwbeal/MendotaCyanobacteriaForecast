#
# This is a Shiny web application to build cyanobacteria biomass forecasts for Lake Mendota
# shiny::runGitHub('MendotaCyanobacteriaForecast', 'mrwbeal')

#Libraries
library(shiny)
library(fitdistrplus)
library(dataRetrieval) 
library(ggplot2)
library(dataRetrieval)
library(dplyr)
library(psych)
library(verification)
library(lubridate)
library(tidyverse)
library(pls)
library(raster)
library(sp)
library(shinythemes)
library(ncdf4)
library(RCurl)

#Data Biomass
basedata<-read.csv(file="https://raw.github.com/mrwbeal/MendotaCyanobacteriaForecast/main/JJA_forecast_github/CbPreds_9518.csv",header=T)
year = basedata$year
basedata<-basedata[,-c(1,2)]
cyano<-basedata$JJAcyano #Dependent variable

#Data Beaches
#basedata_o<-read.csv("~/Desktop/PhD/Forecast Coding/BeachPreds_0520.csv")
basedata_o<-read.csv("https://raw.github.com/mrwbeal/MendotaCyanobacteriaForecast/main/JJA_forecast_github/BeachPreds_0520.csv")[-1,]
Beaches = c("James Madison", "Tenney Park", "Warner")
Predictands = c("Days.closed", "Periods.closed")

aboveOut = c()
nearOut = c()
predOut = c()
beachOut = c()
metricOut = c()

#generate the zscore for each variable
n<-length(cyano)
mamdisc<-scale(basedata$dis,center = TRUE,scale = TRUE)
events<-scale(basedata$EE1.6,center = TRUE,scale = TRUE)
sst<-scale(basedata$sst,center = TRUE,scale = TRUE) 
phos<-scale(basedata$phos,center = TRUE,scale = TRUE) 
precip <-scale(basedata$April_precip,center = TRUE,scale = TRUE)
ones<-c(rep(1,length(cyano))) 


print("Building Model")
#Determine number of PCs with eigenvalues >1
pca<-matrix(c(events,sst,mamdisc,precip),nrow=length(cyano),ncol=4)
dependtpca<-cyano #returns every column of cyano but ncol=i
prim_drop<-prcomp(na.omit(pca))
coeffd<-prim_drop$rotation #egienvector
scored<-prim_drop$x #PC score
##now get eigen values
e_pcad<-eigen(cov(pca))
latentd<-e_pcad$values
pcn<-0
for (j in seq(1,2)){#determine how many PCs should be used
  if (latentd[j]>1){
    pcn=pcn+1;
  }
  else{
    pcn=pcn+0
  }
}

data<-basedata[,c("JJAcyano","EE1.6","April_precip","dis","sst")]

modelpreds <- c()
index = seq(1,length(cyano)) #set up the index

for (i in 1:nrow(data)) {
  ind<-index!=i
  train<-data[ind,]
  test <- data[!ind,]
  pcr_model<-pcr(JJAcyano~., data = train,scale =TRUE, validation = "CV")
  pred<-predict(pcr_model, test, ncomp = pcn)
  modelpreds<-rbind(modelpreds,c(pred,year[!ind]))
}

modelPCA_drop<-modelpreds[,1]

#Error metrics
rsme<-sqrt(mean(cyano-t(modelPCA_drop))^2)
err<-cyano-t(modelPCA_drop)
verr<-as.vector(err)
dist<-fitdist(verr,"norm") 

muhat<-dist$estimate[1]
sigmahat<-dist$estimate[2]
distl<-fitdist(verr,"logis")

lochat<-distl$estimate[1]
shat<-distl$estimate[2]

rfit<-function(n,r,p1,p2,fit){
  mR<-matrix(,ncol = n,nrow = r)
  for (j in seq(1,n)){
    for (i in seq(1,r)){
      mR[i,j]<-fit(1,p1,p2)#generate random normally distributed variables
    }
  }
  return(mR)
}

r=100
#For normal distribution
set.seed(7)
Rnorm<-rfit(n,r,muhat,sigmahat,rnorm)
yhat<-matrix(rep(modelPCA_drop, each = r),ncol=n,nrow=r)
yvarn<-yhat+Rnorm
yvarn[yvarn<0]<-0


buildForecast<-function(year) {
  ######## Download data for forecast ######
  #Get data for the year of interest

  start.date <- paste0(year,"-03-01")
  end.date <- paste0(year,"-05-31")
  
  #Discharge data
  discharge <- readNWISdv(siteNumbers="05427850",parameterCd = "00060",startDate = start.date,endDate = end.date)
  disc<-mean(discharge$X_00060_00003)
  
  #Extreme Events (>40 mm) data (Updated to 40 mm as of 2022) 
  #MRCC measures in inches (40 mm = 1.57 in)
  curlSetOpt(timeout = 200)
  drca = read.csv('https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/USW00014837.csv')
  drca=drca[,c(1:7)]
  drca$YR = as.numeric(substr(drca$DATE,1,4))
  drca$MN = as.numeric(substr(drca$DATE,6,7))
  
  prcp_daily = drca %>% subset(YR==year & MN>=3 & MN<=5)
  prcp_daily$mm=prcp_daily$PRCP/10
  
  #Extreme Events, greater than 40mm
  precip = sum(as.numeric(prcp_daily$mm>40),na.rm=T)
  
  #April Precipitation
  ap_precip=prcp_daily %>% filter(MN==4) %>% summarize(ap_precip=sum(mm,na.rm=TRUE))
  mm_to_inches <- function(mm) {
    inches <- mm / 25.4
    return(inches)
  }
  ap_precip=mm_to_inches(ap_precip$ap_precip)
  
  
  #Sea Surface Temperature
  #Get SST data
  replace_year <- function(url, new_year) {
    new_url <- gsub("\\d{6}", paste0("20",new_year), url)
    return(new_url)
  }
  
  #OpenDAP from IRI Columbia, get ERSSTv5
  urlhold<-'http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCDC/.ERSST/.version5/.sst/T/%28Mar%202023%29%28May%202023%29RANGEEDGES/dods'
  url=replace_year(urlhold,year)
  
  nc_data = nc_open(url)
  
  lon <- ncvar_get(nc_data, "X", verbose = F)
  lat <- ncvar_get(nc_data, "Y", verbose = F)
  t <- ncvar_get(nc_data, "T")
  sst <- ncvar_get(nc_data, "sst")
  nc_close(nc_data)
  
  sstm = apply(sst, c(1,2), mean)
  
  latm = (lat >= -20) & (lat <= -10)
  lonm = (lon >= 210) & (lon <= 235)
  
  SST=mean(sstm[lonm,latm],na.rm=T)
  
  #Put all predictors from current year together
  pca = rbind(basedata[,c("dis","EE1.6","sst","April_precip")],data.frame("dis"=disc,"EE1.6"=precip,"sst"=SST,"April_precip"=ap_precip))
  forecast = data.frame("EE1.6"=precip,"April_precip"=ap_precip,"dis"=disc,"sst"=SST) 
  
  
  #####Biomass Forecast Code####
  pred<-predict(pcr_model, forecast, ncomp = pcn)
  PredVal = pred[1]
  yhat19<-matrix(rep(PredVal, each = 100),ncol=1,nrow=100)
  set.seed(12)
  Rnorm19<-rfit(1,100,muhat,sigmahat,rnorm)
  
  yvar19<-yhat19+Rnorm19
  sim<-yvar19
  obs<-cyano
  
  #creating 3 categories
  lower<-quantile(cyano, prob=.33)
  upper<-quantile(cyano, prob=.67)
  
  #setting up categories
  B<-0
  N<-0
  A<-0
  
  i=length(sim)
  #Loop to determine the probability of each category
  for (i in seq(1,length(sim))){
    simhold<-sim[i,]
    if (simhold>upper){A=A+1}
    else if (simhold<lower){ B=B+1}
    else {N=N+1}
  }
  
  pf_B=B/100;
  pf_N=N/100;
  pf_A=A/100;
  
  forecast<-c(pf_B,pf_N,pf_A)
  high<-max(forecast)
  {if (high == pf_B){fcstrank = "B"}
    else if (high == pf_N) {fcstrank = "N"}
    else{fcstrank = "A"}}
  
  kdensobs<-density(obs)
  kdenssim<-density(sim)
  lowery<-approx(kdenssim$x,kdenssim$y,xout=lower)
  uppery<-approx(kdenssim$x,kdenssim$y,xout=upper)
  forecast_percent<-forecast * 100
  
  density_fig<-data.frame("Data"=c(basedata$JJAcyano,sim),"Group"=NA)
  density_fig[1:23,2] = "June-August Long-term Average"
  density_fig[24:124,2] = "Forecast"
  
  #Probability density graph of Seasonal, sub seasonal, and observations from 1995-2017
  density_plot<-ggplot(density_fig, aes(Data, group=fct_rev(Group))) + 
    geom_density(aes(x=Data, fill=Group, alpha=0.5)) +
    theme_bw() + 
    theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
    theme(legend.position="bottom")+ xlab("Average cyanobacteria abundance (mg/L)") + 
    ylab("Probability Density") + theme(legend.title=element_blank())+
    geom_vline(slope=0,xintercept=lower,linetype="dashed")+
    geom_vline(slope=0,xintercept=upper,linetype="dashed")+ guides(size = "legend", alpha="none") +
    scale_x_continuous(breaks = c(0,2.5,5,7.5,10,12.5), limits = c(0,12.5)) 
  
  density_plot=density_plot + 
    annotate(geom="text", x=6.5, y=0.25, label=paste0(forecast_percent[1],"% probability of Below Normal algae"),hjust=0,color="black")  +
    annotate(geom="text", x=6.5, y=0.23, label=paste0(forecast_percent[2],"% probability of Near Normal algae"),hjust=0,color="black")  +
    annotate(geom="text", x=6.5, y=0.21, label=paste0(forecast_percent[3],"% probability of Above Normal algae"),hjust=0)  +
    annotate(geom="text",x=0.75,y=0.3,label="Below",fontface="bold") +
    annotate(geom="text",x=3,y=0.3,label="Near",fontface="bold")+
    annotate(geom="text",x=5,y=0.3,label="Above",fontface="bold")
  
  predictor = c("Precip Events > 40mm", "Total April Precip (in)", "Avg. Discharge Hwy 113 (cfs)", "Avg. Equatorial Pacific SST (C)")
  lta<-c(round(mean(basedata$EE1.6),1),round(mean(basedata$April_precip),1),round(mean(basedata$dis),1),round(mean(basedata$sst),1))
  current = c(precip,ap_precip,disc,SST)
  current=round(current,1)
  tbl = data.frame("Predictor (Mar-May)"=predictor, "1995-2018 Avg." = lta, "2023"=current,check.names=FALSE)
  
  
  ################ Beach Forecast Code ####################
  for (Bch in Beaches) {
    for (Ptd in Predictands) {
      
      Beach = Bch #1 = JM, 2=T, 3=W
      Predictand = Ptd #1 = days closed, 2 = periods closed
      
      basedata<-basedata_o[basedata_o$Beach==Beach,]
      
      disc<-basedata[,"dis"] #Discharge
      precip<-basedata[,"EE1.6"] #Precip events >1.6"
      phos<-basedata[,"phos"] #Phosphorus
      sst<-basedata[,"sst"] #SSt
      closed<-basedata[,Predictand] #Dependent variable
      #generate the zscore for each variable
      mamdisc<-scale(disc,center = TRUE,scale = TRUE)
      events<-scale(precip,center = TRUE,scale = TRUE)
      phos<-scale(phos,center = TRUE,scale = TRUE)
      ones<-c(rep(1,length(closed))) 
      #build the independent matrix
      rhsvar<-matrix(c(ones,mamdisc,events),nrow=length(closed),ncol=3)
      n<-length(closed)
      m<-mean(closed)
      year = basedata$year
      
      
      index = seq(1,n) #set up the index
      betasPCA_drop<-matrix(,nrow = 3,ncol = n)#create empty matrices for loop
      modelPCA_drop<-matrix(,nrow=1,ncol=n)
      
      #Determine number of PCs with eigenvalues >1
      pca<-matrix(c(mamdisc,events),nrow=length(closed),ncol=4)
      dependtpca<-closed #returns every column of cyano but ncol=i
      prim_drop<-prcomp(na.omit(pca))
      coeffd<-prim_drop$rotation #egienvector
      scored<-prim_drop$x #PC score
      ##now get eigen values
      e_pcad<-eigen(cov(pca))
      latentd<-e_pcad$values
      pcn<-0
      for (j in seq(1,2)){#determine how many PCs should be used
        if (latentd[j]>1){
          pcn=pcn+1;
        }
        else{
          pcn=pcn+0
        }
      }
      
      data<-data.frame(closed,"dis"=mamdisc,"EE"=events)
      
      modelpreds <- c()
      index = seq(1,length(closed)) #set up the index
      
      for (i in 1:nrow(data)) {
        ind<-index!=i
        train<-data[ind,]
        test <- data[!ind,]
        pcr_model<-pcr(closed~., data = train,scale =TRUE, validation = "CV")
        pred<-predict(pcr_model, test, ncomp = pcn)
        modelpreds<-rbind(modelpreds,c(pred,year[!ind]))
      }
      
      modelPCA_drop<-modelpreds[,1]
      
      
      
      rsme<-sqrt(mean(closed-t(modelPCA_drop))^2)
      err<-closed-t(modelPCA_drop)
      verr<-as.vector(err)
      dist<-fitdist(verr,"norm") 
      muhat<-dist$estimate[1]
      sigmahat<-dist$estimate[2]
      
      rfit<-function(n,r,p1,p2,fit){
        mR<-matrix(,ncol = n,nrow = r)
        for (j in seq(1,n)){
          for (i in seq(1,r)){
            mR[i,j]<-fit(1,p1,p2)#generate random normally distributed variables
          }
        }
        return(mR)
      }
      
      
      #For normal distribution
      Rnorm<-rfit(n,100,muhat,sigmahat,rnorm)
      yhat<-matrix(rep(modelPCA_drop, each = 100),ncol=n,nrow=100)
      yvarn<-yhat+Rnorm
      yvarn[yvarn<0]<-0
      
      
      
      ################Predictions for current year#####################
  
      
      forecastB = data.frame("dis"=disc,"EE"=precip)
      
      pca<-rbind(basedata[,c("dis","EE")],forecastB)
      pca<-scale(pca,center = TRUE,scale = TRUE)
      pca<-as.data.frame(pca)
      
      forecastB=pca[nrow(pca),]
      
      pred<-predict(pcr_model, forecastB, ncomp = pcn)
      predval = pred[1]
      yhat19<-matrix(rep(PredVal, each = 100),ncol=1,nrow=100)
      set.seed(12)
      Rnorm19<-rfit(1,100,muhat,sigmahat,rnorm)
      
      yvar19<-yhat19+Rnorm19
      sim<-yvar19
      obs<-closed
      
      
      N<-0
      A<-0
      for (k in seq(1,length(yvar19))){
        if (yvar19[k]>m){
          A=A+1
        }
        else{
          N=N+1
        }
      }
      
      Npercent<-N/100
      Apercent<-A/100
      
      sim<-yvar19
      obs<-closed
      kdensobs<-density(obs)
      kdenssim<-density(sim)
      plot(kdensobs,xlim=c(1,25),ylim=c(0,.2),xlab="Beach closed") + lines(kdenssim,col="red")
      
      print(paste0(A,"% Above Normal"))
      print(paste0(N,"% Near Normal"))
      print(paste0(round(predval,2)," ",Predictand, " ", Beach))
      
      aboveOut=rbind(aboveOut,A)
      nearOut = rbind(nearOut,N)
      predOut = rbind(predOut,predval)
      beachOut = rbind(beachOut,Beach)
      metricOut = rbind(metricOut,Predictand)
      
    }
  }
  
  #Prediction Table
  beachtableinput = data.frame("Beach"=beachOut,"Metric"=metricOut,"Closure Prediction"=round(predOut,1),check.names = FALSE)
  rownames(beachtableinput) <- NULL
  beach_predictions<-data.frame(aboveOut,nearOut,predOut,beachOut,metricOut,check.names = FALSE)
  beach_predictions=data.frame("Percent"=rbind(aboveOut,nearOut),beachOut,metricOut,check.names=FALSE)
  beach_predictions$Category=substr(rownames(beach_predictions),1,1)
  beach_predictions<-beach_predictions %>% mutate(Category = str_replace_all(Category,c("N"="Near Normal","A"="Above Normal")))
  beach_predictions<-beach_predictions %>% mutate(metricOut = gsub("[.]"," ",metricOut))
  
  BeachForecastPlot=ggplot(data=beach_predictions,aes(fill=Category,y=Percent,x=beachOut)) + 
    geom_bar(stat="identity") + facet_grid(~metricOut) + xlab("")+ theme_bw() + scale_fill_manual(values=c("darkred","skyblue")) + 
    theme(text=element_text(size=14),legend.position="bottom")
  
  output=list(density_plot, tbl, forecast_percent, BeachForecastPlot,beachtableinput)
  return(output)
} 


biomass_cat<-function(percent) {
  
  if (percent[1]==max(percent)) {
    out="The Jun-Aug model predicts average cyanobacteria biomass will likely be below normal (CB<1.9 mg/L)"
  }
  
  if (percent[2]==max(percent)) {
    out="The Jun-Aug model predicts average cyanobacteria biomass will likely be near normal (1.9<CB<3.5 mg/L)"
  }
  
  if (percent[3]==max(percent)) {
    out="The Jun-Aug model predicts average cyanobacteria biomass will likely be above normal (CB>3.5 mg/L)"
  }
  return(out)
}




##### Define UI for application #####
ui <- fluidPage(theme = shinytheme("united"),

    # Application title
    titlePanel("Mendota June-August Cyanobacteria Biomass Forecast"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      
      sidebarPanel(
        helpText("Select a year to build a forecast"),
        
        selectInput(
          inputId =  "year", 
          label = "Select time period:", 
          choices = c("Select",2022:2023)
        ),
        
        DT::dataTableOutput("table")
      ),
      
        # Show a plot of the generated distribution
      mainPanel(
        img(src="wss_logo.png",height=90,width=450),
        
        tabsetPanel(type = "tabs",
                    tabPanel("Biomass Forecast", plotOutput("forecastPlot"),textOutput("biomass_category")),
                    tabPanel("Beach Forecast", plotOutput("beachPlot"), DT::dataTableOutput("beachtable")),
        )
      ) #end main panel
    
  )
)

##### Define server logic ####
server <- function(input, output) {
  
  
  observeEvent(input$year,{
      if (input$year=="Select") {
        
      }
      else{
      showModal(modalDialog("Downloading Data...", footer=NULL))
      fcst<-buildForecast(input$year)
      
      output$table <- DT::renderDataTable(DT::datatable({
        fcst[[2]]
      }))
      
      output$forecastPlot <- renderPlot({
        print(fcst[[1]])
        removeModal()
      })
      
      output$beachPlot <- renderPlot({
        print(fcst[[4]])
        removeModal()
      })
      
      output$beachtable <- DT::renderDataTable(DT::datatable({
        print(fcst[[5]],row.names=FALSE)
      }))
      
      
      output$biomass_category <- renderText({ 
        biomass_cat(fcst[[3]])
      })
      
      }
    })
  
  test <- reactive({
    input$year
  })
  
  #Stop the app timing out
  autoInvalidate <- reactiveTimer(10000)
  observe({
    autoInvalidate()
    cat(".")
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
