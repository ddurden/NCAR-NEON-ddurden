######################################################################################################################

#' @title Workflow to download Ameriflux data and format for CLM

#' @author
#' David Durden \email{eddy4R.info@gmail.com}


#' @description
#' Workflow. Read Ameriflux data from API and format for CLM

#' @param Currently none

#' @return Currently none

#' @references
#' License: Terms of use of the NEON FIU algorithm repository dated 2015-01-16.


#' @keywords eddy-covariance, energy balance

#' @examples Currently none

#' @seealso Currently none

# changelog and author contributions / copyrights
#   David Durden (2024-06-04) - original creation
######################################################################################################################

#Install amerifluxr
if(!require(devtools)){install.packages("devtools")}
devtools::install_github("chuhousen/amerifluxr")

#Call the R HDF5 Library
packReq <- c("corrplot", "neonUtilities", "dplyr", "readr",'lubridate', "ggplot2","ggpubr","Hmisc","tidyverse","ggpmisc","jsonlite","factoextra")

#Install and load all required packages
lapply(packReq, function(x) {
  print(x)
  if(require(x, character.only = TRUE) == FALSE) {
    install.packages(x)
    library(x, character.only = TRUE)
  }})


#Download directory
DirDnld <- c(tempdir(),"/home/ddurden/eddy/tmp/clmAmf")[2]
DirExtr <- paste0(DirDnld,"/DirExtr")
#Check if directory exists and create if not
if(!dir.exists(DirExtr)) dir.create(DirExtr, recursive = TRUE)

#The version data for the FP standard conversion processing
ver <- paste0("vAmf/",format(Sys.time(), "%Y%m%d"))
DirOutBase <-paste0("~/eddy/data/CLM/",ver)




#Type of data
DataType <- "FULLSET"
#Time aggregation
TimeAgg <- c("HH","DD","WW","MM","YY")[5]
#NEON site
Site <- c("US-xAE","US-xBR","US-xCL","US-xCP","US-xDC","US-xDL","US-xDS",
          "US-xGR","US-xHA","US-xHE","US-xJE","US-xKA","US-xKZ", "US-xMB",
          "US-xML","US-xNG","US-xNQ","US-xRM","US-xSB","US-xSE","US-xSR",
          "US-xST","US-xTA","US-xTR","US-xUK","US-xUN","US-xYE")[2]

#Append the site to the base output directory
if(DirOutBase == "tmp") DirOutBase <- tempdir()
DirOut <- paste0(DirOutBase, "/", Site)
DirOutAtm <- paste0(DirOutBase, "/", Site, "/atm")
DirOutEval <- paste0(DirOutBase, "/", Site, "/eval")

NEONSiteMeta <- read.csv("/home/ddurden/eddy/code/random_scripts/EB_AGU/NEON_Field_Site_Metadata_20201204.csv")

#Get file via amerifluxr
fileInp <- amerifluxr::amf_download_fluxnet(
  user_id = "ddurden",
  user_email = "daviddurden27@gmail.com",
  site_id = Site,
  data_product = "FLUXNET",
  data_variant = DataType,
  agree_policy = TRUE,
  intended_use = "other",
  intended_use_text = "testing download",
  out_dir = DirDnld)

#Unzip file
sapply(fileInp, function(x) utils::unzip(zipfile = x, exdir = DirExtr))


#Read in file
FileInp <- list.files(DirExtr, pattern = paste0(DataType,"_",TimeAgg), full.names = TRUE)
#Name data frame by site
names(FileInp) <- Site
dataFlux <- lapply(FileInp, read.csv)

#Convert to long dataframe with Site variable
dfDataFlux <- purrr::map_df(dataFlux, ~as.data.frame(.x), .id="Site")

#Replace -9999 with NaN
dfDataFlux <- replace(dfDataFlux, dfDataFlux == -9999, NaN)
#Convert data to allow plotting
DateConv <- c("HH" = "ymdHM","DD" = "ymd","WW" = "ymd","MM" = "ym","YY" = "y")
names(dfDataFlux)[which(names(dfDataFlux) == "TIMESTAMP_START")] <- "TIMESTAMP"
dfDataFlux$TIMESTAMP <- lubridate::parse_date_time(dfDataFlux$TIMESTAMP, orders = DateConv[TimeAgg])




#Define missing value fill
mv <- -9999.  
# startStep <- 1

#Set of year/month combinations for netCDF output
setYearMon <- unique(strftime(dataClm$DateTime, "%Y-%m", tz='UTC'))

for (m in setYearMon) {
  #m <- setYearMon[1] #for testing
  Data.mon <- dataClm[dataClm$yearMon == m,]
  timeStep <- seq(0,nrow(Data.mon)-1,1)
  time     <- timeStep/48
  #endStep  <- startStep + nsteps[m]-1
  # not sure why DateTime[1] doesn't include h:m:s
  tempTime <- Data.mon$DateTime[1]
  tempTime <- format(tempTime,'%Y-%m-%d %H:%M:%S')
  
  print(paste(m,"Data date =",tempTime))
  names(Data.mon)
  
  #NetCDF output filenames
  fileOutAtm <- paste(DirOutAtm,"/",Site,"_atm_",m,".nc", sep = "")
  fileOutEval <- paste(DirOutEval,"/",Site,"_eval_",m,".nc", sep = "")
  #sub(pattern = ".txt", replacement = ".nc", fileOut)
  
  DirOut
  
  # define the netcdf coordinate variables (name, units, type)
  lat  <- ncdf4::ncdim_def("lat","degrees_north", as.double(latSite), create_dimvar=TRUE)
  lon <- ncdf4::ncdim_def("lon","degrees_east", as.double(lonSite), create_dimvar=TRUE)
  
  #Character dimension for flags
  #dimnchar <- ncdim_def("nchar", "", 1:17, create_dimvar=FALSE )
  
  #Variables to output to netCDF
  time <- ncdf4::ncdim_def("time", paste("days since",tempTime), calendar = "gregorian",
                           vals=as.double(time),unlim=FALSE, create_dimvar=TRUE )
  LATIXY  <- ncdf4::ncvar_def("LATIXY", "degrees N", list(lat), mv,
                              longname="latitude", prec="double")
  LONGXY  <- ncdf4::ncvar_def("LONGXY", "degrees E", list(lon), mv,
                              longname="longitude", prec="double")
  FLDS  <- ncdf4::ncvar_def("FLDS", "W/m^2", list(lon,lat,time), mv,
                            longname="incident longwave (FLDS)", prec="double")
  FSDS  <- ncdf4::ncvar_def("FSDS", "W/m^2", list(lon,lat,time), mv,
                            longname="incident shortwave (FSDS)", prec="double")
  PRECTmms <- ncdf4::ncvar_def("PRECTmms", "mm/s", list(lon,lat,time), mv,
                               longname="precipitation (PRECTmms)", prec="double")
  PSRF  <- ncdf4::ncvar_def("PSRF", "Pa", list(lon,lat,time), mv,
                            longname="pressure at the lowest atmospheric level (PSRF)", prec="double")
  RH    <- ncdf4::ncvar_def("RH", "%", list(lon,lat,time), mv,
                            longname="relative humidity at lowest atm level (RH)", prec="double")
  TBOT  <- ncdf4::ncvar_def("TBOT", "K", list(lon,lat,time), mv,
                            longname="temperature at lowest atm level (TBOT)", prec="double")
  WIND  <- ncdf4::ncvar_def("WIND", "m/s", list(lon,lat,time), mv,
                            longname="wind at lowest atm level (WIND)", prec="double")
  ZBOT  <- ncdf4::ncvar_def("ZBOT", "m", list(lon,lat,time), mv,
                            longname="observational height", prec="double")
  NEE <- ncdf4::ncvar_def("NEE", "umolm-2s-1", list(lon,lat,time), mv,
                          longname="net ecosystem exchange", prec="double")
  FC <- ncdf4::ncvar_def("FC", "umolm-2s-1", list(lon,lat,time), mv,
                         longname="turbulent CO2 flux", prec="double")
  FSH  <- ncdf4::ncvar_def("FSH", "Wm-2", list(lon,lat,time), mv,
                           longname="sensible heat flux", prec="double")
  EFLX_LH_TOT  <- ncdf4::ncvar_def("EFLX_LH_TOT", "Wm-2", list(lon,lat,time), mv,
                                   longname="latent heat flux", prec="double")
  FSH_NSAE  <- ncdf4::ncvar_def("FSH_NSAE", "Wm-2", list(lon,lat,time), mv,
                                longname="sensible heat flux net surface atmosphere exchange (turbulent + storage)", prec="double")
  EFLX_LH_TOT_NSAE  <- ncdf4::ncvar_def("EFLX_LH_TOT_NSAE", "Wm-2", list(lon,lat,time), mv, longname="latent heat flux net surface atmosphere exchange (turbulent + storage)", prec="double")
  GPP <- ncdf4::ncvar_def("GPP", "umolm-2s-1", list(lon,lat,time), mv,
                          longname="gross primary productivity calculated using nighttime partitioning method following Reichstein et al. (2005)", prec="double")
  GPP_DT <- ncdf4::ncvar_def("GPP_DT", "umolm-2s-1", list(lon,lat,time), mv,
                             longname="gross primary productivity calculated using daytime partitioning method following Lasslop et al. (2010)", prec="double")
  Reco_DT <- ncdf4::ncvar_def("Reco_DT", "umolm-2s-1", list(lon,lat,time), mv,
                              longname="Ecosystem respiration calculated using daytime partitioning method following Lasslop et al. (2010)", prec="double")
  Ustar <- ncdf4::ncvar_def("Ustar", "m/s", list(lon,lat,time), mv,
                            longname="friction velocity", prec="double")
  VPD <- ncdf4::ncvar_def("VPD", "m/s", list(lon,lat,time), mv,
                          longname="vapor pressure deficit", prec="double")
  Rnet  <- ncdf4::ncvar_def("Rnet", "W/m^2", list(lon,lat,time), mv,
                            longname="net radiation", prec="double")
  RadDif <- ncdf4::ncvar_def("RadDif", "W/m^2", list(lon,lat,time), mv,
                             longname="diffuse radiation", prec="double")
  RadDir  <- ncdf4::ncvar_def("RadDir", "W/m^2", list(lon,lat,time), mv,
                              longname="direct radiation", prec="double")
  #quality flag variables
  qfNEE <- ncdf4::ncvar_def("qfNEE", "NA", list(lon,lat,time), mv,
                            longname="net ecosystem exchange final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  qfFC <- ncdf4::ncvar_def("qfFC", "NA", list(lon,lat,time), mv,
                           longname="turbulent CO2 flux final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  qfFSH  <- ncdf4::ncvar_def("qfFSH", "NA", list(lon,lat,time), mv,
                             longname="sensible heat flux final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  qfEFLX_LH_TOT  <- ncdf4::ncvar_def("qfEFLX_LH_TOT", "NA", list(lon,lat,time), mv,
                                     longname="latent heat flux final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  qfFSH_NSAE  <- ncdf4::ncvar_def("qfFSH_NSAE", "NA", list(lon,lat,time), mv,
                                  longname="sensible heat flux net surface atmosphere exchange (turbulent + storage) final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  
  qfEFLX_LH_TOT_NSAE  <- ncdf4::ncvar_def("qfEFLX_LH_TOT_NSAE", "NA", list(lon,lat,time), mv, longname="latent heat flux net surface atmosphere exchange (turbulent + storage) final quality flag, 0 = good data and 1 = flagged data", prec="integer")
  
  #gap-filling quality flag variables
  FLDS_fqc  <- ncdf4::ncvar_def("FLDS_fqc", "NA", list(lon,lat,time),
                                longname="incident longwave (FLDS) gap-filling flag", prec="integer")
  FSDS_fqc  <- ncdf4::ncvar_def("FSDS_fqc", "NA", list(lon,lat,time),
                                longname="incident shortwave (FSDS) gap-filling flag", prec="integer")
  PRECTmms_fqc <- ncdf4::ncvar_def("PRECTmms_fqc", "NA", list(lon,lat,time), 
                                   longname="precipitation (PRECTmms) gap-filling flag", prec="integer")
  PSRF_fqc  <- ncdf4::ncvar_def("PSRF_fqc", "NA", list(lon,lat,time),
                                longname="pressure at the lowest atmospheric level (PSRF) gap-filling flag", prec="integer")
  RH_fqc    <- ncdf4::ncvar_def("RH_fqc", "NA", list(lon,lat,time), 
                                longname="relative humidity at lowest atm level (RH) gap-filling flag", prec="integer")
  TBOT_fqc  <- ncdf4::ncvar_def("TBOT_fqc", "NA", list(lon,lat,time),
                                longname="temperature at lowest atm level (TBOT) gap-filling flag", prec="integer")
  WIND_fqc  <- ncdf4::ncvar_def("WIND_fqc", "NA", list(lon,lat,time), 
                                longname="wind at lowest atm level (WIND) gap-filling flag", prec="integer")
  NEE_fqc <- ncdf4::ncvar_def("NEE_fqc", "NA", list(lon,lat,time), 
                              longname="net ecosystem exchange gap-filling flag", prec="integer")
  FSH_fqc  <- ncdf4::ncvar_def("FSH_fqc", "NA", list(lon,lat,time),
                               longname="sensible heat flux gap-filling flag", prec="integer")
  EFLX_LH_TOT_fqc  <- ncdf4::ncvar_def("EFLX_LH_TOT_fqc", "NA", list(lon,lat,time),
                                       longname="latent heat flux gap-filling flag", prec="integer")
  GPP_fqc <- ncdf4::ncvar_def("GPP_fqc", "NA", list(lon,lat,time),
                              longname="gross primary productivity gap-filling flag", prec="integer")
  Ustar_fqc <- ncdf4::ncvar_def("Ustar_fqc", "NA", list(lon,lat,time),
                                longname="friction velocity gap-filling flag", prec="integer")
  VPD_fqc <- ncdf4::ncvar_def("VPD_fqc", "NA", list(lon,lat,time),
                              longname="vapor pressure deficit gap-filling flag", prec="integer")
  Rnet_fqc  <- ncdf4::ncvar_def("Rnet_fqc", "NA", list(lon,lat,time),
                                longname="net radiation gap-filling flag", prec="integer")
  
  #Create the output file
  ncAtm <- ncdf4::nc_create(fileOutAtm, list(LATIXY,LONGXY,FLDS,FSDS,PRECTmms,RH,PSRF,TBOT,WIND,ZBOT,FLDS_fqc,FSDS_fqc,PRECTmms_fqc,RH_fqc,PSRF_fqc,TBOT_fqc,WIND_fqc))
  
  ncEval <- ncdf4::nc_create(fileOutEval, list(LATIXY,LONGXY,NEE,FC,FSH,EFLX_LH_TOT,FSH_NSAE,EFLX_LH_TOT_NSAE,GPP,GPP_DT,Reco_DT,Ustar,VPD,Rnet,RadDif,RadDir,ZBOT,NEE_fqc,FSH_fqc,EFLX_LH_TOT_fqc,GPP_fqc,Ustar_fqc,VPD_fqc,Rnet_fqc,qfNEE,qfFC,qfFSH,qfEFLX_LH_TOT,qfFSH_NSAE,qfEFLX_LH_TOT_NSAE))
  
  
  # Write some values to this variable on disk.
  ncdf4::ncvar_put(ncAtm, LATIXY, latSite)
  ncdf4::ncvar_put(ncAtm, LONGXY, lonSite)
  ncdf4::ncvar_put(ncAtm, FLDS, Data.mon$FLDS)
  ncdf4::ncvar_put(ncAtm, FSDS, Data.mon$FSDS)
  ncdf4::ncvar_put(ncAtm, RH,   Data.mon$RH)
  ncdf4::ncvar_put(ncAtm, PRECTmms, Data.mon$PRECTmms)
  ncdf4::ncvar_put(ncAtm, PSRF, Data.mon$PSRF)
  ncdf4::ncvar_put(ncAtm, TBOT, Data.mon$TBOT)
  ncdf4::ncvar_put(ncAtm, WIND, Data.mon$WIND)
  ncdf4::ncvar_put(ncAtm, ZBOT, Data.mon$ZBOT)
  ncdf4::ncvar_put(ncAtm, FLDS_fqc, Data.mon$FLDS_fqc)
  ncdf4::ncvar_put(ncAtm, FSDS_fqc, Data.mon$FSDS_fqc)
  ncdf4::ncvar_put(ncAtm, RH_fqc,   Data.mon$RH_fqc)
  ncdf4::ncvar_put(ncAtm, PRECTmms_fqc, Data.mon$PRECTmms_fqc)
  ncdf4::ncvar_put(ncAtm, PSRF_fqc, Data.mon$PSRF_fqc)
  ncdf4::ncvar_put(ncAtm, TBOT_fqc, Data.mon$TBOT_fqc)
  ncdf4::ncvar_put(ncAtm, WIND_fqc, Data.mon$WIND_fqc)
  
  #add attributes
  #ncdf4::ncatt_put(ncAtm, time,"calendar", "gregorian" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, FLDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, FSDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, RH  ,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, PRECTmms,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, PSRF,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, TBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, WIND,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, ZBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, FLDS_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, FSDS_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, RH_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, PRECTmms_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, PSRF_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, TBOT_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, WIND_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, FLDS_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, FSDS_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, RH_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, PRECTmms_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, PSRF_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, TBOT_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncAtm, WIND_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  
  ncdf4::ncatt_put(ncAtm, 0, "created_on",date()      ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "created_by",user,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "created_from",fileOut   ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "NEON site",Site         ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "TimeDiffUtcLt",metaSite$TimeDiffUtcLst          ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "created_with", "flow.api.clm.R",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncAtm, 0, "supported_by", "This data development was funded by the National Science Foundation (NSF) solicitation PD 20-7684",prec=NA,verbose=FALSE,definemode=FALSE )
  
  
  # Write some values to this variable on disk.
  ncdf4::ncvar_put(ncEval, LATIXY, latSite)
  ncdf4::ncvar_put(ncEval, LONGXY, lonSite) 
  ncdf4::ncvar_put(ncEval, NEE, Data.mon$NEE)
  ncdf4::ncvar_put(ncEval, FC, Data.mon$FC)
  ncdf4::ncvar_put(ncEval, FSH, Data.mon$H)
  ncdf4::ncvar_put(ncEval, FSH_NSAE, Data.mon$H_NSAE)
  ncdf4::ncvar_put(ncEval, EFLX_LH_TOT, Data.mon$LE)
  ncdf4::ncvar_put(ncEval, EFLX_LH_TOT_NSAE, Data.mon$LE_NSAE)
  ncdf4::ncvar_put(ncEval, GPP, Data.mon$GPP)
  ncdf4::ncvar_put(ncEval, GPP_DT, Data.mon$GPP_DT)
  ncdf4::ncvar_put(ncEval, Reco_DT, Data.mon$Reco_DT)
  ncdf4::ncvar_put(ncEval, Ustar, Data.mon$Ustar)
  ncdf4::ncvar_put(ncEval, VPD, Data.mon$VPD)
  ncdf4::ncvar_put(ncEval, Rnet, Data.mon$radNet)
  ncdf4::ncvar_put(ncEval, RadDir, Data.mon$RadDir)
  ncdf4::ncvar_put(ncEval, RadDif, Data.mon$RadDif)
  ncdf4::ncvar_put(ncEval, ZBOT, Data.mon$ZBOT)
  ncdf4::ncvar_put(ncEval, NEE_fqc, Data.mon$NEE_fqc)
  ncdf4::ncvar_put(ncEval, FSH_fqc, Data.mon$H_fqc)
  ncdf4::ncvar_put(ncEval, EFLX_LH_TOT_fqc, Data.mon$LE_fqc)
  ncdf4::ncvar_put(ncEval, GPP_fqc, Data.mon$GPP_fqc)
  ncdf4::ncvar_put(ncEval, Ustar_fqc, Data.mon$Ustar_fqc)
  ncdf4::ncvar_put(ncEval, VPD_fqc, Data.mon$VPD_fqc)
  ncdf4::ncvar_put(ncEval, Rnet_fqc, Data.mon$radNet_fqc)
  ncdf4::ncvar_put(ncEval, qfNEE, Data.mon$qfNEE)
  ncdf4::ncvar_put(ncEval, qfFC, Data.mon$qfFC)
  ncdf4::ncvar_put(ncEval, qfFSH, Data.mon$qfH)
  ncdf4::ncvar_put(ncEval, qfFSH_NSAE, Data.mon$qfH_NSAE)
  ncdf4::ncvar_put(ncEval, qfEFLX_LH_TOT, Data.mon$qfLE)
  ncdf4::ncvar_put(ncEval, qfEFLX_LH_TOT_NSAE, Data.mon$qfLE_NSAE)
  
  
  ncdf4::ncatt_put(ncEval, NEE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, FC,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, FSH,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, FSH_NSAE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, EFLX_LH_TOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, EFLX_LH_TOT_NSAE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, GPP,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, GPP_DT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, Reco_DT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, Ustar,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, VPD,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, Rnet,"mode","time-dependent",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, RadDif,"mode","time-dependent",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, RadDir,"mode","time-dependent",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, ZBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, NEE_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, FSH_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, EFLX_LH_TOT_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, GPP_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, Ustar_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, VPD_fqc,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, Rnet_fqc,"mode","time-dependent",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, NEE_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, FSH_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, EFLX_LH_TOT_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, GPP_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, Ustar_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, VPD_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, Rnet_fqc,"method_gap-fill","0=no gap-filling, 1=regression, 2=ReddyProc_methA, 3=ReddyProc_methB, 4=ReddyProc_methC" ,prec=NA,verbose=FALSE,definemode=FALSE)
  
  #Quality flag attributes
  ncdf4::ncatt_put(ncEval, qfNEE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, qfFC,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, qfFSH,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, qfFSH_NSAE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, qfEFLX_LH_TOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, qfEFLX_LH_TOT_NSAE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  
  #Global attributes
  ncdf4::ncatt_put(ncEval, 0, "created_on",date()      ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, 0, "created_by",user,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, 0, "created_from",fileOut   ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, 0, "NEON site",Site         ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, 0, "TimeDiffUtcLt",metaSite$TimeDiffUtcLst         ,prec=NA,verbose=FALSE,definemode=FALSE)
  ncdf4::ncatt_put(ncEval, 0, "created_with", "flow.api.clm.R",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncEval, 0, "supported_by", "This data development was funded by the National Science Foundation (NSF) solicitation PD 20-7684",prec=NA,verbose=FALSE,definemode=FALSE )
  
  #Close Netcdf file connection
  ncdf4::nc_close(ncAtm)
  ncdf4::nc_close(ncEval)
  
  remove(time, timeStep, fileOutAtm, fileOutEval, ncAtm, ncEval, Data.mon,
         FLDS,FSDS,RH,PRECTmms,PSRF,TBOT,WIND,ZBOT, NEE,FC, FSH, FSH_NSAE, EFLX_LH_TOT, EFLX_LH_TOT_NSAE,GPP,GPP_DT,Reco_DT, Ustar,VPD, Rnet,RadDif,RadDir,
         FLDS_fqc,FSDS_fqc,RH_fqc,PRECTmms_fqc,PSRF_fqc,TBOT_fqc,WIND_fqc, NEE_fqc, FSH_fqc,EFLX_LH_TOT_fqc, GPP_fqc, Ustar_fqc,VPD_fqc, Rnet_fqc,qfNEE,qfFC,qfFSH,qfFSH_NSAE,qfEFLX_LH_TOT,qfEFLX_LH_TOT_NSAE)
} #End of monthloop
