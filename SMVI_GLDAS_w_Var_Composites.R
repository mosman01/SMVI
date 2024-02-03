### ### ### ### ### ### ### ### ### 
#######FLASH DROUGHTS - SMVI#######
### ### ### ### ### ### ### ### ###
##########NLDAS/GLDAS based##############
#####LAST UPDATED 02/01/2023#######
######################################################
##Author: Mahmoud Osman, PhD##########################
##Johns Hopkins University  ##########################
##Email: mahmoud.osman@jhu.edu; mahosman01@gmail.com##
######################################################

######!#!#!####################!#!##!##!##!##!##!##!###########################!#!#!######
##!!!!Modification to this script is only permitted by the author's permission!!!!!#######
##!!!!It is not permitted to distribute the script without the author's permission!!!!####
######!#!#!####################!#!##!##!##!##!##!##!###########################!#!#!######

##References:
###Osman, M., Zaitchik, B. F., Badr, H. S., Christian, J. I., Tadesse, T., Otkin, J. A., & Anderson, M. C. (2021). Flash drought onset over the contiguous United States: sensitivity of inventories and trends to quantitative definitions. Hydrology and Earth System Sciences, 25(2), 565–581. https://doi.org/10.5194/hess-25-565-2021
###Osman, M., Zaitchik, B. F., Badr, H. S., Otkin, J., Zhong, Y., Lorenz, D., et al. (2022). Diagnostic Classification of Flash Drought Events Reveals Distinct Classes of Forcings and Impacts. Journal of Hydrometeorology, 23(2), 275–289. https://doi.org/10.1175/JHM-D-21-0134.1
###Osman, M., B. Zaitchik, J. Otkin, M. Anderson (2024). SMVI Global Flash Droughts Dataset, HydroShare, https://doi.org/10.4211/hs.080002bd7cc44242bb37c02b049ed532
##This research has been supported by the National Science Foundation (grant no. 1854902)


###OUTPUT will be saved at the directory: SMVI_<Dataset_Name>_out within the defined working directory. Refer to the README file for more details.

###############How to pre-process input files################
###Note! CDO tool can do all the needed pre-processing#####
#Get the Meteorological data in daily format, if it's hourly (just like GLDAS_NOAH025) then daily averaging is required.
##For example after downloading all hourly data for 1 month:
####cdo mergetime infiles outfile1
##Create daily data:
####cdo -daymean outfile1 outfile2
##Create Climatology/Stdv/Min/Max from all processed data of all years (even if case study is 1 month, climatology is based on long-term mean):
####cdo -ydaymean outfile2 outfile31
####cdo -ydaystd outfile2 outfile32
####cdo -ydaymin outfile2 outfile33
####cdo -ydaymax outfile2 outfile34
##Create Anomalies and Std-Anomalies:
####cdo -ydaysub outfile2 outfile31 outfile4
####cdo -ydaydiv outfile4 outfile5      #(This is the standardized anomaly)
##Create Percentiles files (e.g. 20th pct):
####cdo -ydaypctl,20 outfile2 outfile33 outfile34 outfile6
##To select a specific variable:
####cdo -selvar,RZSM infile outfile
##To select a specific region:
####cdo -sellonlatbox,-125.0,-67.0,25.0,53.0 infile outfile
##To select a specific region:
####cdo -sellonlatbox,-125.0,-67.0,25.0,53.0 infile outfile
#
#No need to put VAR used in identifying FD using SMVI method in std-anomaly since the method is relative to the change. i.e. keep percentiles and daily data in absolute values.
#However, for composites, std-anomalies are required!
#
##Compute windspeed from U and V:
####cdo -P 12 -v -z zip_9 -expr,'WS=sqrt(UGRD*UGRD+VGRD*VGRD)' NLDAS_FORA0125_D.1979_2021_UV.nc NLDAS_FORA0125_D.A21979_2021_WS.nc
##OR::
####cdo -z zip_5 chname,u10,ws -sqrt -add -sqr -selname,u10 "$input_u_file" -sqr -selname,v10 "$input_v_file" "$output_file"
##To add timeaxis if missing:
####cdo -v -z zip_4 -setreftime,1990-01-01,00:00:00,1day -settaxis,1990-01-01,00:00:00,1day -setcalendar,standard ERA5_1990_ALL.nc ERA5_1990_ALL_TEST.nc

#Note some flags to compress and parallel process

#####

################Beginning of script###################


##Clean up R environment and memory 
#.rs.restartR()
rm(list = ls())
if(is.null(dev.list())==F){
  dev.off()  
}
gc()
ptm <- proc.time()  ##To compute time of script

#######LIBRARIES##########################
#Check  & Load libraries:
##Some package may be extra. Will clean up the extra ones in the next update
#list.of.packages <- c("maps","maptools","ncdf4","zoo","data.table","xts","pracma","TTR","RColorBrewer","abind","Cairo","birk",'latticeExtra','lattice',"foreach","doParallel","doSNOW")
list.of.packages <- c("maps","maptools","ncdf4","zoo","data.table","xts","pracma","TTR","RColorBrewer","abind","Cairo","birk",'latticeExtra','lattice')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,dependencies = T)
lapply(list.of.packages, library, character.only = TRUE)
#update.packages(checkBuilt=T,ask=F)

#rasterTmpFile('mytemp_')
#showTmpFiles()
#removeTmpFiles(h=0.1)

##########################################


#####USER INPUT#####
############################

#Dataset name
DS <- 'GLDAS'  #(NLDAS or GLDAS) -- Edit files names in lines 190-243 to match yours.
#Set to parent Working Dir:
#setwd('/mnt/sahara/data1/mahmoud/FD/')
setwd('/mnt/local_drive/Mahmoud/FD/GLOBAL_FD/')


#cores=detectCores()
#cl <- makeCluster(cores[1]/3) #not to overload your computer
#registerDoParallel(cl)

#Set Output to screen or PNG
##Only if plotting - Currently commented out
out_fig<- T #(T for PNG or F for no figure)

transparent_BG<- F  #Transparent background

BR_msk<- T #filter out regions by Bowen's ratio
BR_min<- 0.2
BR_max<- 7

TS_msk<- T #filter out regions by Surface Temperature - Usefull to avoid on-snow detected FD
TS_thd<- 273	#min surface temperature in K to detect FD

#YYYY<-2012
YEAR_BEGIN<- 2021
YEAR_END<- 2021

VAR<-"RZSM"		#(Analysis Variable. RZSM is the default)

COMPOSITES <- T #Compute composites (T), Skip composites (F). Composites calculation is a lengthy process, and requires all variables data.
#Anom_dir <- "/mnt/local_drive/Mahmoud/ERA5/Merged/"
Anom_dir <- "ERA5"
NCompAv <- 1  	#change for a different running averageing window of variables compsites - default=1

#Select analysis method (! each method has specific options)
#(R: Bukovsky regions, S: Raster shapefile, C: Coordinates, ST: US states, A: All gridded points)
RSC<-'A'	 ##This script works only for 'A'

#If Coordinates are used fill in coordinates
#If shapefile is used: make sure it's WGS84 and saved as .img raster format inside Shape_Rstr directory
#(Do that in ArcMap "Polygon to Raster' after creating the shapefile)


#####FD Def. Method########

#Threshold: (20, 30, 40) 'th' percentiles or Min, Mean (Climatology based) or 20days runnung average
#Time Varying Variable: 5days runnung average or Original input variable

##Predefined methods:
# "SMVI" = 5 days less than 20 days and below 20th pctl in 4 pen at least. (Osman et al. 2021)
FD_mthd<- "SMVI" #DO NOT CHANGE THIS!!
# Changing following lines means modifying/tuning SMVI method, which is acceptable based on the understanding of the data and conditions.
SM_Th<- 20        #Soil moisture threshold (20th Pctl is recommended based on USDM D1 & Otkin et al., use 'MAX' to disable this)
FD_Th <- '20DR'   # choose the threshold from (20,30,40,50,60,90,'MIN','MEAN','20DR') - 20DR means 20 days running average, 20 means 20th percentile
FD_V <- 5         # choose (1 for original variable or 5 for 5days running average)
pentads<- 4       #min length of FD event
FDL<- pentads*5   #min length of FD event, just multiplying
FDgap <- 15       #Neglected gap between FD events

nlon=1440	#Use 464 for NLDAS, use 1440 for GLDAS
nlat=600	#Use 224 for NLDAS, use 600 for GLDAS
#ntime=275	#daily from 01Mar-30Nov. Make sure to change this if different period is used
DAY_FIRST= "01/01"	#first day of the caclulation season
DAY_LAST= "12/31"	#last day of the calculation season 

nc01x <- nc_open(paste('GLDAS/GLDAS_CLSM025_1980_2022_SEL_CLIM_RZSM.nc',sep = '')) #use a small filre from your dataset to grab lon and lat from it

##Dates are MAR 1 to NOV 30 ------ That may require changing to fit the calculation season
##Consequently the timesries length will need recalculation
#time_d<-seq(as.Date(paste(YYYY,"/01/01",sep='')),as.Date(paste(YYYY,"/12/31",sep='')), by = "1 day")

###########################################
## ## ## END OF USER INPUT## ## ## ##

lon <- ncvar_get(nc01x,var="lon")
lon[lon > 180] <- lon[lon > 180] - 360
lat <- ncvar_get(nc01x,var="lat")

#LONLAT= expand.grid(LON = lon, LAT = lat)
#write.csv(LONLAT  ,file = paste('LONLAT',".csv",sep=''))

rm(nc01x)

for (YYYY in YEAR_BEGIN:YEAR_END){

time_d<-seq(as.Date(paste(YYYY,"/01/01",sep='')),as.Date(paste(YYYY,"/12/31",sep='')), by = "1 day")


if (COMPOSITES == T) {

	print("READING COMPOSITE FILES...")
	#This section can be highly improved if all variables are stored in 1 file, then only first "in_var" is used and delete all similar lines
	
	#Composite data reading:
	#VAR= PEVPR, TMP, ARAIN, PRESS, EVP, RZSM, WS, VPD
	###ERA5:  pev, t2m, tp, sp, slhf, RZSM, ws, vpd	#RZSM is still taken from GLDAS
	PET_var = 'pev'
	TMP_var = 't2m'
	PRCP_var = 'tp'
	SP_var = 'sp'
	ET_var = 'slhf'
	RZSM_var = 'RZSM'
	WS_var = 'ws'
	VPD_var = 'vpd'
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir, '/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),TMP_var)
	ntime <- dim(in_var)[3]
	TMP_arr<-array(in_var,c((nlon*nlat),ntime)) #Temperature

	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),PET_var)
	PEVPR_arr<-array(in_var,c((nlon*nlat),ntime)) #Potential Evapotranspiration
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),PRCP_var)
	ARAIN_arr<-array(in_var,c((nlon*nlat),ntime)) #Precipitation
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),SP_var)
	PRES_arr<-array(in_var,c((nlon*nlat),ntime)) #Surface Pressure
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),ET_var)
	EVP_arr<-array(in_var,c((nlon*nlat),ntime)) #Evapotranspiration
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),WS_var)
	WS_arr<-array(in_var,c((nlon*nlat),ntime)) #Wind Speed
	
	in_var<- ncvar_get(nc_open(paste(Anom_dir,'/Years/ERA_ANOM_STD_',YYYY,'.nc',sep='')),VPD_var)
	VPD_arr<-array(in_var,c((nlon*nlat),ntime)) #Vapor Pressure Deficit
}

rm(in_var)

##nrow is equal to nlon*nlat
#var_arr<-matrix(in_var,nrow = 103936, ncol = 275, byrow = F)


#######ALL POINTS######
print('READING RZSM FILES...')

#nc01 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_',YYYY,'.nc',sep = '')),'RZSM')
nc01<- ncvar_get(nc_open(paste('GLDAS/YEARS/GLDAS_CLSM025_ANOM_STD_',YYYY,'.nc',sep='')),'SoilMoist_RZ_tavg') #GLDAS
n_reg<-dim(nc01)[1]*dim(nc01)[2]   #=103936
ntime <- dim(nc01)[3]

nc01.2 <- array(nc01,c((n_reg),ntime))
rm(nc01)
##Could be improved by saveing all percentiles in 1 netcdf file

print("READING VEG FILE...")

	if (ntime==365) {
		in_var<- ncvar_get(nc_open(paste('GLDAS/GLDAS_CLSM025_dveg_365.nc',sep='')),'dveg')
	} else {
		in_var<- ncvar_get(nc_open(paste('GLDAS/GLDAS_CLSM025_dveg_366.nc',sep='')),'dveg')
	}
	
	VEG_arr<-array(in_var,c((nlon*nlat),ntime)) #Vegetation
	

if (SM_Th != 20){ #need to be used with caution!!

	nc02 <- ncvar_get(nc_open(paste('GLDAS/GLDAS_CLSM025_1980_2022_SEL_20PCTL.nc',sep = '')),'SoilMoist_RZ_tavg')
	nc03 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_1979_2018_10PCTL.nc',sep = '')),'RZSM')
	nc04 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_1979_2018_90PCTL.nc',sep = '')),'RZSM')
	nc05 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_1979_2018_40PCTL.nc',sep = '')),'RZSM')
	#nc02 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_1979_2018_20PCTL.nc',sep = '')),'RZSM')
	
	nc02.2 <- array(nc02,c((n_reg),ntime))
	nc03.2 <- array(nc03,c((n_reg),ntime))
	nc04.2 <- array(nc04,c((n_reg),ntime))
	nc05.2 <- array(nc05,c((n_reg),ntime))
	
	rm(nc02, nc03, nc04, nc05)

} else {
	print("READING PCTL FILE...")

	nc02 <- ncvar_get(nc_open(paste('GLDAS/GLDAS_CLSM025_1980_2022_SEL_20PCTL.nc',sep = '')),'SoilMoist_RZ_tavg')
	#nc02 <- ncvar_get(nc_open(paste('NLDAS/RZSM_NLDAS_NOAH0125_D_1979_2018_20PCTL.nc',sep = '')),'RZSM')
	nc02.2 <- array(nc02,c((n_reg),ntime))
	rm(nc02)
}


if (BR_msk == T){ #Bowen's ratio condition
	print("READING BR FILE...")
	ncbr<- ncvar_get(nc_open(paste('GLDAS/YEARS/GLDAS_CLSM025_BR_',YYYY,'.nc',sep='')),'BR') #GLDAS
	ncbr.2 <- array(ncbr,c((n_reg),ntime))
	rm(ncbr)
}
 if (TS_msk == T){ #Surface temperature condition
	print("READING T_SURF FILE...")
	ncts<- ncvar_get(nc_open(paste('GLDAS/YEARS/GLDAS_CLSM025_ST_',YYYY,'.nc',sep='')),'AvgSurfT_tavg') #GLDAS
	ncts.2 <- array(ncts,c((n_reg),ntime))
	rm(ncts)
}  

###### FULL MAP######
#Here we create empty data-frames to store the output     
print("CREATING DATAFRAMES")   
    dm<-rep(NA,n_reg)
    dm.df <- data.frame(L_MX1=as.double(rep(NA,n_reg)),
                        L_MX2=as.double(rep(NA,n_reg)),
                        L_MX3=as.double(rep(NA,n_reg)),
                        L_MX4=as.double(rep(NA,n_reg)),
                        L_MX5=as.double(rep(NA,n_reg)),
                        L_MX6=as.double(rep(NA,n_reg)),
                        fstdate1=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate1=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        fstdate2=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate2=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        fstdate3=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate3=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        fstdate4=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate4=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        fstdate5=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate5=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        fstdate6=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        lstdate6=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                        FD_Max=as.double(rep(NA,n_reg)))
    
	dm.df.2 <- data.frame(fstdate1=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate1=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          fstdate2=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate2=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          fstdate3=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate3=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          fstdate4=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate4=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          fstdate5=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate5=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          fstdate6=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdate6=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          FD_Max=as.double(rep(NA,n_reg)),
						  fstdateMx=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
                          lstdateMx=as.Date(rep('1900-01-01',n_reg),'%y-%m-%d'),
						  SV1=as.double(rep(NA,n_reg)),
						  SV2=as.double(rep(NA,n_reg)),
						  SV3=as.double(rep(NA,n_reg)),
						  SV4=as.double(rep(NA,n_reg)),
						  SV5=as.double(rep(NA,n_reg)),
						  SV6=as.double(rep(NA,n_reg)),
						  SVMx=as.double(rep(NA,n_reg)),
						  VEGID=as.double(rep(NA,n_reg)))
						  
	if (COMPOSITES == T) { 
	
	print("COMPOSITES DATAFRAMES")

		   DF_X <- data.frame(On17=as.double(rep(NA,n_reg)),
								On16=as.double(rep(NA,n_reg)),
								On15=as.double(rep(NA,n_reg)),
								On14=as.double(rep(NA,n_reg)),
								On13=as.double(rep(NA,n_reg)),
								On12=as.double(rep(NA,n_reg)),
								On11=as.double(rep(NA,n_reg)),
								On10=as.double(rep(NA,n_reg)),
								On09=as.double(rep(NA,n_reg)),
								On08=as.double(rep(NA,n_reg)),
								On07=as.double(rep(NA,n_reg)),
								On06=as.double(rep(NA,n_reg)),
								On05=as.double(rep(NA,n_reg)),
								On04=as.double(rep(NA,n_reg)),
								On03=as.double(rep(NA,n_reg)),
								On02=as.double(rep(NA,n_reg)),
								On01=as.double(rep(NA,n_reg)),
								On00=as.double(rep(NA,n_reg)),
								Re00=as.double(rep(NA,n_reg)),
								Re01=as.double(rep(NA,n_reg)),
								Re02=as.double(rep(NA,n_reg)),
								Re03=as.double(rep(NA,n_reg)),
								Re04=as.double(rep(NA,n_reg)),
								Re05=as.double(rep(NA,n_reg)),
								Re06=as.double(rep(NA,n_reg)),
								Re07=as.double(rep(NA,n_reg)),
								Re08=as.double(rep(NA,n_reg)),
								Re09=as.double(rep(NA,n_reg)),
								Re10=as.double(rep(NA,n_reg)),
								Re11=as.double(rep(NA,n_reg)),
								Re12=as.double(rep(NA,n_reg)),
								Re13=as.double(rep(NA,n_reg)),
								Re14=as.double(rep(NA,n_reg)),
								Re15=as.double(rep(NA,n_reg)),
								Re16=as.double(rep(NA,n_reg)),
								Re17=as.double(rep(NA,n_reg)))

		#Note that there are 6 data-frames for each variables as we count for up to 6 events of FD per year (it never happens anyway!) sorted by longest.
		DF_X_TMP  <- DF_X_PEVPR<- DF_X_ARAIN<- DF_X_PRES <- DF_X_RZSM <- DF_X_EVP  <- DF_X_WS   <- DF_X_VPD  <- DF_X

			DF_X6_TMP  <- DF_X5_TMP  <- DF_X4_TMP  <- DF_X3_TMP  <- DF_X2_TMP   <- DF_X_TMP  
			DF_X6_PEVPR<- DF_X5_PEVPR<- DF_X4_PEVPR<- DF_X3_PEVPR<- DF_X2_PEVPR <- DF_X_PEVPR
			DF_X6_ARAIN<- DF_X5_ARAIN<- DF_X4_ARAIN<- DF_X3_ARAIN<- DF_X2_ARAIN <- DF_X_ARAIN
			DF_X6_PRES <- DF_X5_PRES <- DF_X4_PRES <- DF_X3_PRES <- DF_X2_PRES  <- DF_X_PRES 
			DF_X6_RZSM <- DF_X5_RZSM <- DF_X4_RZSM <- DF_X3_RZSM <- DF_X2_RZSM  <- DF_X_RZSM 
			DF_X6_EVP  <- DF_X5_EVP  <- DF_X4_EVP  <- DF_X3_EVP  <- DF_X2_EVP   <- DF_X_EVP  
			DF_X6_WS   <- DF_X5_WS   <- DF_X4_WS   <- DF_X3_WS   <- DF_X2_WS    <- DF_X_WS   
			DF_X6_VPD  <- DF_X5_VPD  <- DF_X4_VPD  <- DF_X3_VPD  <- DF_X2_VPD   <- DF_X_VPD  
			
			
			DF_Y <- data.frame(On01=as.double(rep(NA,n_reg)),
								On02=as.double(rep(NA,n_reg)),
								On03=as.double(rep(NA,n_reg)),
								On00=as.double(rep(NA,n_reg)),
								Re00=as.double(rep(NA,n_reg)),
								Re01=as.double(rep(NA,n_reg)),
								Re02=as.double(rep(NA,n_reg)),
								Re03=as.double(rep(NA,n_reg)),
								On3Pn=as.double(rep(NA,n_reg)))
			 
			DF_Y_TMP <- DF_Y_PEVPR <- DF_Y_ARAIN <- DF_Y_PRES <- DF_Y_RZSM <- DF_Y_EVP <- DF_Y_WS <- DF_Y_VPD <- DF_Y
			
			DF_Y6_TMP  <- DF_Y5_TMP  <- DF_Y4_TMP   <- DF_Y3_TMP   <- DF_Y2_TMP   <- DF_Y_TMP  
			DF_Y6_PEVPR<- DF_Y5_PEVPR<- DF_Y4_PEVPR <- DF_Y3_PEVPR <- DF_Y2_PEVPR <- DF_Y_PEVPR
			DF_Y6_ARAIN<- DF_Y5_ARAIN<- DF_Y4_ARAIN <- DF_Y3_ARAIN <- DF_Y2_ARAIN <- DF_Y_ARAIN
			DF_Y6_PRES <- DF_Y5_PRES <- DF_Y4_PRES  <- DF_Y3_PRES  <- DF_Y2_PRES  <- DF_Y_PRES 
			DF_Y6_RZSM <- DF_Y5_RZSM <- DF_Y4_RZSM  <- DF_Y3_RZSM  <- DF_Y2_RZSM  <- DF_Y_RZSM 
			DF_Y6_EVP  <- DF_Y5_EVP  <- DF_Y4_EVP   <- DF_Y3_EVP   <- DF_Y2_EVP   <- DF_Y_EVP  
			DF_Y6_WS   <- DF_Y5_WS   <- DF_Y4_WS    <- DF_Y3_WS    <- DF_Y2_WS    <- DF_Y_WS   
			DF_Y6_VPD  <- DF_Y5_VPD  <- DF_Y4_VPD   <- DF_Y3_VPD   <- DF_Y2_VPD   <- DF_Y_VPD  
  	} 
    
	print("STARTING FOR LOOP")

	for (ic in 1:n_reg){
    #foreach (ic=1:n_reg, .packages=list.of.packages) %dopar% {
	#sink("log.txt", append=TRUE)
	#cat(paste("Starting iteration",ic,"\n"))
	#cat(paste("\n","grid# ",ic," ",ic/n_reg*100,"%","\n"),file=paste("log.SMVI.",YYYY,".txt",sep=""), append=F)
	   
	print(paste(ic,"-",round(ic/n_reg*100,2),"%"))

	  msk_Day <- xts(nc01.2[ic,],time_d)
	  if(all(is.na(msk_Day))){
	  } else {
		
        #progress(ic)
        #Sys.sleep(0.01)
        if (ic == n_reg) { cat("Done!\n") }
		
			msk_5Day<-SMA(msk_Day,n = 5)
			msk_20Day<-SMA(msk_Day,n = 20)
			
			if (SM_Th != 20){
				msk_10th<-xts(nc03.2[ic,],time_d)
				msk_40th<-xts(nc05.2[ic,],time_d)
				msk_90th<-xts(nc04.2[ic,],time_d)
			} else {
				msk_20th<-xts(nc02.2[ic,],time_d)
			}
			
			if (BR_msk == T){
			msk_BR<-xts(ncbr.2[ic,],time_d)
			}

			if (TS_msk == T){
			msk_TS<-xts(ncts.2[ic,],time_d)
			}

      ##########FM Find FD!!#########
      #e.g. If below 20th pctl, choose msk_20pctl as the threshold value for Drought (THD) and Comparing value (Cvar) is by default set to daily data
      print("Finding FDs")
      
      if (FD_Th==10) { THD<-msk_10th 
		} else if (FD_Th==20) { THD<-msk_20th 
			} else if (FD_Th==30) { THD<-msk_30th
				} else if (FD_Th==40) { THD<-msk_40th
					} else if (FD_Th==50) { print("Choose a correct First SM threshold from list")
						} else if (FD_Th==60) { print("Choose a correct First SM threshold from list") 
							} else if (FD_Th==90) { THD<-msk_90th
								} else if (FD_Th=='MIN') { THD<-msk_MIN
									} else if (FD_Th=="MEAN") { THD<-msk_MEAN
										} else if (FD_Th=="20DR") { THD<-msk_20Day
											} else if (FD_Th=='1_50') { THD<-msk_Day50
												} else { print("Choose a correct First SM threshold from list") }
         
      #Setting a second threshold based on SM percentile
      if (SM_Th==10) { THD2<-msk_10th
		} else if (SM_Th==20) { THD2<-msk_20th
			} else if (SM_Th==30) { THD2<-msk_30th
				} else if (SM_Th==40) { THD2<-msk_40th 
					} else if (SM_Th==50) { print("Choose a correct First SM threshold from list") 
						} else if (SM_Th==60) { print("Choose a correct First SM threshold from list")
							} else if (SM_Th==90) { THD2<-msk_90th
								} else if (SM_Th=='MAX') { THD2<-msk_MAX
									} else if (SM_Th=="MEAN") { THD2<-msk_MEAN 
										} else if (SM_Th=="20DR") { THD2<-msk_20Day
											} else if (SM_Th=='1_50') { THD2<-msk_Day50
												} else {print("Choose a correct second SM threshold from list") }
  
      
      if (FD_V==1) {cvar<-msk_Day
		} else if (FD_V==5) { cvar<-msk_5Day
			} else { print("Choose Either 1 or 5 in FD_V")}

      
      #cvar<-msk_Day #Change to other comparing variable if needed
	  #Uncomment for testing only
      #THD<-msk_30th  #change when testing
      
      #Finding Event begining and end
      THD_diff<-cvar-THD
    
      fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
      lstNeg_dt<-index(THD_diff[THD_diff<0])[length(THD_diff[THD_diff<0])]
      
      x<-1
      if(all(is.na(fstNeg_dt),is.na(lstNeg_dt))==FALSE) {
        x<-seq(as.Date(fstNeg_dt)-1,as.Date(lstNeg_dt)+1,by="1 day")
      }
      #Setting FD period to continous -ve diff days only
      THD_diff[length(THD_diff)]<-1
      
      L_MX<-0;    L_MX2<-0;    L_MX3<-0;    L_MX4<-0;    L_MX5<-0; L_MX6<-0  
      
      ##Only 1 -ve period
      if (length(which(THD_diff[x]>0))==1) {
        L_MX<-which(THD_diff[x]>0) #number  of -ve days
        L_Neg<- which(THD_diff[x]>0) #array of number  of -ve days
        L_MX_ind<- which.max(L_Neg) #index of maximum length of drought
        
        if (L_MX>=FDL) {
          
          x1.1<- gsub("-",'',fstNeg_dt)
          x2.1<- gsub("-",'',lstNeg_dt)
          y1.1<- cvar[paste(x1.1,'/',x2.1,sep='')]
          y2.1<- THD2[paste(x1.1,'/',x2.1,sep='')]
          THD_diff2.1<- y1.1-y2.1
          
          T1<- THD_diff[paste(x1.1,'/',x2.1,sep='')]
          T2<-THD_diff2.1
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          fstNeg_dt1<- index(T1[CT])[1]
          lstNeg_dt1<- index(T1[CT])[length(T1[CT])]
          
          L_MX<-ifelse(!any(CT)==TRUE,0,lstNeg_dt1 - fstNeg_dt1)
		  
		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX<-0
		  	}
		  }		  
		  
        } else {
          print("No FD events with the selected criteria")
        }
        
        ###More than 1 -ve period    
      } else if ((length(which(THD_diff[x]>0))>1)) {
        L_MX<-max(diff(which(THD_diff[x]>0))) #number  of -ve days
        L_Neg<- diff(which(THD_diff[x]>0)) #array of number  of -ve days
        L_MX_ind<- which.max(L_Neg) #index of maximum length of drought
        
        if (length(L_Neg)>1) {
          L_MX2<-sort(L_Neg,decreasing = T)[2] #length of second -ve event
          L_MX_ind2<-which(L_Neg==L_MX2)[which(L_Neg==L_MX2)!=L_MX_ind][1]   #index of second -ve event
        }
        if (length(L_Neg)>2) {
          L_MX3<-sort(L_Neg,decreasing = T)[3] #length of third -ve event
          L_MX_ind3<-which(L_Neg==L_MX3)[which(L_Neg==L_MX3)!=L_MX_ind2][1]  #index of third -ve event
        }
        if (length(L_Neg)>3) {
          L_MX4<-sort(L_Neg,decreasing = T)[4] #length of fourth -ve event  
          L_MX_ind4<-which(L_Neg==L_MX4)[which(L_Neg==L_MX4)!=L_MX_ind3][1]  #index of fourth -ve event
        }
        if (length(L_Neg)>4) {
          L_MX5<-sort(L_Neg,decreasing = T)[5] #length of fifth -ve event  
          L_MX_ind5<-which(L_Neg==L_MX5)[which(L_Neg==L_MX5)!=L_MX_ind4][1]  #index of fifth -ve event
        }
        if (length(L_Neg)>5) {
          L_MX6<-sort(L_Neg,decreasing = T)[6] #length of sixth -ve event
          L_MX_ind6<-which(L_Neg==L_MX6)[which(L_Neg==L_MX6)!=L_MX_ind4][1]  #index of sixth -ve event
        }
      }
        #Biggest FD Event
        if (L_MX>=FDL) {
          st_dum<- which(THD_diff[x]>0)[L_MX_ind][1]
          fstNeg_dt1<-fstNeg_dt+st_dum-0
          lstNeg_dt1<-fstNeg_dt1+L_MX-2
          
          fstNeg_dt.dum<-fstNeg_dt1
          
          x1.1<- gsub("-",'',fstNeg_dt1)
          x2.1<- gsub("-",'',lstNeg_dt1)
          y1.1<- cvar[paste(x1.1,'/',x2.1,sep='')]
          y2.1<- THD2[paste(x1.1,'/',x2.1,sep='')]
          THD_diff2.1<- y1.1-y2.1
          
          T1<- THD_diff[paste(x1.1,'/',x2.1,sep='')]
          T2<-THD_diff2.1
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          
          fstNeg_dt1<- index(T1[CT])[1]
          lstNeg_dt1<- index(T1[CT])[length(T1[CT])]  
          
          L_MX<-ifelse(!any(CT)==TRUE,0,lstNeg_dt1 - fstNeg_dt1)
		  
		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX<-0
		  	}
		  }
  
		  
        } else {
          print("No FD events with the selected criteria")
        }
        
        #Second FD Event (if available)
        if (L_MX2>=FDL) {
          st_dum2<- which(THD_diff[x]>0)[L_MX_ind2][1]
          fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
          fstNeg_dt2<-fstNeg_dt+st_dum2-1
          lstNeg_dt2<-fstNeg_dt2+L_MX2-2
          
          fstNeg_dt.dum<-fstNeg_dt2
          
          #Shading and event begining and end
          x1.2<- gsub("-",'',fstNeg_dt2)
          x2.2<- gsub("-",'',lstNeg_dt2)
          y1.2<- cvar[paste(x1.2,'/',x2.2,sep='')]        
          y2.2<- THD2[paste(x1.2,'/',x2.2,sep='')]
          THD_diff2.2<- y1.2-y2.2
          
          T1<- THD_diff[paste(x1.2,'/',x2.2,sep='')]
          T2<-THD_diff2.2
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          fstNeg_dt2<- index(T1[CT])[1]
          lstNeg_dt2<- index(T1[CT])[length(T1[CT])]  
          
          L_MX2<-ifelse(!any(CT)==TRUE,0,lstNeg_dt2 - fstNeg_dt2)

		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX2<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX2<-0
		  	}
		  }
		  
        }
        
        #Third FD Event (if available
        if (L_MX3>=FDL) {
          st_dum3<- which(THD_diff[x]>0)[L_MX_ind3][1]
          fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
          fstNeg_dt3<-fstNeg_dt+st_dum3-1
          lstNeg_dt3<-fstNeg_dt3+L_MX3-2
          
          fstNeg_dt.dum<-fstNeg_dt3
          
          #Shading and event begining and end
          x1.3<- gsub("-",'',fstNeg_dt3)
          x2.3<- gsub("-",'',lstNeg_dt3)
          y1.3<- cvar[paste(x1.3,'/',x2.3,sep='')]
          y2.3<- THD2[paste(x1.3,'/',x2.3,sep='')]
          THD_diff2.3<- y1.3-y2.3
          
          T1<- THD_diff[paste(x1.3,'/',x2.3,sep='')]
          T2<-THD_diff2.3
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
		  fstNeg_dt3<- index(T1[CT])[1]
          lstNeg_dt3<- index(T1[CT])[length(T1[CT])]  
                    
          L_MX3<-ifelse(!any(CT)==TRUE,0,lstNeg_dt3 - fstNeg_dt3)
 
 		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX3<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX3<-0
		  	}
		  }

		}
		
        #fourth FD Event (if available)
        if (L_MX4>=FDL) {
          st_dum4<- which(THD_diff[x]>0)[L_MX_ind4][1]
          fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
          fstNeg_dt4<-fstNeg_dt+st_dum4-1
          lstNeg_dt4<-fstNeg_dt4+L_MX4-2
          
          fstNeg_dt.dum<-fstNeg_dt4
          
          x1.4<- gsub("-",'',fstNeg_dt4)
          x2.4<- gsub("-",'',lstNeg_dt4)
          y1.4<- cvar[paste(x1.4,'/',x2.4,sep='')]
          y2.4<- THD2[paste(x1.4,'/',x2.4,sep='')]
          THD_diff2.4<- y1.4-y2.4
          
          T1<- THD_diff[paste(x1.4,'/',x2.4,sep='')]
          T2<-THD_diff2.4
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          fstNeg_dt4<- index(T1[CT])[1]
          lstNeg_dt4<- index(T1[CT])[length(T1[CT])]  
                    
          L_MX4<-ifelse(!any(CT)==TRUE,0,lstNeg_dt4 - fstNeg_dt4)

		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX4<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX4<-0
		  	}
		  }

        }
        
        #fifth FD Event (if available
        if (L_MX5>=FDL) {
          st_dum5<- which(THD_diff[x]>0)[L_MX_ind5][1]
          fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
          fstNeg_dt5<-fstNeg_dt+st_dum5-1
          lstNeg_dt5<-fstNeg_dt5+L_MX5-2
          
          fstNeg_dt.dum<-fstNeg_dt5
          
          x1.5<- gsub("-",'',fstNeg_dt5)
          x2.5<- gsub("-",'',lstNeg_dt5)
          y1.5<- cvar[paste(x1.5,'/',x2.5,sep='')]
          y2.5<- THD2[paste(x1.5,'/',x2.5,sep='')]
          THD_diff2.5<- y1.5-y2.5
          
          T1<- THD_diff[paste(x1.5,'/',x2.5,sep='')]
          T2<-THD_diff2.5
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          fstNeg_dt5<- index(T1[CT])[1]
          lstNeg_dt5<- index(T1[CT])[length(T1[CT])]  
                    
          L_MX5<-ifelse(!any(CT)==TRUE,0,lstNeg_dt5 - fstNeg_dt5)

		  if (BR_msk == T){
		  z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX5<-0
		  	}
		  }
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX5<-0
		  	}
		  }

        }
	
        #sixth FD Event (if available)
        if (L_MX6>=FDL) {
          st_dum6<- which(THD_diff[x]>0)[L_MX_ind6][1]
          fstNeg_dt<-index(THD_diff[THD_diff<0])[1]
          fstNeg_dt6<-fstNeg_dt+st_dum6-1
          lstNeg_dt6<-fstNeg_dt6+L_MX6-2
          
          fstNeg_dt.dum<-fstNeg_dt6
          
          x1.6<- gsub("-",'',fstNeg_dt6)
          x2.6<- gsub("-",'',lstNeg_dt6)
          y1.6<- cvar[paste(x1.6,'/',x2.6,sep='')]
          y2.6<- THD2[paste(x1.6,'/',x2.6,sep='')]
          THD_diff2.6<- y1.6-y2.6
          
          T1<- THD_diff[paste(x1.6,'/',x2.6,sep='')]
          T2<-THD_diff2.6
          T1[T1>0]=NA
          T2[T2>0]=NA
          CT<-complete.cases(T1,T2)
          fstNeg_dt6<- index(T1[CT])[1]
          lstNeg_dt6<- index(T1[CT])[length(T1[CT])]  
          
          L_MX6<-ifelse(!any(CT)==TRUE,0,lstNeg_dt6 - fstNeg_dt6)
	
			if (BR_msk == T){
		  	z1<- mean(msk_BR[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z1>BR_max | z1<BR_min) {
		  		L_MX6<-0
		  	}
		  }
			
			if (TS_msk == T){
		  	z2<- mean(msk_TS[paste(x1.1,'/',x2.1,sep='')],na.rm=T,trim=0.1)
		  	if (z2<TS_thd) {
		  		L_MX6<-0
		  	}
		  }


        }

      
      if (max(L_MX,L_MX2,L_MX3,L_MX4,L_MX5,L_MX6,na.rm = T)<FDL) {
        print('No FD Events with selected definition for the selected region')
        
      }
      
      #########FM Summarizing FDs#######
      #only 1 event
      if (length(which(THD_diff[x]>0))==1) { 
        if (L_MX>=FDL) {
          print(paste('FD Length =', L_MX))
          print(paste('from', fstNeg_dt,'to',lstNeg_dt))
          
          dm.df$L_MX1[ic]<- L_MX
          dm.df$fstdate1[ic]<- fstNeg_dt
          dm.df$lstdate1[ic]<- lstNeg_dt
          
        }
      }
      #More than 1 event
      ALL_dates<-c()
      if (L_MX>=FDL) {
        print(paste('FD-1 Length =', L_MX))
        print(paste('from', fstNeg_dt1,'to',lstNeg_dt1))
        ALL_dates<-c(fstNeg_dt1,lstNeg_dt1)
        
        dm.df$L_MX1[ic]<- L_MX
        dm.df$fstdate1[ic]<- fstNeg_dt1
        dm.df$lstdate1[ic]<- lstNeg_dt1
        
      }
      if (L_MX2>=FDL) {
        print(paste('FD-2 Length =', L_MX2))
        print(paste('from', fstNeg_dt2,'to',lstNeg_dt2))
        ALL_dates<-c(ALL_dates,fstNeg_dt2,lstNeg_dt2)
        
        dm.df$L_MX2[ic]<- L_MX2
        dm.df$fstdate2[ic]<- fstNeg_dt2
        dm.df$lstdate2[ic]<- lstNeg_dt2
        
      }
      if (L_MX3>=FDL) {
        print(paste('FD-3 Length =', L_MX3))
        print(paste('from', fstNeg_dt3,'to',lstNeg_dt3))
        ALL_dates<-c(ALL_dates,fstNeg_dt3,lstNeg_dt3)
        
        dm.df$L_MX3[ic]<- L_MX3
        dm.df$fstdate3[ic]<- fstNeg_dt3
        dm.df$lstdate3[ic]<- lstNeg_dt3
        
      }
      if (L_MX4>=FDL) {
        print(paste('FD-4 Length =', L_MX4))
        print(paste('from', fstNeg_dt4,'to',lstNeg_dt4))
        ALL_dates<-c(ALL_dates,fstNeg_dt4,lstNeg_dt4)
        
        dm.df$L_MX4[ic]<- L_MX4
        dm.df$fstdate4[ic]<- fstNeg_dt4
        dm.df$lstdate4[ic]<- lstNeg_dt4
        
      }
      if (L_MX5>=FDL) {
        print(paste('FD-5 Length =', L_MX5))
        print(paste('from', fstNeg_dt5,'to',lstNeg_dt5))
        ALL_dates<-c(ALL_dates,fstNeg_dt5,lstNeg_dt5)
        
        dm.df$L_MX5[ic]<- L_MX5
        dm.df$fstdate5[ic]<- fstNeg_dt5
        dm.df$lstdate5[ic]<- lstNeg_dt5
        
      }
      if (L_MX6>=FDL) {
        print(paste('FD-6 Length =', L_MX6))
        print(paste('from', fstNeg_dt6,'to',lstNeg_dt6))
        ALL_dates<-c(ALL_dates,fstNeg_dt6,lstNeg_dt6)
        
        dm.df$L_MX6[ic]<- L_MX6
        dm.df$fstdate6[ic]<- fstNeg_dt6
        dm.df$lstdate6[ic]<- lstNeg_dt6
        
      }
      
      
      if (max(L_MX,L_MX2,L_MX3,L_MX4,L_MX5,L_MX6,na.rm = T)<FDL) {
        print('No FD Events with selected definition for the selected region')
      }
      
      dm.df$FD_Max[ic]<-max(L_MX,L_MX2,L_MX3,L_MX4,L_MX5,L_MX6,na.rm=T)
	  
	 #Filling gaps
	 
    dm.df.2[ic,1:12]<-sort(as.matrix(dm.df[ic,7:18]),na.last = T)
	print("FILLING GAPS...")
for (ifill in 1:5) {
        if (!is.na(abs((as.numeric(dm.df.2$lstdate1[ic])-(as.numeric(dm.df.2$fstdate2[ic])))))){
          if (abs((as.numeric(dm.df.2$lstdate1[ic])-(as.numeric(dm.df.2$fstdate2[ic]))))<=FDgap) {
            dm.df.2$lstdate1[ic]<-dm.df.2$lstdate2[ic]
			dm.df.2$fstdate2[ic]<-dm.df.2$fstdate3[ic]
			dm.df.2$lstdate2[ic]<-dm.df.2$lstdate3[ic]
			dm.df.2$fstdate3[ic]<-dm.df.2$fstdate4[ic]
			dm.df.2$lstdate3[ic]<-dm.df.2$lstdate4[ic]
			dm.df.2$fstdate4[ic]<-dm.df.2$fstdate5[ic]
			dm.df.2$lstdate4[ic]<-dm.df.2$lstdate5[ic]
			dm.df.2$fstdate5[ic]<-dm.df.2$fstdate6[ic]
			dm.df.2$lstdate5[ic]<-dm.df.2$lstdate6[ic]
			dm.df.2$fstdate6[ic]<-dm.df.2$lstdate6[ic]<-NA
          }
        }
        
		if (!is.na(abs((as.numeric(dm.df.2$lstdate2[ic])-(as.numeric(dm.df.2$fstdate3[ic])))))){
          if (abs((as.numeric(dm.df.2$lstdate2[ic])-(as.numeric(dm.df.2$fstdate3[ic]))))<=FDgap) {
            dm.df.2$lstdate2[ic]<-dm.df.2$lstdate3[ic]
			dm.df.2$fstdate3[ic]<-dm.df.2$fstdate4[ic]
			dm.df.2$lstdate3[ic]<-dm.df.2$lstdate4[ic]
			dm.df.2$fstdate4[ic]<-dm.df.2$fstdate5[ic]
			dm.df.2$lstdate4[ic]<-dm.df.2$lstdate5[ic]
			dm.df.2$fstdate5[ic]<-dm.df.2$fstdate6[ic]
			dm.df.2$lstdate5[ic]<-dm.df.2$lstdate6[ic]
			dm.df.2$fstdate6[ic]<-dm.df.2$lstdate6[ic]<-NA
          }
        }
		
		if (!is.na(abs((as.numeric(dm.df.2$lstdate3[ic])-(as.numeric(dm.df.2$fstdate4[ic])))))){
          if (abs((as.numeric(dm.df.2$lstdate3[ic])-(as.numeric(dm.df.2$fstdate4[ic]))))<=FDgap) {
            dm.df.2$lstdate3[ic]<-dm.df.2$lstdate4[ic]
			dm.df.2$fstdate4[ic]<-dm.df.2$fstdate5[ic]
			dm.df.2$lstdate4[ic]<-dm.df.2$lstdate5[ic]
			dm.df.2$fstdate5[ic]<-dm.df.2$fstdate6[ic]
			dm.df.2$lstdate5[ic]<-dm.df.2$lstdate6[ic]
			dm.df.2$fstdate6[ic]<-dm.df.2$lstdate6[ic]<-NA
          }
        }
         
		if (!is.na(abs((as.numeric(dm.df.2$lstdate4[ic])-(as.numeric(dm.df.2$fstdate5[ic])))))){
          if (abs((as.numeric(dm.df.2$lstdate4[ic])-(as.numeric(dm.df.2$fstdate5[ic]))))<=FDgap) {
            dm.df.2$lstdate4[ic]<-dm.df.2$lstdate5[ic]
			dm.df.2$fstdate5[ic]<-dm.df.2$fstdate6[ic]
			dm.df.2$lstdate5[ic]<-dm.df.2$lstdate6[ic]
			dm.df.2$fstdate6[ic]<-dm.df.2$lstdate6[ic]<-NA
          }
        }
		if (!is.na(abs((as.numeric(dm.df.2$lstdate5[ic])-(as.numeric(dm.df.2$fstdate6[ic])))))){
          if (abs((as.numeric(dm.df.2$lstdate5[ic])-(as.numeric(dm.df.2$fstdate6[ic]))))<=FDgap) {
            dm.df.2$lstdate5[ic]<-dm.df.2$lstdate6[ic]
			dm.df.2$fstdate6[ic]<-dm.df.2$lstdate6[ic]<-NA
          }
        }
}
  
    #dm.df.2$FD_Max[ic]<- max(diff(as.numeric(dm.df.2[ic,1:2])),diff(as.numeric(dm.df.2[ic,3:4])),diff(as.numeric(dm.df.2[ic,5:6])),diff(as.numeric(dm.df.2[ic,7:8])),diff(as.numeric(dm.df.2[ic,9:10])),diff(as.numeric(dm.df.2[ic,11:12])),na.rm = T)
	#max_diff <- dm.df.2$FD_Max[ic]

		## TESTING NEW METHOD FOR MAX
    column_pairs <- list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12)
	for (pair in column_pairs) {
	  diff_values <- diff(as.numeric(dm.df.2[ic, pair]), na.rm = TRUE)
	  max_diff_in_pair <- max(diff_values, na.rm = TRUE)
	  
	  # Check if the current pair has a greater maximum difference
	  #if (max_diff_in_pair > max_diff) {
		max_diff <- max_diff_in_pair
		max_pair <- pair
	  #}
	}
	dm.df.2$FD_Max[ic] <- max_diff
	dm.df.2$fstdateMx[ic] <- dm.df.2[ic,max_pair[1]]
 	dm.df.2$lstdateMx[ic] <- dm.df.2[ic,max_pair[2]]
	
	
 ###Computing Severity:
 ##
 	print("COMPUTING SEVERITY...")
		#Event 1
		if (!is.na(dm.df.2$fstdate1[ic])) {
		
		x1.s<- gsub("-",'',dm.df.2$fstdate1[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate1[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s	#It's the opposite here. We need the +ve diff. between 20th and 5days RZSM.
        SV1 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV1[ic] <- SV1
		}
		
		#Event 2
		if (!is.na(dm.df.2$fstdate2[ic])) {
		
		x1.s<- gsub("-",'',dm.df.2$fstdate2[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate2[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SV2 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV2[ic] <- SV2
		}

		#Event 3
		if (!is.na(dm.df.2$fstdate3[ic])) {

		x1.s<- gsub("-",'',dm.df.2$fstdate3[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate3[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SV3 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV3[ic] <- SV3
		}

		#Event 4
		if (!is.na(dm.df.2$fstdate4[ic])) {

		x1.s<- gsub("-",'',dm.df.2$fstdate4[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate4[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SV4 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV4[ic] <- SV4
		}

		#Event 5
		if (!is.na(dm.df.2$fstdate5[ic])) {

		x1.s<- gsub("-",'',dm.df.2$fstdate5[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate5[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SV5 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV5[ic] <- SV5
		}

		#Event 6
		if (!is.na(dm.df.2$fstdate6[ic])) {

		x1.s<- gsub("-",'',dm.df.2$fstdate6[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdate6[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SV6 <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SV6[ic] <- SV6
		}
 		#Event MX
		if (!is.na(dm.df.2$fstdateMx[ic])) {

		x1.s<- gsub("-",'',dm.df.2$fstdateMx[ic])
        x2.s<- gsub("-",'',dm.df.2$lstdateMx[ic])
        y1.s<- cvar[paste(x1.s,'/',x2.s,sep='')]
        y2.s<- THD2[paste(x1.s,'/',x2.s,sep='')]
        THD_diff2.s<- y2.s-y1.s
        SVMx <- sum(THD_diff2.s[THD_diff2.s>0])
		dm.df.2$SVMx[ic] <- SVMx
		}

	#vegetation class
	msk_VEG  <-  xts(VEG_arr[ic,],time_d)
	VEGID <- as.numeric(msk_VEG[1])
	#VEGID<- msk_VEG[gsub("-",'',dm.df.2$fstdate1[ic])]
	dm.df.2$VEGID[ic] <- VEGID

 
 
 if (COMPOSITES == T) {
		 		 	 print("PROCESSING COMPOSITES...")

		 #Compute composites:
		 
		 
		 			#msk_VAR<-  SMA((xts(var_arr[ic,],time_d)),n = 1) #change n for a different running average window
			
			msk_TMP  <-  xts(TMP_arr  [ic,],time_d)
  			msk_PEVPR<-  xts(PEVPR_arr[ic,],time_d)
  			msk_ARAIN<-  xts(ARAIN_arr[ic,],time_d)
  			msk_PRES <-  xts(PRES_arr [ic,],time_d)
  			msk_RZSM <-  xts(nc01.2 [ic,],time_d)		#FROM GLDAS
  			msk_EVP  <-  xts(EVP_arr  [ic,],time_d)
  			msk_WS   <-  xts(WS_arr   [ic,],time_d)
  			msk_VPD  <-  xts(VPD_arr  [ic,],time_d)

		 
		 
		 if (!is.na(dm.df.2$fstdate1[ic])) {
		 		 	 print("E1 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_TMP[ic,18]<- X_Onset
					DF_X_TMP[ic,19]<- X_End
					DF_X_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 00]){ DF_X_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 01]){ DF_X_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 02]){ DF_X_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 03]){ DF_X_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 04]){ DF_X_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 05]){ DF_X_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 06]){ DF_X_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 07]){ DF_X_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 08]){ DF_X_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 09]){ DF_X_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 10]){ DF_X_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 11]){ DF_X_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 12]){ DF_X_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 13]){ DF_X_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 14]){ DF_X_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 15]){ DF_X_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[ntime - 16]){ DF_X_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_TMP[ic,04:18]),na.rm=T)
			
					DF_Y_TMP[ic,1]<-L1N
					DF_Y_TMP[ic,2]<-L2N
					DF_Y_TMP[ic,3]<-L3N
					DF_Y_TMP[ic,4]<-X_Onset
					DF_Y_TMP[ic,5]<-X_End
					DF_Y_TMP[ic,6]<-L1P
					DF_Y_TMP[ic,7]<-L2P
					DF_Y_TMP[ic,8]<-L3P		  
					DF_Y_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_PEVPR[ic,18]<- X_Onset
					DF_X_PEVPR[ic,19]<- X_End
					DF_X_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y_PEVPR[ic,1]<-L1N
					DF_Y_PEVPR[ic,2]<-L2N
					DF_Y_PEVPR[ic,3]<-L3N
					DF_Y_PEVPR[ic,4]<-X_Onset
					DF_Y_PEVPR[ic,5]<-X_End
					DF_Y_PEVPR[ic,6]<-L1P
					DF_Y_PEVPR[ic,7]<-L2P
					DF_Y_PEVPR[ic,8]<-L3P		  
					DF_Y_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_ARAIN[ic,18]<- X_Onset
					DF_X_ARAIN[ic,19]<- X_End
					DF_X_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y_ARAIN[ic,1]<-L1N
					DF_Y_ARAIN[ic,2]<-L2N
					DF_Y_ARAIN[ic,3]<-L3N
					DF_Y_ARAIN[ic,4]<-X_Onset
					DF_Y_ARAIN[ic,5]<-X_End
					DF_Y_ARAIN[ic,6]<-L1P
					DF_Y_ARAIN[ic,7]<-L2P
					DF_Y_ARAIN[ic,8]<-L3P		  
					DF_Y_ARAIN[ic,9]<-On3Pn

					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_PRES[ic,18]<- X_Onset
					DF_X_PRES[ic,19]<- X_End
					DF_X_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_PRES[ic,04:18]),na.rm=T)
			
					DF_Y_PRES[ic,1]<-L1N
					DF_Y_PRES[ic,2]<-L2N
					DF_Y_PRES[ic,3]<-L3N
					DF_Y_PRES[ic,4]<-X_Onset
					DF_Y_PRES[ic,5]<-X_End
					DF_Y_PRES[ic,6]<-L1P
					DF_Y_PRES[ic,7]<-L2P
					DF_Y_PRES[ic,8]<-L3P		  
					DF_Y_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_RZSM[ic,18]<- X_Onset
					DF_X_RZSM[ic,19]<- X_End
					DF_X_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y_RZSM[ic,1]<-L1N
					DF_Y_RZSM[ic,2]<-L2N
					DF_Y_RZSM[ic,3]<-L3N
					DF_Y_RZSM[ic,4]<-X_Onset
					DF_Y_RZSM[ic,5]<-X_End
					DF_Y_RZSM[ic,6]<-L1P
					DF_Y_RZSM[ic,7]<-L2P
					DF_Y_RZSM[ic,8]<-L3P		  
					DF_Y_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_EVP[ic,18]<- X_Onset
					DF_X_EVP[ic,19]<- X_End
					DF_X_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X[ic,04:18]),na.rm=T)
			
					DF_Y_EVP[ic,1]<-L1N
					DF_Y_EVP[ic,2]<-L2N
					DF_Y_EVP[ic,3]<-L3N
					DF_Y_EVP[ic,4]<-X_Onset
					DF_Y_EVP[ic,5]<-X_End
					DF_Y_EVP[ic,6]<-L1P
					DF_Y_EVP[ic,7]<-L2P
					DF_Y_EVP[ic,8]<-L3P		  
					DF_Y_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_WS[ic,18]<- X_Onset
					DF_X_WS[ic,19]<- X_End
					DF_X_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_WS[ic,04:18]),na.rm=T)
			
					DF_Y_WS[ic,1]<-L1N
					DF_Y_WS[ic,2]<-L2N
					DF_Y_WS[ic,3]<-L3N
					DF_Y_WS[ic,4]<-X_Onset
					DF_Y_WS[ic,5]<-X_End
					DF_Y_WS[ic,6]<-L1P
					DF_Y_WS[ic,7]<-L2P
					DF_Y_WS[ic,8]<-L3P		  
					DF_Y_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic])]
					DF_X_VPD[ic,18]<- X_Onset
					DF_X_VPD[ic,19]<- X_End
					DF_X_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-1)]
					DF_X_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-2)]
					DF_X_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-3)]
					DF_X_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-4)]
					DF_X_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-5)]
					DF_X_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-6)]
					DF_X_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-7)]
					DF_X_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-8)]
					DF_X_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-9)]
					DF_X_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-10)]
					DF_X_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-11)]
					DF_X_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-12)]
					DF_X_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-13)]
					DF_X_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-14)]
					DF_X_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-15)]
					DF_X_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-16)]
					DF_X_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate1[ic]-17)]
					
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 00)]){ DF_X_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+1)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 01)]){ DF_X_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+2)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 02)]){ DF_X_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+3)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 03)]){ DF_X_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+4)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 04)]){ DF_X_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+5)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 05)]){ DF_X_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+6)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 06)]){ DF_X_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+7)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 07)]){ DF_X_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+8)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 08)]){ DF_X_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+9)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 09)]){ DF_X_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+10)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 10)]){ DF_X_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+11)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 11)]){ DF_X_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+12)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 12)]){ DF_X_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+13)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 13)]){ DF_X_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+14)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 14)]){ DF_X_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+15)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 15)]){ DF_X_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+16)] }  
					if (dm.df.2$lstdate1[ic]<time_d[(ntime - 16)]){ DF_X_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate1[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X_VPD[ic,04:18]),na.rm=T)
			
					DF_Y_VPD[ic,1]<-L1N
					DF_Y_VPD[ic,2]<-L2N
					DF_Y_VPD[ic,3]<-L3N
					DF_Y_VPD[ic,4]<-X_Onset
					DF_Y_VPD[ic,5]<-X_End
					DF_Y_VPD[ic,6]<-L1P
					DF_Y_VPD[ic,7]<-L2P
					DF_Y_VPD[ic,8]<-L3P		  
					DF_Y_VPD[ic,9]<-On3Pn
					
			}

		if (!is.na(dm.df.2$fstdate2[ic])) {
		 		 	 print("E2 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_TMP[ic,18]<- X_Onset
					DF_X2_TMP[ic,19]<- X_End
					DF_X2_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_TMP[ic,04:18]),na.rm=T)
			
					DF_Y2_TMP[ic,1]<-L1N
					DF_Y2_TMP[ic,2]<-L2N
					DF_Y2_TMP[ic,3]<-L3N
					DF_Y2_TMP[ic,4]<-X_Onset
					DF_Y2_TMP[ic,5]<-X_End
					DF_Y2_TMP[ic,6]<-L1P
					DF_Y2_TMP[ic,7]<-L2P
					DF_Y2_TMP[ic,8]<-L3P		  
					DF_Y2_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_PEVPR[ic,18]<- X_Onset
					DF_X2_PEVPR[ic,19]<- X_End
					DF_X2_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y2_PEVPR[ic,1]<-L1N
					DF_Y2_PEVPR[ic,2]<-L2N
					DF_Y2_PEVPR[ic,3]<-L3N
					DF_Y2_PEVPR[ic,4]<-X_Onset
					DF_Y2_PEVPR[ic,5]<-X_End
					DF_Y2_PEVPR[ic,6]<-L1P
					DF_Y2_PEVPR[ic,7]<-L2P
					DF_Y2_PEVPR[ic,8]<-L3P		  
					DF_Y2_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_ARAIN[ic,18]<- X_Onset
					DF_X2_ARAIN[ic,19]<- X_End
					DF_X2_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y2_ARAIN[ic,1]<-L1N
					DF_Y2_ARAIN[ic,2]<-L2N
					DF_Y2_ARAIN[ic,3]<-L3N
					DF_Y2_ARAIN[ic,4]<-X_Onset
					DF_Y2_ARAIN[ic,5]<-X_End
					DF_Y2_ARAIN[ic,6]<-L1P
					DF_Y2_ARAIN[ic,7]<-L2P
					DF_Y2_ARAIN[ic,8]<-L3P		  
					DF_Y2_ARAIN[ic,9]<-On3Pn

					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_PRES[ic,18]<- X_Onset
					DF_X2_PRES[ic,19]<- X_End
					DF_X2_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_PRES[ic,04:18]),na.rm=T)
			
					DF_Y2_PRES[ic,1]<-L1N
					DF_Y2_PRES[ic,2]<-L2N
					DF_Y2_PRES[ic,3]<-L3N
					DF_Y2_PRES[ic,4]<-X_Onset
					DF_Y2_PRES[ic,5]<-X_End
					DF_Y2_PRES[ic,6]<-L1P
					DF_Y2_PRES[ic,7]<-L2P
					DF_Y2_PRES[ic,8]<-L3P		  
					DF_Y2_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_RZSM[ic,18]<- X_Onset
					DF_X2_RZSM[ic,19]<- X_End
					DF_X2_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y2_RZSM[ic,1]<-L1N
					DF_Y2_RZSM[ic,2]<-L2N
					DF_Y2_RZSM[ic,3]<-L3N
					DF_Y2_RZSM[ic,4]<-X_Onset
					DF_Y2_RZSM[ic,5]<-X_End
					DF_Y2_RZSM[ic,6]<-L1P
					DF_Y2_RZSM[ic,7]<-L2P
					DF_Y2_RZSM[ic,8]<-L3P		  
					DF_Y2_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_EVP[ic,18]<- X_Onset
					DF_X2_EVP[ic,19]<- X_End
					DF_X2_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_EVP[ic,04:18]),na.rm=T)
			
					DF_Y2_EVP[ic,1]<-L1N
					DF_Y2_EVP[ic,2]<-L2N
					DF_Y2_EVP[ic,3]<-L3N
					DF_Y2_EVP[ic,4]<-X_Onset
					DF_Y2_EVP[ic,5]<-X_End
					DF_Y2_EVP[ic,6]<-L1P
					DF_Y2_EVP[ic,7]<-L2P
					DF_Y2_EVP[ic,8]<-L3P		  
					DF_Y2_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_WS[ic,18]<- X_Onset
					DF_X2_WS[ic,19]<- X_End
					DF_X2_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_WS[ic,04:18]),na.rm=T)
			
					DF_Y2_WS[ic,1]<-L1N
					DF_Y2_WS[ic,2]<-L2N
					DF_Y2_WS[ic,3]<-L3N
					DF_Y2_WS[ic,4]<-X_Onset
					DF_Y2_WS[ic,5]<-X_End
					DF_Y2_WS[ic,6]<-L1P
					DF_Y2_WS[ic,7]<-L2P
					DF_Y2_WS[ic,8]<-L3P		  
					DF_Y2_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic])]
					DF_X2_VPD[ic,18]<- X_Onset
					DF_X2_VPD[ic,19]<- X_End
					DF_X2_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-1)]
					DF_X2_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-2)]
					DF_X2_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-3)]
					DF_X2_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-4)]
					DF_X2_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-5)]
					DF_X2_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-6)]
					DF_X2_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-7)]
					DF_X2_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-8)]
					DF_X2_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-9)]
					DF_X2_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-10)]
					DF_X2_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-11)]
					DF_X2_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-12)]
					DF_X2_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-13)]
					DF_X2_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-14)]
					DF_X2_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-15)]
					DF_X2_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-16)]
					DF_X2_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate2[ic]-17)]
					
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 00)]){ DF_X2_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+1)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 01)]){ DF_X2_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+2)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 02)]){ DF_X2_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+3)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 03)]){ DF_X2_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+4)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 04)]){ DF_X2_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+5)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 05)]){ DF_X2_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+6)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 06)]){ DF_X2_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+7)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 07)]){ DF_X2_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+8)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 08)]){ DF_X2_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+9)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 09)]){ DF_X2_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+10)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 10)]){ DF_X2_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+11)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 11)]){ DF_X2_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+12)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 12)]){ DF_X2_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+13)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 13)]){ DF_X2_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+14)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 14)]){ DF_X2_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+15)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 15)]){ DF_X2_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+16)] }  
					if (dm.df.2$lstdate2[ic]<time_d[(ntime - 16)]){ DF_X2_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate2[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X2_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X2_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X2_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X2_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X2_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X2_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X2_VPD[ic,04:18]),na.rm=T)
			
					DF_Y2_VPD[ic,1]<-L1N
					DF_Y2_VPD[ic,2]<-L2N
					DF_Y2_VPD[ic,3]<-L3N
					DF_Y2_VPD[ic,4]<-X_Onset
					DF_Y2_VPD[ic,5]<-X_End
					DF_Y2_VPD[ic,6]<-L1P
					DF_Y2_VPD[ic,7]<-L2P
					DF_Y2_VPD[ic,8]<-L3P		  
					DF_Y2_VPD[ic,9]<-On3Pn

			}

		if (!is.na(dm.df.2$fstdate3[ic])) {
		 		 	 print("E3 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_TMP[ic,18]<- X_Onset
					DF_X3_TMP[ic,19]<- X_End
					DF_X3_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_TMP[ic,04:18]),na.rm=T)
			
					DF_Y3_TMP[ic,1]<-L1N
					DF_Y3_TMP[ic,2]<-L2N
					DF_Y3_TMP[ic,3]<-L3N
					DF_Y3_TMP[ic,4]<-X_Onset
					DF_Y3_TMP[ic,5]<-X_End
					DF_Y3_TMP[ic,6]<-L1P
					DF_Y3_TMP[ic,7]<-L2P
					DF_Y3_TMP[ic,8]<-L3P		  
					DF_Y3_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_PEVPR[ic,18]<- X_Onset
					DF_X3_PEVPR[ic,19]<- X_End
					DF_X3_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y3_PEVPR[ic,1]<-L1N
					DF_Y3_PEVPR[ic,2]<-L2N
					DF_Y3_PEVPR[ic,3]<-L3N
					DF_Y3_PEVPR[ic,4]<-X_Onset
					DF_Y3_PEVPR[ic,5]<-X_End
					DF_Y3_PEVPR[ic,6]<-L1P
					DF_Y3_PEVPR[ic,7]<-L2P
					DF_Y3_PEVPR[ic,8]<-L3P		  
					DF_Y3_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_ARAIN[ic,18]<- X_Onset
					DF_X3_ARAIN[ic,19]<- X_End
					DF_X3_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y3_ARAIN[ic,1]<-L1N
					DF_Y3_ARAIN[ic,2]<-L2N
					DF_Y3_ARAIN[ic,3]<-L3N
					DF_Y3_ARAIN[ic,4]<-X_Onset
					DF_Y3_ARAIN[ic,5]<-X_End
					DF_Y3_ARAIN[ic,6]<-L1P
					DF_Y3_ARAIN[ic,7]<-L2P
					DF_Y3_ARAIN[ic,8]<-L3P		  
					DF_Y3_ARAIN[ic,9]<-On3Pn


					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_PRES[ic,18]<- X_Onset
					DF_X3_PRES[ic,19]<- X_End
					DF_X3_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_PRES[ic,04:18]),na.rm=T)
			
					DF_Y3_PRES[ic,1]<-L1N
					DF_Y3_PRES[ic,2]<-L2N
					DF_Y3_PRES[ic,3]<-L3N
					DF_Y3_PRES[ic,4]<-X_Onset
					DF_Y3_PRES[ic,5]<-X_End
					DF_Y3_PRES[ic,6]<-L1P
					DF_Y3_PRES[ic,7]<-L2P
					DF_Y3_PRES[ic,8]<-L3P		  
					DF_Y3_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_RZSM[ic,18]<- X_Onset
					DF_X3_RZSM[ic,19]<- X_End
					DF_X3_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y3_RZSM[ic,1]<-L1N
					DF_Y3_RZSM[ic,2]<-L2N
					DF_Y3_RZSM[ic,3]<-L3N
					DF_Y3_RZSM[ic,4]<-X_Onset
					DF_Y3_RZSM[ic,5]<-X_End
					DF_Y3_RZSM[ic,6]<-L1P
					DF_Y3_RZSM[ic,7]<-L2P
					DF_Y3_RZSM[ic,8]<-L3P		  
					DF_Y3_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_EVP[ic,18]<- X_Onset
					DF_X3_EVP[ic,19]<- X_End
					DF_X3_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_EVP[ic,04:18]),na.rm=T)
			
					DF_Y3_EVP[ic,1]<-L1N
					DF_Y3_EVP[ic,2]<-L2N
					DF_Y3_EVP[ic,3]<-L3N
					DF_Y3_EVP[ic,4]<-X_Onset
					DF_Y3_EVP[ic,5]<-X_End
					DF_Y3_EVP[ic,6]<-L1P
					DF_Y3_EVP[ic,7]<-L2P
					DF_Y3_EVP[ic,8]<-L3P		  
					DF_Y3_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_WS[ic,18]<- X_Onset
					DF_X3_WS[ic,19]<- X_End
					DF_X3_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_WS[ic,04:18]),na.rm=T)
			
					DF_Y3_WS[ic,1]<-L1N
					DF_Y3_WS[ic,2]<-L2N
					DF_Y3_WS[ic,3]<-L3N
					DF_Y3_WS[ic,4]<-X_Onset
					DF_Y3_WS[ic,5]<-X_End
					DF_Y3_WS[ic,6]<-L1P
					DF_Y3_WS[ic,7]<-L2P
					DF_Y3_WS[ic,8]<-L3P		  
					DF_Y3_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic])]
					DF_X3_VPD[ic,18]<- X_Onset
					DF_X3_VPD[ic,19]<- X_End
					DF_X3_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-1)]
					DF_X3_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-2)]
					DF_X3_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-3)]
					DF_X3_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-4)]
					DF_X3_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-5)]
					DF_X3_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-6)]
					DF_X3_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-7)]
					DF_X3_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-8)]
					DF_X3_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-9)]
					DF_X3_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-10)]
					DF_X3_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-11)]
					DF_X3_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-12)]
					DF_X3_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-13)]
					DF_X3_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-14)]
					DF_X3_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-15)]
					DF_X3_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-16)]
					DF_X3_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate3[ic]-17)]
					
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 00)]){ DF_X3_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+1)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 01)]){ DF_X3_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+2)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 02)]){ DF_X3_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+3)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 03)]){ DF_X3_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+4)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 04)]){ DF_X3_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+5)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 05)]){ DF_X3_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+6)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 06)]){ DF_X3_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+7)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 07)]){ DF_X3_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+8)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 08)]){ DF_X3_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+9)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 09)]){ DF_X3_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+10)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 10)]){ DF_X3_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+11)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 11)]){ DF_X3_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+12)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 12)]){ DF_X3_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+13)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 13)]){ DF_X3_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+14)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 14)]){ DF_X3_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+15)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 15)]){ DF_X3_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+16)] }  
					if (dm.df.2$lstdate3[ic]<time_d[(ntime - 16)]){ DF_X3_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate3[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X3_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X3_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X3_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X3_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X3_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X3_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X3_VPD[ic,04:18]),na.rm=T)
			
					DF_Y3_VPD[ic,1]<-L1N
					DF_Y3_VPD[ic,2]<-L2N
					DF_Y3_VPD[ic,3]<-L3N
					DF_Y3_VPD[ic,4]<-X_Onset
					DF_Y3_VPD[ic,5]<-X_End
					DF_Y3_VPD[ic,6]<-L1P
					DF_Y3_VPD[ic,7]<-L2P
					DF_Y3_VPD[ic,8]<-L3P		  
					DF_Y3_VPD[ic,9]<-On3Pn
					
			}

		if (!is.na(dm.df.2$fstdate4[ic])) {
		 		 	 print("E4 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_TMP[ic,18]<- X_Onset
					DF_X4_TMP[ic,19]<- X_End
					DF_X4_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_TMP[ic,04:18]),na.rm=T)
			
					DF_Y4_TMP[ic,1]<-L1N
					DF_Y4_TMP[ic,2]<-L2N
					DF_Y4_TMP[ic,3]<-L3N
					DF_Y4_TMP[ic,4]<-X_Onset
					DF_Y4_TMP[ic,5]<-X_End
					DF_Y4_TMP[ic,6]<-L1P
					DF_Y4_TMP[ic,7]<-L2P
					DF_Y4_TMP[ic,8]<-L3P		  
					DF_Y4_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_PEVPR[ic,18]<- X_Onset
					DF_X4_PEVPR[ic,19]<- X_End
					DF_X4_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y4_PEVPR[ic,1]<-L1N
					DF_Y4_PEVPR[ic,2]<-L2N
					DF_Y4_PEVPR[ic,3]<-L3N
					DF_Y4_PEVPR[ic,4]<-X_Onset
					DF_Y4_PEVPR[ic,5]<-X_End
					DF_Y4_PEVPR[ic,6]<-L1P
					DF_Y4_PEVPR[ic,7]<-L2P
					DF_Y4_PEVPR[ic,8]<-L3P		  
					DF_Y4_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_ARAIN[ic,18]<- X_Onset
					DF_X4_ARAIN[ic,19]<- X_End
					DF_X4_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y4_ARAIN[ic,1]<-L1N
					DF_Y4_ARAIN[ic,2]<-L2N
					DF_Y4_ARAIN[ic,3]<-L3N
					DF_Y4_ARAIN[ic,4]<-X_Onset
					DF_Y4_ARAIN[ic,5]<-X_End
					DF_Y4_ARAIN[ic,6]<-L1P
					DF_Y4_ARAIN[ic,7]<-L2P
					DF_Y4_ARAIN[ic,8]<-L3P		  
					DF_Y4_ARAIN[ic,9]<-On3Pn

					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_PRES[ic,18]<- X_Onset
					DF_X4_PRES[ic,19]<- X_End
					DF_X4_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_PRES[ic,04:18]),na.rm=T)
			
					DF_Y4_PRES[ic,1]<-L1N
					DF_Y4_PRES[ic,2]<-L2N
					DF_Y4_PRES[ic,3]<-L3N
					DF_Y4_PRES[ic,4]<-X_Onset
					DF_Y4_PRES[ic,5]<-X_End
					DF_Y4_PRES[ic,6]<-L1P
					DF_Y4_PRES[ic,7]<-L2P
					DF_Y4_PRES[ic,8]<-L3P		  
					DF_Y4_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_RZSM[ic,18]<- X_Onset
					DF_X4_RZSM[ic,19]<- X_End
					DF_X4_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y4_RZSM[ic,1]<-L1N
					DF_Y4_RZSM[ic,2]<-L2N
					DF_Y4_RZSM[ic,3]<-L3N
					DF_Y4_RZSM[ic,4]<-X_Onset
					DF_Y4_RZSM[ic,5]<-X_End
					DF_Y4_RZSM[ic,6]<-L1P
					DF_Y4_RZSM[ic,7]<-L2P
					DF_Y4_RZSM[ic,8]<-L3P		  
					DF_Y4_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_EVP[ic,18]<- X_Onset
					DF_X4_EVP[ic,19]<- X_End
					DF_X4_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_EVP[ic,04:18]),na.rm=T)
			
					DF_Y4_EVP[ic,1]<-L1N
					DF_Y4_EVP[ic,2]<-L2N
					DF_Y4_EVP[ic,3]<-L3N
					DF_Y4_EVP[ic,4]<-X_Onset
					DF_Y4_EVP[ic,5]<-X_End
					DF_Y4_EVP[ic,6]<-L1P
					DF_Y4_EVP[ic,7]<-L2P
					DF_Y4_EVP[ic,8]<-L3P		  
					DF_Y4_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_WS[ic,18]<- X_Onset
					DF_X4_WS[ic,19]<- X_End
					DF_X4_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_WS[ic,04:18]),na.rm=T)
			
					DF_Y4_WS[ic,1]<-L1N
					DF_Y4_WS[ic,2]<-L2N
					DF_Y4_WS[ic,3]<-L3N
					DF_Y4_WS[ic,4]<-X_Onset
					DF_Y4_WS[ic,5]<-X_End
					DF_Y4_WS[ic,6]<-L1P
					DF_Y4_WS[ic,7]<-L2P
					DF_Y4_WS[ic,8]<-L3P		  
					DF_Y4_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic])]
					DF_X4_VPD[ic,18]<- X_Onset
					DF_X4_VPD[ic,19]<- X_End
					DF_X4_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-1)]
					DF_X4_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-2)]
					DF_X4_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-3)]
					DF_X4_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-4)]
					DF_X4_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-5)]
					DF_X4_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-6)]
					DF_X4_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-7)]
					DF_X4_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-8)]
					DF_X4_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-9)]
					DF_X4_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-10)]
					DF_X4_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-11)]
					DF_X4_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-12)]
					DF_X4_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-13)]
					DF_X4_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-14)]
					DF_X4_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-15)]
					DF_X4_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-16)]
					DF_X4_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate4[ic]-17)]
					
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 00)]){ DF_X4_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+1)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 01)]){ DF_X4_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+2)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 02)]){ DF_X4_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+3)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 03)]){ DF_X4_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+4)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 04)]){ DF_X4_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+5)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 05)]){ DF_X4_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+6)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 06)]){ DF_X4_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+7)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 07)]){ DF_X4_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+8)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 08)]){ DF_X4_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+9)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 09)]){ DF_X4_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+10)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 10)]){ DF_X4_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+11)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 11)]){ DF_X4_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+12)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 12)]){ DF_X4_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+13)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 13)]){ DF_X4_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+14)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 14)]){ DF_X4_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+15)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 15)]){ DF_X4_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+16)] }  
					if (dm.df.2$lstdate4[ic]<time_d[(ntime - 16)]){ DF_X4_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate4[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X4_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X4_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X4_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X4_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X4_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X4_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X4_VPD[ic,04:18]),na.rm=T)
			
					DF_Y4_VPD[ic,1]<-L1N
					DF_Y4_VPD[ic,2]<-L2N
					DF_Y4_VPD[ic,3]<-L3N
					DF_Y4_VPD[ic,4]<-X_Onset
					DF_Y4_VPD[ic,5]<-X_End
					DF_Y4_VPD[ic,6]<-L1P
					DF_Y4_VPD[ic,7]<-L2P
					DF_Y4_VPD[ic,8]<-L3P		  
					DF_Y4_VPD[ic,9]<-On3Pn
					
			}

		if (!is.na(dm.df.2$fstdate5[ic])) {
		 		 	 print("E5 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_TMP[ic,18]<- X_Onset
					DF_X5_TMP[ic,19]<- X_End
					DF_X5_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_TMP[ic,04:18]),na.rm=T)
			
					DF_Y5_TMP[ic,1]<-L1N
					DF_Y5_TMP[ic,2]<-L2N
					DF_Y5_TMP[ic,3]<-L3N
					DF_Y5_TMP[ic,4]<-X_Onset
					DF_Y5_TMP[ic,5]<-X_End
					DF_Y5_TMP[ic,6]<-L1P
					DF_Y5_TMP[ic,7]<-L2P
					DF_Y5_TMP[ic,8]<-L3P		  
					DF_Y5_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_PEVPR[ic,18]<- X_Onset
					DF_X5_PEVPR[ic,19]<- X_End
					DF_X5_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y5_PEVPR[ic,1]<-L1N
					DF_Y5_PEVPR[ic,2]<-L2N
					DF_Y5_PEVPR[ic,3]<-L3N
					DF_Y5_PEVPR[ic,4]<-X_Onset
					DF_Y5_PEVPR[ic,5]<-X_End
					DF_Y5_PEVPR[ic,6]<-L1P
					DF_Y5_PEVPR[ic,7]<-L2P
					DF_Y5_PEVPR[ic,8]<-L3P		  
					DF_Y5_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_ARAIN[ic,18]<- X_Onset
					DF_X5_ARAIN[ic,19]<- X_End
					DF_X5_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y5_ARAIN[ic,1]<-L1N
					DF_Y5_ARAIN[ic,2]<-L2N
					DF_Y5_ARAIN[ic,3]<-L3N
					DF_Y5_ARAIN[ic,4]<-X_Onset
					DF_Y5_ARAIN[ic,5]<-X_End
					DF_Y5_ARAIN[ic,6]<-L1P
					DF_Y5_ARAIN[ic,7]<-L2P
					DF_Y5_ARAIN[ic,8]<-L3P		  
					DF_Y5_ARAIN[ic,9]<-On3Pn

					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_PRES[ic,18]<- X_Onset
					DF_X5_PRES[ic,19]<- X_End
					DF_X5_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_PRES[ic,04:18]),na.rm=T)
			
					DF_Y5_PRES[ic,1]<-L1N
					DF_Y5_PRES[ic,2]<-L2N
					DF_Y5_PRES[ic,3]<-L3N
					DF_Y5_PRES[ic,4]<-X_Onset
					DF_Y5_PRES[ic,5]<-X_End
					DF_Y5_PRES[ic,6]<-L1P
					DF_Y5_PRES[ic,7]<-L2P
					DF_Y5_PRES[ic,8]<-L3P		  
					DF_Y5_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_RZSM[ic,18]<- X_Onset
					DF_X5_RZSM[ic,19]<- X_End
					DF_X5_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y5_RZSM[ic,1]<-L1N
					DF_Y5_RZSM[ic,2]<-L2N
					DF_Y5_RZSM[ic,3]<-L3N
					DF_Y5_RZSM[ic,4]<-X_Onset
					DF_Y5_RZSM[ic,5]<-X_End
					DF_Y5_RZSM[ic,6]<-L1P
					DF_Y5_RZSM[ic,7]<-L2P
					DF_Y5_RZSM[ic,8]<-L3P		  
					DF_Y5_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_EVP[ic,18]<- X_Onset
					DF_X5_EVP[ic,19]<- X_End
					DF_X5_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_EVP[ic,04:18]),na.rm=T)
			
					DF_Y5_EVP[ic,1]<-L1N
					DF_Y5_EVP[ic,2]<-L2N
					DF_Y5_EVP[ic,3]<-L3N
					DF_Y5_EVP[ic,4]<-X_Onset
					DF_Y5_EVP[ic,5]<-X_End
					DF_Y5_EVP[ic,6]<-L1P
					DF_Y5_EVP[ic,7]<-L2P
					DF_Y5_EVP[ic,8]<-L3P		  
					DF_Y5_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_WS[ic,18]<- X_Onset
					DF_X5_WS[ic,19]<- X_End
					DF_X5_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_WS[ic,04:18]),na.rm=T)
			
					DF_Y5_WS[ic,1]<-L1N
					DF_Y5_WS[ic,2]<-L2N
					DF_Y5_WS[ic,3]<-L3N
					DF_Y5_WS[ic,4]<-X_Onset
					DF_Y5_WS[ic,5]<-X_End
					DF_Y5_WS[ic,6]<-L1P
					DF_Y5_WS[ic,7]<-L2P
					DF_Y5_WS[ic,8]<-L3P		  
					DF_Y5_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic])]
					DF_X5_VPD[ic,18]<- X_Onset
					DF_X5_VPD[ic,19]<- X_End
					DF_X5_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-1)]
					DF_X5_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-2)]
					DF_X5_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-3)]
					DF_X5_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-4)]
					DF_X5_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-5)]
					DF_X5_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-6)]
					DF_X5_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-7)]
					DF_X5_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-8)]
					DF_X5_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-9)]
					DF_X5_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-10)]
					DF_X5_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-11)]
					DF_X5_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-12)]
					DF_X5_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-13)]
					DF_X5_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-14)]
					DF_X5_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-15)]
					DF_X5_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-16)]
					DF_X5_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate5[ic]-17)]
					
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 00)]){ DF_X5_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+1)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 01)]){ DF_X5_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+2)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 02)]){ DF_X5_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+3)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 03)]){ DF_X5_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+4)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 04)]){ DF_X5_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+5)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 05)]){ DF_X5_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+6)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 06)]){ DF_X5_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+7)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 07)]){ DF_X5_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+8)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 08)]){ DF_X5_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+9)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 09)]){ DF_X5_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+10)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 10)]){ DF_X5_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+11)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 11)]){ DF_X5_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+12)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 12)]){ DF_X5_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+13)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 13)]){ DF_X5_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+14)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 14)]){ DF_X5_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+15)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 15)]){ DF_X5_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+16)] }  
					if (dm.df.2$lstdate5[ic]<time_d[(ntime - 16)]){ DF_X5_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate5[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X5_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X5_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X5_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X5_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X5_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X5_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X5_VPD[ic,04:18]),na.rm=T)
			
					DF_Y5_VPD[ic,1]<-L1N
					DF_Y5_VPD[ic,2]<-L2N
					DF_Y5_VPD[ic,3]<-L3N
					DF_Y5_VPD[ic,4]<-X_Onset
					DF_Y5_VPD[ic,5]<-X_End
					DF_Y5_VPD[ic,6]<-L1P
					DF_Y5_VPD[ic,7]<-L2P
					DF_Y5_VPD[ic,8]<-L3P		  
					DF_Y5_VPD[ic,9]<-On3Pn

			}

		if (!is.na(dm.df.2$fstdate6[ic])) {
		 	 print("E6 - COMPOSITES...")

					X_Onset <- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_TMP[ic,18]<- X_Onset
					DF_X6_TMP[ic,19]<- X_End
					DF_X6_TMP[ic,17]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_TMP[ic,16]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_TMP[ic,15]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_TMP[ic,14]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_TMP[ic,13]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_TMP[ic,12]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_TMP[ic,11]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_TMP[ic,10]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_TMP[ic,09]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_TMP[ic,08]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_TMP[ic,07]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_TMP[ic,06]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_TMP[ic,05]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_TMP[ic,04]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_TMP[ic,03]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_TMP[ic,02]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_TMP[ic,01]<- msk_TMP[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_TMP[ic,20]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_TMP[ic,21]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_TMP[ic,22]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_TMP[ic,23]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_TMP[ic,24]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_TMP[ic,25]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_TMP[ic,26]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_TMP[ic,27]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_TMP[ic,28]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_TMP[ic,29]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_TMP[ic,30]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_TMP[ic,31]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_TMP[ic,32]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_TMP[ic,33]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_TMP[ic,34]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_TMP[ic,35]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_TMP[ic,36]<- msk_TMP[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_TMP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_TMP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_TMP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_TMP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_TMP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_TMP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_TMP[ic,04:18]),na.rm=T)
			
					DF_Y6_TMP[ic,1]<-L1N
					DF_Y6_TMP[ic,2]<-L2N
					DF_Y6_TMP[ic,3]<-L3N
					DF_Y6_TMP[ic,4]<-X_Onset
					DF_Y6_TMP[ic,5]<-X_End
					DF_Y6_TMP[ic,6]<-L1P
					DF_Y6_TMP[ic,7]<-L2P
					DF_Y6_TMP[ic,8]<-L3P		  
					DF_Y6_TMP[ic,9]<-On3Pn
				
					X_Onset <- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_PEVPR[ic,18]<- X_Onset
					DF_X6_PEVPR[ic,19]<- X_End
					DF_X6_PEVPR[ic,17]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_PEVPR[ic,16]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_PEVPR[ic,15]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_PEVPR[ic,14]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_PEVPR[ic,13]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_PEVPR[ic,12]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_PEVPR[ic,11]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_PEVPR[ic,10]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_PEVPR[ic,09]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_PEVPR[ic,08]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_PEVPR[ic,07]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_PEVPR[ic,06]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_PEVPR[ic,05]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_PEVPR[ic,04]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_PEVPR[ic,03]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_PEVPR[ic,02]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_PEVPR[ic,01]<- msk_PEVPR[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_PEVPR[ic,20]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_PEVPR[ic,21]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_PEVPR[ic,22]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_PEVPR[ic,23]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_PEVPR[ic,24]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_PEVPR[ic,25]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_PEVPR[ic,26]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_PEVPR[ic,27]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_PEVPR[ic,28]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_PEVPR[ic,29]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_PEVPR[ic,30]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_PEVPR[ic,31]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_PEVPR[ic,32]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_PEVPR[ic,33]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_PEVPR[ic,34]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_PEVPR[ic,35]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_PEVPR[ic,36]<- msk_PEVPR[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_PEVPR[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_PEVPR[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_PEVPR[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_PEVPR[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_PEVPR[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_PEVPR[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_PEVPR[ic,04:18]),na.rm=T)
			
					DF_Y6_PEVPR[ic,1]<-L1N
					DF_Y6_PEVPR[ic,2]<-L2N
					DF_Y6_PEVPR[ic,3]<-L3N
					DF_Y6_PEVPR[ic,4]<-X_Onset
					DF_Y6_PEVPR[ic,5]<-X_End
					DF_Y6_PEVPR[ic,6]<-L1P
					DF_Y6_PEVPR[ic,7]<-L2P
					DF_Y6_PEVPR[ic,8]<-L3P		  
					DF_Y6_PEVPR[ic,9]<-On3Pn

					X_Onset <- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_ARAIN[ic,18]<- X_Onset
					DF_X6_ARAIN[ic,19]<- X_End
					DF_X6_ARAIN[ic,17]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_ARAIN[ic,16]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_ARAIN[ic,15]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_ARAIN[ic,14]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_ARAIN[ic,13]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_ARAIN[ic,12]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_ARAIN[ic,11]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_ARAIN[ic,10]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_ARAIN[ic,09]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_ARAIN[ic,08]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_ARAIN[ic,07]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_ARAIN[ic,06]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_ARAIN[ic,05]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_ARAIN[ic,04]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_ARAIN[ic,03]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_ARAIN[ic,02]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_ARAIN[ic,01]<- msk_ARAIN[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_ARAIN[ic,20]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_ARAIN[ic,21]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_ARAIN[ic,22]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_ARAIN[ic,23]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_ARAIN[ic,24]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_ARAIN[ic,25]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_ARAIN[ic,26]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_ARAIN[ic,27]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_ARAIN[ic,28]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_ARAIN[ic,29]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_ARAIN[ic,30]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_ARAIN[ic,31]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_ARAIN[ic,32]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_ARAIN[ic,33]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_ARAIN[ic,34]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_ARAIN[ic,35]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_ARAIN[ic,36]<- msk_ARAIN[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_ARAIN[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_ARAIN[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_ARAIN[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_ARAIN[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_ARAIN[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_ARAIN[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_ARAIN[ic,04:18]),na.rm=T)
			
					DF_Y6_ARAIN[ic,1]<-L1N
					DF_Y6_ARAIN[ic,2]<-L2N
					DF_Y6_ARAIN[ic,3]<-L3N
					DF_Y6_ARAIN[ic,4]<-X_Onset
					DF_Y6_ARAIN[ic,5]<-X_End
					DF_Y6_ARAIN[ic,6]<-L1P
					DF_Y6_ARAIN[ic,7]<-L2P
					DF_Y6_ARAIN[ic,8]<-L3P		  
					DF_Y6_ARAIN[ic,9]<-On3Pn

					X_Onset <- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_PRES[ic,18]<- X_Onset
					DF_X6_PRES[ic,19]<- X_End
					DF_X6_PRES[ic,17]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_PRES[ic,16]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_PRES[ic,15]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_PRES[ic,14]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_PRES[ic,13]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_PRES[ic,12]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_PRES[ic,11]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_PRES[ic,10]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_PRES[ic,09]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_PRES[ic,08]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_PRES[ic,07]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_PRES[ic,06]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_PRES[ic,05]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_PRES[ic,04]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_PRES[ic,03]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_PRES[ic,02]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_PRES[ic,01]<- msk_PRES[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_PRES[ic,20]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_PRES[ic,21]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_PRES[ic,22]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_PRES[ic,23]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_PRES[ic,24]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_PRES[ic,25]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_PRES[ic,26]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_PRES[ic,27]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_PRES[ic,28]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_PRES[ic,29]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_PRES[ic,30]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_PRES[ic,31]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_PRES[ic,32]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_PRES[ic,33]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_PRES[ic,34]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_PRES[ic,35]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_PRES[ic,36]<- msk_PRES[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_PRES[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_PRES[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_PRES[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_PRES[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_PRES[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_PRES[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_PRES[ic,04:18]),na.rm=T)
			
					DF_Y6_PRES[ic,1]<-L1N
					DF_Y6_PRES[ic,2]<-L2N
					DF_Y6_PRES[ic,3]<-L3N
					DF_Y6_PRES[ic,4]<-X_Onset
					DF_Y6_PRES[ic,5]<-X_End
					DF_Y6_PRES[ic,6]<-L1P
					DF_Y6_PRES[ic,7]<-L2P
					DF_Y6_PRES[ic,8]<-L3P		  
					DF_Y6_PRES[ic,9]<-On3Pn
					
					X_Onset <- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_RZSM[ic,18]<- X_Onset
					DF_X6_RZSM[ic,19]<- X_End
					DF_X6_RZSM[ic,17]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_RZSM[ic,16]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_RZSM[ic,15]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_RZSM[ic,14]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_RZSM[ic,13]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_RZSM[ic,12]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_RZSM[ic,11]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_RZSM[ic,10]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_RZSM[ic,09]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_RZSM[ic,08]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_RZSM[ic,07]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_RZSM[ic,06]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_RZSM[ic,05]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_RZSM[ic,04]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_RZSM[ic,03]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_RZSM[ic,02]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_RZSM[ic,01]<- msk_RZSM[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_RZSM[ic,20]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_RZSM[ic,21]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_RZSM[ic,22]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_RZSM[ic,23]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_RZSM[ic,24]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_RZSM[ic,25]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_RZSM[ic,26]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_RZSM[ic,27]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_RZSM[ic,28]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_RZSM[ic,29]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_RZSM[ic,30]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_RZSM[ic,31]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_RZSM[ic,32]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_RZSM[ic,33]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_RZSM[ic,34]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_RZSM[ic,35]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_RZSM[ic,36]<- msk_RZSM[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_RZSM[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_RZSM[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_RZSM[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_RZSM[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_RZSM[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_RZSM[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_RZSM[ic,04:18]),na.rm=T)
			
					DF_Y6_RZSM[ic,1]<-L1N
					DF_Y6_RZSM[ic,2]<-L2N
					DF_Y6_RZSM[ic,3]<-L3N
					DF_Y6_RZSM[ic,4]<-X_Onset
					DF_Y6_RZSM[ic,5]<-X_End
					DF_Y6_RZSM[ic,6]<-L1P
					DF_Y6_RZSM[ic,7]<-L2P
					DF_Y6_RZSM[ic,8]<-L3P		  
					DF_Y6_RZSM[ic,9]<-On3Pn
					
					X_Onset <- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_EVP[ic,18]<- X_Onset
					DF_X6_EVP[ic,19]<- X_End
					DF_X6_EVP[ic,17]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_EVP[ic,16]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_EVP[ic,15]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_EVP[ic,14]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_EVP[ic,13]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_EVP[ic,12]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_EVP[ic,11]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_EVP[ic,10]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_EVP[ic,09]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_EVP[ic,08]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_EVP[ic,07]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_EVP[ic,06]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_EVP[ic,05]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_EVP[ic,04]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_EVP[ic,03]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_EVP[ic,02]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_EVP[ic,01]<- msk_EVP[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_EVP[ic,20]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_EVP[ic,21]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_EVP[ic,22]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_EVP[ic,23]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_EVP[ic,24]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_EVP[ic,25]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_EVP[ic,26]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_EVP[ic,27]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_EVP[ic,28]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_EVP[ic,29]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_EVP[ic,30]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_EVP[ic,31]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_EVP[ic,32]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_EVP[ic,33]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_EVP[ic,34]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_EVP[ic,35]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_EVP[ic,36]<- msk_EVP[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_EVP[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_EVP[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_EVP[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_EVP[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_EVP[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_EVP[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_EVP[ic,04:18]),na.rm=T)
			
					DF_Y6_EVP[ic,1]<-L1N
					DF_Y6_EVP[ic,2]<-L2N
					DF_Y6_EVP[ic,3]<-L3N
					DF_Y6_EVP[ic,4]<-X_Onset
					DF_Y6_EVP[ic,5]<-X_End
					DF_Y6_EVP[ic,6]<-L1P
					DF_Y6_EVP[ic,7]<-L2P
					DF_Y6_EVP[ic,8]<-L3P		  
					DF_Y6_EVP[ic,9]<-On3Pn
					
					X_Onset <- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_WS[ic,18]<- X_Onset
					DF_X6_WS[ic,19]<- X_End
					DF_X6_WS[ic,17]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_WS[ic,16]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_WS[ic,15]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_WS[ic,14]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_WS[ic,13]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_WS[ic,12]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_WS[ic,11]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_WS[ic,10]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_WS[ic,09]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_WS[ic,08]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_WS[ic,07]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_WS[ic,06]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_WS[ic,05]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_WS[ic,04]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_WS[ic,03]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_WS[ic,02]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_WS[ic,01]<- msk_WS[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_WS[ic,20]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_WS[ic,21]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_WS[ic,22]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_WS[ic,23]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_WS[ic,24]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_WS[ic,25]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_WS[ic,26]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_WS[ic,27]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_WS[ic,28]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_WS[ic,29]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_WS[ic,30]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_WS[ic,31]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_WS[ic,32]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_WS[ic,33]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_WS[ic,34]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_WS[ic,35]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_WS[ic,36]<- msk_WS[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_WS[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_WS[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_WS[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_WS[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_WS[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_WS[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_WS[ic,04:18]),na.rm=T)
			
					DF_Y6_WS[ic,1]<-L1N
					DF_Y6_WS[ic,2]<-L2N
					DF_Y6_WS[ic,3]<-L3N
					DF_Y6_WS[ic,4]<-X_Onset
					DF_Y6_WS[ic,5]<-X_End
					DF_Y6_WS[ic,6]<-L1P
					DF_Y6_WS[ic,7]<-L2P
					DF_Y6_WS[ic,8]<-L3P		  
					DF_Y6_WS[ic,9]<-On3Pn

					X_Onset <- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic])]
					X_End   <- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic])]
					DF_X6_VPD[ic,18]<- X_Onset
					DF_X6_VPD[ic,19]<- X_End
					DF_X6_VPD[ic,17]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-1)]
					DF_X6_VPD[ic,16]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-2)]
					DF_X6_VPD[ic,15]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-3)]
					DF_X6_VPD[ic,14]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-4)]
					DF_X6_VPD[ic,13]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-5)]
					DF_X6_VPD[ic,12]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-6)]
					DF_X6_VPD[ic,11]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-7)]
					DF_X6_VPD[ic,10]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-8)]
					DF_X6_VPD[ic,09]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-9)]
					DF_X6_VPD[ic,08]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-10)]
					DF_X6_VPD[ic,07]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-11)]
					DF_X6_VPD[ic,06]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-12)]
					DF_X6_VPD[ic,05]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-13)]
					DF_X6_VPD[ic,04]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-14)]
					DF_X6_VPD[ic,03]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-15)]
					DF_X6_VPD[ic,02]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-16)]
					DF_X6_VPD[ic,01]<- msk_VPD[gsub("-",'',dm.df.2$fstdate6[ic]-17)]
					
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 00)]){ DF_X6_VPD[ic,20]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+1)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 01)]){ DF_X6_VPD[ic,21]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+2)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 02)]){ DF_X6_VPD[ic,22]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+3)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 03)]){ DF_X6_VPD[ic,23]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+4)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 04)]){ DF_X6_VPD[ic,24]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+5)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 05)]){ DF_X6_VPD[ic,25]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+6)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 06)]){ DF_X6_VPD[ic,26]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+7)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 07)]){ DF_X6_VPD[ic,27]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+8)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 08)]){ DF_X6_VPD[ic,28]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+9)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 09)]){ DF_X6_VPD[ic,29]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+10)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 10)]){ DF_X6_VPD[ic,30]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+11)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 11)]){ DF_X6_VPD[ic,31]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+12)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 12)]){ DF_X6_VPD[ic,32]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+13)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 13)]){ DF_X6_VPD[ic,33]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+14)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 14)]){ DF_X6_VPD[ic,34]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+15)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 15)]){ DF_X6_VPD[ic,35]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+16)] }  
					if (dm.df.2$lstdate6[ic]<time_d[(ntime - 16)]){ DF_X6_VPD[ic,36]<- msk_VPD[gsub("-",'',dm.df.2$lstdate6[ic]+17)] }  
		  

					#Lag -1 Pentad Onset
					L1N<- mean(as.double(DF_X6_VPD[ic,14:18]),na.rm=T)
					#Lag -2 Pentad Onset
					L2N<- mean(as.double(DF_X6_VPD[ic,09:13]),na.rm=T)
					#Lag -3 Pentad Onset
					L3N<- mean(as.double(DF_X6_VPD[ic,04:08]),na.rm=T)
					
					#Lag +1 Pentad Recovery
					L1P<- mean(as.double(DF_X6_VPD[ic,19:23]),na.rm=T)
					#Lag +2 Pentad Recovery
					L2P<- mean(as.double(DF_X6_VPD[ic,24:28]),na.rm=T)
					#Lag +3 Pentad Recovery
					L3P<- mean(as.double(DF_X6_VPD[ic,29:33]),na.rm=T)
					
					#3-Pentads pre-onset composite
					On3Pn<- mean(as.double(DF_X6_VPD[ic,04:18]),na.rm=T)
			
					DF_Y6_VPD[ic,1]<-L1N
					DF_Y6_VPD[ic,2]<-L2N
					DF_Y6_VPD[ic,3]<-L3N
					DF_Y6_VPD[ic,4]<-X_Onset
					DF_Y6_VPD[ic,5]<-X_End
					DF_Y6_VPD[ic,6]<-L1P
					DF_Y6_VPD[ic,7]<-L2P
					DF_Y6_VPD[ic,8]<-L3P		  
					DF_Y6_VPD[ic,9]<-On3Pn

			}
	}
	 
	    ##define colors and legends etc...
    #
    #ncolr=length(seq(min(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)]), max(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)])))+2
    #col=rev(colorspace::sequential_hcl(ncolr,palette = "YlOrRd"))
    #at_cmap=seq(min(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)])-1, max(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)])+1)
    #ckey=list( at=at_cmap, col=col,raster=TRUE, interpolate=TRUE)
    #
    #
    #  tp <- matrix(dm.df.2$FD_Max, nrow=dim(nc01)[1], byrow = T)
    #  tp <- apply(tp, 2, rev)
    #  
    #  nc01x <- nc_open(paste('NLDAS/NLDAS_NOAH0125_D0.A',YYYY,'03_',YYYY,'11','.nc',sep = ''))
    #  
    #  lon <- ncvar_get(nc01x,var="lon")
    #  lon[lon > 180] <- lon[lon > 180] - 360
    #  lat <- ncvar_get(nc01x,var="lat")  
    #  
    #  nc_close(nc01x)
    #  LC <- 1
    #        
    #  wld <- maps::map('state',plot=F)
    #  wld <- data.frame(lon=wld$x, lat=wld$y)
    #  
    #  
    #  ncolr=length(seq(min(dm.df.2$FD_Max[which(dm.df.2$FD_Max>0)]), max(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)])))+2
    #  col=rev(colorspace::sequential_hcl(ncolr,palette = "YlOrRd"))
    #  at_cmap=seq(min(dm.df.2$FD_Max[which(dm.df.2$FD_Max>0)])-1, max(dm.df.2$FD_Max[which(dm.df.2$FD_Max>=0)])+1)
    #  ckey=list( at=at_cmap, col=col,raster=TRUE, interpolate=TRUE)
    #  
    #  phc2<- levelplot(t(tp)*LC,row.values = lon,column.values = lat,col.regions=col,xlab='Longitude',ylab='Latitude',main='',at=at_cmap,colorkey=ckey)+xyplot(lat ~ lon, wld, type='l', lty=1, lwd=1, col='black')
    #  png(filename=paste("'SMVI_NLDAS_out/GridedMapFD_",FD_mthd,'_',YYYY,'_n',n_reg,".png",sep=''), width = 2000, height = 1500,res=200,type='cairo') 
    #  if (transparent_BG==T){
    #    par(bg=NA)
    #  }
    #  plot(phc2)
    #  dev.off()
    #  plot(phc2)
	 }
    }  
    
#stopCluster(cl)
	
dm.df.2<-transform(dm.df.2, FD_Max = as.numeric(FD_Max))
dm.df.2$FD_Max[which(!is.finite(dm.df.2$FD_Max))]<-NA
   

###data<- dm.df.2
###data<- data[rowSums(is.na(data[,1:12])) != ncol(data[,1:12]),]

####Write output#####
#####################

#Set and create writing directory
OUT_DIR=paste('SMVI_',DS,'_out',sep='')
dir.create(file.path(OUT_DIR))
#FD events will be saved here:
dir.create(file.path(OUT_DIR,"FD_Events"))
if (COMPOSITES == T) {
	#Composites will be saved here:
	dir.create(file.path(OUT_DIR,"Composites"))
}

write.csv(dm.df.2,file = paste(OUT_DIR,'/FD_Events/SMVI_',DS,'_','Summ_table_',FD_mthd,'_',YYYY,'_n',n_reg,".csv",sep=''))

if (COMPOSITES == T) {
	 print("WRITING OUTPUT...")
	 write.csv(DF_Y_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"TMP",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"PEVPR",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"ARAIN",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"PRES",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"RZSM",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"EVP",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"WS",'_',YYYY,".csv",sep=''))
	 write.csv(DF_Y_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E1_',"VPD",'_',YYYY,".csv",sep=''))
																 
	write.csv(DF_Y2_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"TMP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"PEVPR",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"ARAIN",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"PRES",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"RZSM",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"EVP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"WS",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y2_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E2_',"VPD",'_',YYYY,".csv",sep=''))
																 
	write.csv(DF_Y3_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"TMP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"PEVPR",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"ARAIN",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"PRES",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"RZSM",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"EVP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"WS",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y3_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E3_',"VPD",'_',YYYY,".csv",sep=''))
																 
	write.csv(DF_Y4_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"TMP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"PEVPR",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"ARAIN",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"PRES",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"RZSM",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"EVP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"WS",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y4_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E4_',"VPD",'_',YYYY,".csv",sep=''))
																 
	write.csv(DF_Y5_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"TMP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"PEVPR",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"ARAIN",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"PRES",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"RZSM",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"EVP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"WS",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y5_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E5_',"VPD",'_',YYYY,".csv",sep=''))
																 
	write.csv(DF_Y6_TMP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"TMP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_PEVPR,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"PEVPR",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_ARAIN,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"ARAIN",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_PRES ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"PRES",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_RZSM ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"RZSM",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_EVP  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"EVP",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_WS   ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"WS",'_',YYYY,".csv",sep=''))
	write.csv(DF_Y6_VPD  ,file = paste(OUT_DIR,'/Composites/SMVI_',DS,'_E6_',"VPD",'_',YYYY,".csv",sep=''))
}


if (out_fig == T){

		#############
		##PLOTTING###
		print ("PLOTTING")
		wld <- maps::map('world',plot=F)
		wld <- data.frame(lon=wld$x, lat=wld$y)
		dir.create(file.path(OUT_DIR,"Figures"))
		FIG_DIR=paste(OUT_DIR,'/Figures/',sep='')


		#nc01x <- nc_open(paste('GLDAS/GLDAS_CLSM025_1980_2022_SEL_CLIM_RZSM.nc',sep = ''))

		#lon <- ncvar_get(nc01x,var="lon")
		#lon[lon > 180] <- lon[lon > 180] - 360
		#lat <- ncvar_get(nc01x,var="lat")

		#rm(nc01x)
		color_palette <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00", "gold1",
				   "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70")

		csv_in<- fread(file= paste(OUT_DIR,'/FD_Events/SMVI_',DS,'_','Summ_table_',FD_mthd,'_',YYYY,'_n',n_reg,".csv",sep=''))
		
		###EMx plot
		fd_YY <- as.data.frame(csv_in[,'fstdateMx'],colClasses = c("date"))
		#fd_YY_Mon <- as.POSIXlt(fd_YY*86400, origin= "1970-01-01")$mon +1
		fd_YY_Mon <- month(fd_YY[,])
		month_matrix <- array(fd_YY_Mon, dim = c(nlon, nlat))


		var_MAT <- array(fd_YY, dim = c(nlon, nlat))    #Converts the vector to matrix of saze nlon x nlat
		month_names <- month.abb

		Plt<-levelplot(month_matrix,
						 row.values = lon,
						 column.values = lat,
						 xlab = 'Longitude',
						 ylab = 'Latitude',
						 ylim = c(-60, 90),
						 xlim = c(-180, 180),
						 main = paste(YYYY),
						 col.regions = color_palette[1:12],
						 at = seq(0, 13),
						 colorkey = list(at = seq(1, 13),
										 labels = list(at = seq(1.5, 12.5),
													   labels = month_names)),
						 scales = list(x = list(at = seq(-180, 180, by = 20)),y = list(at = seq(-60, 90, by = 20))),
						 panel.grid = list(col = "lightgray", lty = "dotted"))+
			xyplot(lat ~ lon, wld, type='l', lty=1, lwd=1, col='black')
		png(filename=paste(FIG_DIR,"/SMVI_",DS,"_EMx_OnsetMon_",YYYY,".png",sep=''), width = 7000, height = 3500,res=300,type='cairo')
		print(Plt)
		dev.off()
		
		###E1 plot
		fd_YY <- as.data.frame(csv_in[,'fstdate1'],colClasses = c("date"))
		#fd_YY_Mon <- as.POSIXlt(fd_YY*86400, origin= "1970-01-01")$mon +1
		fd_YY_Mon <- month(fd_YY[,])
		month_matrix <- array(fd_YY_Mon, dim = c(nlon, nlat))


		var_MAT <- array(fd_YY, dim = c(nlon, nlat))    #Converts the vector to matrix of saze nlon x nlat
		month_names <- month.abb

		Plt<-levelplot(month_matrix,
						 row.values = lon,
						 column.values = lat,
						 xlab = 'Longitude',
						 ylab = 'Latitude',
						 ylim = c(-60, 90),
						 xlim = c(-180, 180),
						 main = paste(YYYY),
						 col.regions = color_palette[1:12],
						 at = seq(0, 13),
						 colorkey = list(at = seq(1, 13),
										 labels = list(at = seq(1.5, 12.5),
													   labels = month_names)),
						 scales = list(x = list(at = seq(-180, 180, by = 20)),y = list(at = seq(-60, 90, by = 20))),
						 panel.grid = list(col = "lightgray", lty = "dotted"))+
			xyplot(lat ~ lon, wld, type='l', lty=1, lwd=1, col='black')
		png(filename=paste(FIG_DIR,"/SMVI_",DS,"_E1_OnsetMon_",YYYY,".png",sep=''), width = 7000, height = 3500,res=300,type='cairo')
		print(Plt)
		dev.off()
		
		###E2 plot
		fd_YY <- as.data.frame(csv_in[,'fstdate2'],colClasses = c("date"))
		#fd_YY_Mon <- as.POSIXlt(fd_YY*86400, origin= "1970-01-01")$mon +1
		fd_YY_Mon <- month(fd_YY[,])
		month_matrix <- array(fd_YY_Mon, dim = c(nlon, nlat))


		var_MAT <- array(fd_YY, dim = c(nlon, nlat))    #Converts the vector to matrix of saze nlon x nlat
		month_names <- month.abb

		Plt<-levelplot(month_matrix,
						 row.values = lon,
						 column.values = lat,
						 xlab = 'Longitude',
						 ylab = 'Latitude',
						 ylim = c(-60, 90),
						 xlim = c(-180, 180),
						 main = paste(YYYY),
						 col.regions = color_palette[1:12],
						 at = seq(0, 13),
						 colorkey = list(at = seq(1, 13),
										 labels = list(at = seq(1.5, 12.5),
													   labels = month_names)),
						 scales = list(x = list(at = seq(-180, 180, by = 20)),y = list(at = seq(-60, 90, by = 20))),
						 panel.grid = list(col = "lightgray", lty = "dotted"))+
			xyplot(lat ~ lon, wld, type='l', lty=1, lwd=1, col='black')
		png(filename=paste(FIG_DIR,"/SMVI_",DS,"_E2_OnsetMon_",YYYY,".png",sep=''), width = 7000, height = 3500,res=300,type='cairo')
		print(Plt)
		dev.off()
	}
}
