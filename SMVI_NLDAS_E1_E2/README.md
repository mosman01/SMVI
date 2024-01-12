# Soil Moisture Volatility Index (SMVI)
## SMVI - NLDAS-2
### Corresponding Author: Mahmoud Osman
#### Email-1: mahosman01@gmail.com
#### Emain-2: mahmoud.osman@jhu.edu

Last Update: December 2023

#########################################

SMVI data is in this inventory is based on the published studies of flash droughts over the contiguous United States (Osman et al. 2021 & 2022).


Based on the SMVI methodology, I identify up to 6 flash drought events per year (Mar – Nov), which never actually happens, and the longest event is considered E01, and second is E02. Based on my later post processing, most of the events are E01 and count of E02 are nearly 30% of what is detected in E01.

My recommendation is to go with E01, and you can consider E02 as another set of flash drought events of shorter duration.

#### Regarding the fields names, this is the description:

- Lon : 		Longitude of the gridpoint
- Lat : 		Latitude  of the gridpoint
- Date : 		Identified flash drought event onset date
- Length: 	Length of the event
- FD : 		Flag (0,1) to represent flash drought or not

#### The follwoing are the 15-days average (3-pentads) value for selected variables' anomaly 3-pentad (15-days) back from Onset date (i.e. average from Onset to Onset-15d).

- TMP : 		 Temperature 
- PRCP : 		Same as TMP, for Precipitation
- RZSM : 		Same as TMP, for root-zone soil moisture
- EVP : 		Same as TMP, for evapotranspiration
- PEVP : 		Same as TMP, for potential evapotranspiration
- TCDC : 		Same as TMP, for total cloud cover
- PRESS : 	Same as TMP, for surface pressure
- WS : 		Same as TMP, for 10-m above ground wind speed
- CAPE : 		Same as TMP, for convective available potential energy
- VPD : 		Same as TMP, for vapor pressure deficit

- VEG_mask :	Flag (0,1) to represent regions of short vegetation (e.g. crops and shrubs)
- clust_3 :	Value (1,2,3) represents the clustered group for the detected event based on the classification process.


### References:
-----------
- Osman, M., Zaitchik, B. F., Badr, H. S., Christian, J. I., Tadesse, T., Otkin, J. A., & Anderson, M. C. (2021). Flash drought onset over the contiguous United States: sensitivity of inventories and trends to quantitative definitions. Hydrology and Earth System Sciences, 25(2), 565–581. https://doi.org/10.5194/hess-25-565-2021. 

Osman, M., Zaitchik, B. F., Badr, H. S., Otkin, J., Zhong, Y., Lorenz, D., et al. (2022). Diagnostic Classification of Flash Drought Events Reveals Distinct Classes of Forcings and Impacts. Journal of Hydrometeorology, 23(2), 275–289. https://doi.org/10.1175/JHM-D-21-0134.1.
