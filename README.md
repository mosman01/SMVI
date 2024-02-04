# Soil Moisture Volatility Index (SMVI)
### Corresponding Author: Mahmoud Osman
#### Email-1: mahosman01@gmail.com
#### Emain-2: mahmoud.osman@jhu.edu

Last Update: December 2023

#########################################

SMVI data is in this inventory is based on the published studies of flash droughts over the contiguous United States (Osman et al. 2021 & 2022) with slight modifations and enhancements.

## A flash drought is idenified in a gridpoint when the following conditions are fulfilled:

		1- The 5-day running average (1-Pentad) RZSM fall below the 20-days running average RZSM.
		2- The crossover occurs below the 20th percentile RZSM of the corresponding day's long-term record.
		3- The previous conditions are persistent for at least 4-pentads period.
		4- The grid-point's average temperature during the detected period is more than Zero C (273 K).
		5- The grid-point's computed average Bowen's ratio during the detected period is between 0.2 and 7.

	- In this version multiple flash drought events are allowed to be detected up to 6 events per year, sorted by onset date of event
	- Events 15 days or less apart are considered to be 1 event, and merged into one.

## Global SMVI inventory:
A global SMVI dataset is published in hydroshare at: https://www.hydroshare.org/resource/080002bd7cc44242bb37c02b049ed532/

## References:
- Osman, M., Zaitchik, B. F., Badr, H. S., Christian, J. I., Tadesse, T., Otkin, J. A., & Anderson, M. C. (2021). Flash drought onset over the contiguous United States: sensitivity of inventories and trends to quantitative definitions. Hydrology and Earth System Sciences, 25(2), 565–581. https://doi.org/10.5194/hess-25-565-2021
- Osman, M., Zaitchik, B. F., Badr, H. S., Otkin, J., Zhong, Y., Lorenz, D., et al. (2022). Diagnostic Classification of Flash Drought Events Reveals Distinct Classes of Forcings and Impacts. Journal of Hydrometeorology, 23(2), 275–289. https://doi.org/10.1175/JHM-D-21-0134.1
- Osman, M., Zaitchik, B. F., Otkin, J., Anderson, M. (2024). SMVI Global Flash Droughts Dataset, HydroShare, https://doi.org/10.4211/hs.080002bd7cc44242bb37c02b049ed532

