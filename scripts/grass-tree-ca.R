#grass-tree-ca.R
#this script is used to compare model simulations and obs about tree and grass
#dynamics in CA, main focus is plant phys? 

#load all necessary packages for reading and plotting netCDF files

#library(readr)
library(tidyverse)   #just load tidyverse for dplyr, stringr, and tidyr all at once
library(ggplot2)
#library(ggmap)
#library(viridis)
#library(weathermetrics)
library(ncdf4) 
#library(chron)
library(RColorBrewer)
#library(lattice)
library(bigleaf)    #package for processing flux data
library(lubridate)  #deal with time series 
library(zoo)

#load ggplot theme
source("./ggplot-theme.R")

#read in netcdf file 
ftest7 <- nc_open("../data/fates_clm50_cali_test7_4pfts_test035_bb419bcc_e6cb1a72.fullrun.nc")

#read in flux data, tower site lat: 38.4133 degrees_north, lon: 239.0492 degrees_east

fluxobs <- read.csv("../data/AMF_US-Var_BASE_HH_16-5.csv", sep = ",", skip =2, #skip the first 2 rows
                    stringsAsFactors = FALSE, na.strings = '-9999') #missing value is indicated by -9999 
colnames(fluxobs) #variable names in flux obs, might just focus on NEE & GPP for now
tail(fluxobs$TIMESTAMP_END, n=1) #end date of obs
head(fluxobs$TIMESTAMP_START, n=1) #start date of obs, hmm, it's numeric not time, need to convert
                                   # to time format so can calculate monthly mean for obs
                                    #flux obs are available from 2000/01/01 to 2020/12/31


#summary(ftest7)

#sink(file = "ncdinfo.txt")
#print(ftest7)
#sink() #save output from print() to txt file so can open to read info about netcdf file 
#read in the saved txt 
#ncdinfo <- read.delim("../ncdinfo.txt")
varnames <- attributes(ftest7$var)$names #get all variable names in netCDF 

#get ecosystem production related variables, search for GPP, NPP, NEP, and NEE in var names

varnames[grep('GPP', varnames)] #GPP, GPP_CANOPY, GPP_UNDERSTORY
varnames[grep('NPP', varnames)] #NPP, NPP_CROOT, NPP_FROOT, NPP_SEED, NPP_STEM, NPP_STOR
varnames[grep('NEP', varnames)] #only NEP
varnames[grep('NEE', varnames)] #no NEE found

#latent heat (LH) and sensible heat (SH) flux variables
varnames[grep('LH', varnames)] # EFLX_LH_TOT: total latent heat flux to atm, W/m^2
varnames[grep('FSH', varnames)] #FSH? there are many sensible heat flux variables, difference???


#extract all variables from netcdf
gpp <- ncvar_get(ftest7, varid = 'GPP')
#dim(gpp) #dimension 14 22 1860 (lon, lat, time-step???)
ncatt_get(ftest7, 'GPP') #unit: gC/m^2/s
npp <- ncvar_get(ftest7, varid = 'NPP')
#dim(npp)
ncatt_get(ftest7, 'NPP') #same unit
nep <- ncvar_get(ftest7, varid = 'NEP')
#dim(nep)
ncatt_get(ftest7, 'NEP') #same unit
lheat <- ncvar_get(ftest7, varid = 'EFLX_LH_TOT')  #W/m^2
sheat <- ncvar_get(ftest7, varid = "FSH")
#get longitude and latitude
lon <- ncvar_get(ftest7, "lon"); dim(lon); ncatt_get(ftest7, 'lon') #dimension and attributes
lat <- ncvar_get(ftest7, "lat"); dim(lat); ncatt_get(ftest7, 'lat')
ncatt_get(ftest7, 'time') #date since 1860-01-01
tail(ncvar_get(ftest7, 'time'), n=1) #model simulation ends at 2014/12/31, 56575 (155yrs) days since 1860/01/01
                                     # as there are 1860 obs, so each ob is a monthly mean


#filter out GPP, NEE and RECO (ecosystem respiration) from flux obs,
#units for GPP, NEE and RECO are umol CO2/m^2/s, need to convert to gC

ecoprod <- fluxobs %>% select(TIMESTAMP_START, TIMESTAMP_END, NEE_PI_F, GPP_PI_F, RECO_PI_F, 
                              LE, H)

#unit conversion from umolco2 to gC using bigleaf::umolCO2.to.gC/m^2/d
ecoprod <- ecoprod %>% mutate(nee_flux = umolCO2.to.gC(NEE_PI_F),
                              gpp_flux = umolCO2.to.gC(GPP_PI_F),
                              reco_flux = umolCO2.to.gC(RECO_PI_F)) %>% 
  select(-c(NEE_PI_F, GPP_PI_F, RECO_PI_F))
#rename latent heat and sensible heat flux
ecoprod <- rename(ecoprod, sh_flux = H)
ecoprod <- rename(ecoprod, lh_flux = LE)

#convert *_flux from gC/m^2/d to gC/m^2/s to be consistent with model simulation
ecoprod <- ecoprod %>% mutate(nee_flux = nee_flux/86400,
                              gpp_flux = gpp_flux/86400,
                              reco_flux = reco_flux/86400)

#convert time columns and subset 2000-2015 obs
ecoprod <- ecoprod %>% mutate(time_start = ymd_hm(TIMESTAMP_START, tz='UTC'),
                              time_end = ymd_hm(TIMESTAMP_END, tz = 'UTC')) %>%   
  select(-c(TIMESTAMP_START, TIMESTAMP_END))

#filter out obs between 2000-01-01 00:00 and 2016-01-01 00:00 to match model simulations
ecoprod$time_interval <- interval(ecoprod$time_start, ecoprod$time_end)
start <- ymd_hm(200001010000, tz = 'UTC')
end <- ymd_hm(201501010000, tz = 'UTC')
filter.interval <- interval(start, end)
ecoprod <- ecoprod %>% filter(time_interval%within%filter.interval)

#as model simulations are monthly means, so calculate monthly mean for flux obs
ecoprod <- ecoprod %>% mutate(year = year(time_start), 
                              month = month(time_start)) 
#monthly mean
month_flux <- ecoprod %>% group_by(year, month) %>% 
  summarise_at(c("nee_flux", "gpp_flux", "reco_flux", "lh_flux", "sh_flux"),
               mean, na.rm = TRUE)

#get the GPP simulation closest to the flux tower site (38.4133, 239.0492),
# so the closest grid cell is lat[13]:38.25, and lon[9]: 239.25, and also 
# between 2000-01-01 and 2015-12-31
gpp_slice <- gpp[9,13, 1681:1860] #subset model simulations 
nep_slice <- nep[9,13, 1681:1860]
lh_slice <- lheat[9, 13, 1681:1860]
sh_slice <- sheat[9, 13, 1681:1860]

#convert nep to nee
nee_slice <- -nep_slice

#combine model simulated nee and gpp to flux observed ones
month_flux$gpp_mod <- gpp_slice
month_flux$nee_mod <- nee_slice
month_flux$lh_mod <- lh_slice
month_flux$sh_mod <- sh_slice

# plot model simulation against flux observations for gpp and nee
#GPP
ggplot(month_flux, aes(gpp_mod, gpp_flux)) + geom_point() +
  geom_abline(slope =1, intercept = 0, size = 0.8) + 
  labs(x = expression(Model~simulated~GPP~(gC~m^-2~s^-1)),
       y = expression(Observed~GPP~(gC~m^-2~s^-1))) +
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  prestheme.nogridlines              #over estimated 

ggsave("../results/gpp.png", width = col2, height= 0.7*col2, 
units="cm", dpi = 800)

#NEE
ggplot(month_flux, aes(nee_mod, nee_flux)) + geom_point() +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~NEE~(gC~m^-2~s^-1)),
       y = expression(Observed~NEE~(gC~m^-2~s^-1))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  prestheme.nogridlines             #better predicted than GPP

ggsave("../results/nee.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)
#LH
ggplot(month_flux, aes(lh_mod, lh_flux)) + geom_point() +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~latent~heat~(W~m^-2)),
       y = expression(Observed~latent~heat~(W~m^-2))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  prestheme.nogridlines  
ggsave("../results/lh.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)

#SH
ggplot(month_flux, aes(sh_mod, sh_flux)) + geom_point() +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~sensible~heat~(W~m^-2)),
       y = expression(Observed~sensible~heat~(W~m^-2))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  prestheme.nogridlines
ggsave("../results/sh.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)
## seasonal cycle of GPP and NEE for both simulations and obs
#first paste year + month to make a date column
month_flux <- ungroup(month_flux)
month_flux <- month_flux %>% mutate(date = (paste(month_flux$year, 
                                                 month_flux$month, sep = '')))
                                    
month_flux$dates <- format(parse_date_time(month_flux$date, orders = c("Y/m")), "%Y-%m") 


#as ggplot only works with date class, so convert dates to class date
# with zoo::yearmon()

month_flux $date <- as.yearmon(month_flux$dates)
month_flux$date <- as.Date(month_flux$date)
month_flux <- select(month_flux, -dates)

ggplot(month_flux, aes(x = date)) + 
  geom_line(aes(y = gpp_mod,color = "gpp_mod")) +
  geom_line(aes(y = gpp_flux, color = "gpp_flux")) +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(GPP~(gC~m^-2~s^-1)))+
  scale_colour_manual(breaks=c("gpp_mod", "gpp_flux"),
                      values = c("black", "red")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
  
ggsave("../results/seasonal-cycle-gpp.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

##time series plot for NEE

ggplot(month_flux, aes(x = date)) + 
  geom_line(aes(y = nee_mod,color = "nee_mod")) +
  geom_line(aes(y = nee_flux, color = "nee_flux")) +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(NEE~(gC~m^-2~s^-1)))+
  scale_colour_manual(breaks=c("nee_mod", "nee_flux"),
                      values = c("black", "red")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")

ggsave("../results/seasonal-cycle-nee.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

# time series plot for latent heat flux

ggplot(month_flux, aes(x = date)) + 
  geom_line(aes(y = lh_mod,color = "lh_mod")) +
  geom_line(aes(y = lh_flux, color = "lh_flux")) +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(Latent~heat~flux~(W~m^-2)))+
  scale_colour_manual(breaks=c("lh_mod", "lh_flux"),
                      values = c("black", "red")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-lh.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

#SH

ggplot(month_flux, aes(x = date)) + 
  geom_line(aes(y = sh_mod,color = "sh_mod")) +
  geom_line(aes(y = sh_flux, color = "sh_flux")) +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(Sensible~heat~flux~(W~m^-2)))+
  scale_colour_manual(breaks=c("sh_mod", "sh_flux"),
                      values = c("black", "red")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-sh.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

#close connection to netcdf 
nc_close(ftest7)

#clean env.
rm('fluxobs', 'ecoprod', 'gpp', 'nep', 'npp', 'varnames', 
   'start', 'end', 'filter.interval', 'nep_slice')

