#grass-tree-ca.R
#this script is used to compare model simulations and obs about tree and grass
#dynamics in CA, main focus is plant phys? 

#load all necessary packages for reading and plotting netCDF files

#library(readr)
library(tidyverse)   #just load tidyverse for dplyr, stringr, and tidyr all at once
#library(ggplot2)
library(ggmap)
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

flux_var <- read.csv("../data/AMF_US-Var_BASE_HH_16-5.csv", sep = ",", skip =2, #skip the first 2 rows
                    stringsAsFactors = FALSE, na.strings = '-9999') #missing value is indicated by -9999 
flux_ton <- read.csv("../data/AMF_US-Ton_BASE_HH_13-5.csv", sep = ",", skip =2, #skip the first 2 rows
                     stringsAsFactors = FALSE, na.strings = '-9999')

#read in biological measurements
bio_var <- read.csv("../data/AMF_US-Var_BIF_20210331.csv", sep = ",", 
                    stringsAsFactors = FALSE, na.strings = c('-9999', 'na', 'NA', ''))
bio_ton <- read.csv("../data/AMF_US-Ton_BIF_20210331.csv", sep = ",", 
                    stringsAsFactors = FALSE, na.strings = c('-9999', 'na', 'NA', ''))
bio2_long <- rbind(bio_var, bio_ton)
bio2_wide <- bio2_long %>% group_by(SITE_ID, GROUP_ID, VARIABLE_GROUP) %>% 
  pivot_wider(names_from = VARIABLE, values_from = DATAVALUE)
bio2_wide <- bio2_wide %>% select(c('SITE_ID', 'GROUP_ID', 'VARIABLE_GROUP',
                                    'AG_BIOMASS_GRASS', 'AG_BIOMASS_GRASS_SPATIAL_VARIABILITY',
                                    'AG_BIOMASS_GRASS_UNIT', 'AG_BIOMASS_DATE', 'LAI_DATE',
                                    'LAI_TOT', 'LAI_TOT_SPATIAL_VARIABILITY', 'SWC', 'SWC_DATE',
                                    'AG_BIOMASS_OTHER', 'AG_BIOMASS_OTHER_UNIT', 'AG_BIOMASS_TREE',
                                    'AG_BIOMASS_TREE_UNIT', 'AG_LIT_PROD_TOT', 'AG_LIT_PROD_UNIT',
                                    'LAI_U', 'LAI_U_SPATIAL_VARIABILITY'))
bio2_wide <- bio2_wide %>% mutate_at(c('AG_BIOMASS_GRASS', 'AG_BIOMASS_GRASS_SPATIAL_VARIABILITY',
                                       'LAI_TOT', 'LAI_TOT_SPATIAL_VARIABILITY', 'SWC',
                                       'AG_BIOMASS_OTHER', 'AG_BIOMASS_TREE', 'LAI_U', 'LAI_U_SPATIAL_VARIABILITY',
                                       'AG_LIT_PROD_TOT'),  list(~as.numeric(.)))
bio2_wide <- ungroup(bio2_wide)

# as obs are measured on diff date, might just have separate df for diff. avriable at each site 
# then combine
var_biomass <- bio2_wide %>% select('SITE_ID', 'AG_BIOMASS_DATE',
                                    'AG_BIOMASS_GRASS') %>% 
  mutate(DATE = ymd(AG_BIOMASS_DATE)) %>% filter(!is.na(AG_BIOMASS_GRASS))  %>% 
  select(-AG_BIOMASS_DATE)
var_biomass <- var_biomass %>%  mutate(year = year(DATE))
#data.table::setnames(var_biomass, old = "AG_BIOMASS_GRASS", new = "measurement")
#var_biomass$type <- "grass_ag"


var_lai <- bio2_wide %>% select('SITE_ID', 'LAI_DATE', 'LAI_TOT') %>% 
  mutate(DATE = ymd(LAI_DATE)) %>% filter(!is.na(LAI_TOT)) %>% select(-LAI_DATE)
var_lai <- var_lai %>% mutate(year = year(DATE), month = month(DATE)) %>% 
  group_by(year, month) %>% summarise(lai_var = mean(LAI_TOT, na.rm = TRUE))
#data.table::setnames(var_lai, old = 'LAI_TOT', new = 'measurement')
#var_lai$type <- 'lai_tot'

#tonzi site does not have good biomass measurements, ignore biomass
ton_lai <- bio2_wide %>% select('SITE_ID', 'LAI_DATE', 'LAI_U') %>% 
  mutate(DATE = ymd(LAI_DATE)) %>% filter(!is.na(LAI_U)) %>% select(-LAI_DATE)
ton_lai <- ton_lai %>% mutate(year = year(DATE), month = month(DATE)) %>% 
  group_by(year, month) %>% summarise(lai_ton = mean(LAI_U, na.rm = TRUE))


#laiobs <- rbind(var_lai, ton_lai)
#laiobs <- laiobs %>% mutate(date = as.yearmon(paste(year, month, sep="-")))
#laiobs <- laiobs %>% mutate(date = as.Date(date))
#data.table::setnames(ton_lai, old = 'LAI_U', new = 'measurement')
#ton_lai$type <- 'lai_u'
#bioall_long <- rbind(var_biomass, var_lai, ton_lai)
#bioall_wide <- pivot_wider(bioall_long, names_from = c(type, SITE_ID), values_from = measurement)

## select GPP, NEE, LE (latent heat), and H (sensible heat) and combine two flux datasets
var <- flux_var %>% select(TIMESTAMP_START, TIMESTAMP_END, LE_PI_F, H_PI_F, NDVI, NEE_PI_F,
                           GPP_PI_F)
##for ton data, *_1_1_1 measured at 23m height, *_1_2_1 measured at 2m
ton <- flux_ton %>% select(TIMESTAMP_START, TIMESTAMP_END, LE_PI_F_1_1_1, H_PI_F_1_1_1,
                           NDVI_1_1_1, NEE_PI_F_1_1_1, GPP_PI_F_1_1_1)
#tail(fluxobs$TIMESTAMP_END, n=1) #end date of obs
#head(fluxobs$TIMESTAMP_START, n=1) #start date of obs, hmm, it's numeric not time, need to convert
                                   # to time format so can calculate monthly mean for obs
                                    #flux obs are available from 2000/01/01 to 2020/12/31


#summary(ftest7)

#sink(file = "ncdinfo.txt")
#print(ftest7)
#sink() #save output from print() to txt file so can open to read info about netcdf file 
#read in the saved txt 
#ncdinfo <- read.delim("../ncdinfo.txt")
#varnames <- attributes(ftest7$var)$names #get all variable names in netCDF 
#dimnames <- attributes(ftest7$dim)$names

##get ecosystem production related variables, search for GPP, NPP, NEP, and NEE in var names

#varnames[grep('GPP', varnames)] #GPP, GPP_CANOPY, GPP_UNDERSTORY
#varnames[grep('NPP', varnames)] #NPP, NPP_CROOT, NPP_FROOT, NPP_SEED, NPP_STEM, NPP_STOR
#varnames[grep('NEP', varnames)] #only NEP
#varnames[grep('NEE', varnames)] #no NEE found

##latent heat (LH) and sensible heat (SH) flux variables
#varnames[grep('LH', varnames)] # EFLX_LH_TOT: total latent heat flux to atm, W/m^2
#varnames[grep('FSH', varnames)] #FSH? there are many sensible heat flux variables, difference???

## mortality variables
#varnames[grep('M', varnames)]  #M1-M10 (_SCLS (by size class) OR _SCPF (by pft/size class))
#M1: background mortality, M2: hydraulic mortality
#M3: carbon starvation mortality, M4: impact mortality
#M5: fire mortality, M6:termination mortality
#M7: logging mortality, M8: freezing mortality
#M9_SCLS: senescence mortality by size 
#M1_SCPF, M3_SCPF, M5_SCPF


#extract all variables from netcdf
gpp <- ncvar_get(ftest7, varid = 'GPP')
#dim(gpp) #dimension 14 22 1860 (lon, lat, time-step???)
#ncatt_get(ftest7, 'GPP') #unit: gC/m^2/s
npp <- ncvar_get(ftest7, varid = 'NPP')
#dim(npp)
#ncatt_get(ftest7, 'NPP') #same unit
nep <- ncvar_get(ftest7, varid = 'NEP')
#dim(nep)
#ncatt_get(ftest7, 'NEP') #same unit
lheat <- ncvar_get(ftest7, varid = 'EFLX_LH_TOT')  #W/m^2
sheat <- ncvar_get(ftest7, varid = "FSH")
lai <- ncvar_get(ftest7, varid = "ELAI")
#mortality
## mortality should be in %, so averaged over number of plants per ha by pft/size (NPLANT_SCPF)
bgloss <- ncvar_get(ftest7, varid = "M1_SCPF")
#dim(bgloss) # 4 dimensions: 14 (lat), 22 (lon), 39 (fates_levscpf), 1860 (month since 1860-01-01)
hydloss <- ncvar_get(ftest7, varid = "M2_SCPF")
carbloss <- ncvar_get(ftest7, varid = "M3_SCPF")
fireloss <- ncvar_get(ftest7, varid = "M5_SCPF")
termloss <- ncvar_get(ftest7, varid = "M6_SCPF") #when >2 canopy layers, understory veg will be terminated
nplants <- ncvar_get(ftest7, varid = "NPLANT_SCPF")



#get longitude and latitude
lon <- ncvar_get(ftest7, "lon"); dim(lon); ncatt_get(ftest7, 'lon') #dimension and attributes
lat <- ncvar_get(ftest7, "lat"); dim(lat); ncatt_get(ftest7, 'lat')
#ncatt_get(ftest7, 'time') #date since 1860-01-01
#tail(ncvar_get(ftest7, 'time'), n=1) #model simulation ends at 2014/12/31, 56575 (155yrs) days since 1860/01/01
                                     # as there are 1860 obs, so each ob is a monthly mean


#units for GPP, NEE and RECO from flux obs are umol CO2/m^2/s, need to convert to gC

#unit conversion from umolco2 to gC using bigleaf::umolCO2.to.gC/m^2/d
var <- var %>% mutate(nee = umolCO2.to.gC(NEE_PI_F),
                              gpp = umolCO2.to.gC(GPP_PI_F),
                      site = "var") %>% 
  select(-c(NEE_PI_F, GPP_PI_F))
ton <- ton %>% mutate(nee = umolCO2.to.gC(NEE_PI_F_1_1_1),
                      gpp = umolCO2.to.gC(GPP_PI_F_1_1_1),
                      site = "ton") %>% 
  select(-c(NEE_PI_F_1_1_1, GPP_PI_F_1_1_1))
#rename columns using data.table::setnames()
data.table::setnames(var, old = c("LE_PI_F", "H_PI_F", "NDVI"),
                     new = c("le", "h", "ndvi"))
data.table::setnames(ton, old = c("LE_PI_F_1_1_1", "H_PI_F_1_1_1",
                                  "NDVI_1_1_1"), new = c("le", 
                                                         "h", "ndvi"))
fluxobs <- rbind(var, ton)

#convert *_flux from gC/m^2/d to gC/m^2/s to be consistent with model simulation
fluxobs <- fluxobs %>% mutate(nee = nee/86400,
                              gpp = gpp/86400)
#convert time columns and subset 2000-2015 obs
fluxobs <- fluxobs %>% mutate(time_start = ymd_hm(TIMESTAMP_START, tz='UTC'),
                              time_end = ymd_hm(TIMESTAMP_END, tz = 'UTC')) %>%   
  select(-c(TIMESTAMP_START, TIMESTAMP_END))

#filter out obs between 2000-01-01 00:00 and 2016-01-01 00:00 to match model simulations
fluxobs$time_interval <- interval(fluxobs$time_start, fluxobs$time_end)
start <- ymd_hm(200001010000, tz = 'UTC')
end <- ymd_hm(201501010000, tz = 'UTC')
filter.interval <- interval(start, end)
fluxobs <- fluxobs %>% filter(time_interval%within%filter.interval)

#as model simulations are monthly means, so calculate monthly mean for flux obs
fluxobs <- fluxobs %>% mutate(year = year(time_start), 
                              month = month(time_start)) 
#monthly mean
month_flux <- fluxobs %>% group_by(year, month, site) %>% 
  summarise_at(c("nee", "gpp", "le", "h", "ndvi"),
               mean, na.rm = TRUE)

#get the GPP simulation closest to the flux tower site (38.4133, 239.0492),
# so the closest grid cell is lat[13]:38.25, and lon[9]: 239.25, and also 
# between 2000-01-01 and 2015-12-31

## store model outputs as time series 
modelsims <- month_flux %>% filter(site =="var") %>%  select(year, month)
modelsims$site <- "mod"


modelsims$gpp <- gpp[9,13, 1681:1860] #subset model simulations 2000-2014
modelsims$nee <- -1* nep[9,13, 1681:1860]
modelsims$le <- lheat[9, 13, 1681:1860]
modelsims$h <- sheat[9, 13, 1681:1860]
modelsims$lai <- lai[9, 13, 1681:1860]



#combine model simulated nee and gpp to flux observed ones
# both wide and long data formats for different plotting purposes
all_long <- rbind(month_flux, modelsims)
all_wide <- pivot_wider(all_long, id_cols = c(year, month),
                        names_from = site, 
                        values_from = c(nee, gpp, le, h, ndvi, lai)) 
all_wide <- select(all_wide, -c(lai_var, lai_ton))
#bind lai to all_wide
all_wide <- left_join(all_wide, var_lai, by = c("year", "month"))
all_wide <- left_join(all_wide, ton_lai, by = c("year", "month"))


# plot model simulation against flux observations for gpp and nee
#GPP
ggplot(all_wide, aes(x = gpp_mod)) + 
  geom_point(aes(y = gpp_var, color = "gpp_var")) +
  geom_point(aes(y = gpp_ton, color = "gpp_ton")) +
  geom_abline(slope =1, intercept = 0, size = 0.8) + 
  labs(x = expression(Model~simulated~GPP~(gC~m^-2~s^-1)),
       y = expression(Observed~GPP~(gC~m^-2~s^-1))) +
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  scale_color_manual(breaks=c("gpp_var", "gpp_ton"),
                     values = c("black", "red")) +
  prestheme.nogridlines   +
  theme(legend.title = element_blank(), 
        legend.position = "bottom")      #over estimated 

ggsave("../results/gpp-comparison.png", width = col2, height= 0.7*col2, 
units="cm", dpi = 800)

#NEE
ggplot(all_wide, aes(x = nee_mod)) + 
  geom_point(aes(y = nee_var, color = "nee_var")) +
  geom_point(aes(y = nee_ton, color = "nee_ton")) +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~NEE~(gC~m^-2~s^-1)),
       y = expression(Observed~NEE~(gC~m^-2~s^-1))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  scale_color_manual(breaks = c("nee_var", "nee_ton"),
                     values = c("black", "red")) +
  prestheme.nogridlines  +
  theme(legend.title = element_blank(), 
              legend.position = "bottom")#better predicted than GPP

ggsave("../results/nee-comparison.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)
#LH
ggplot(all_wide, aes(x = le_mod)) + 
  geom_point(aes(y = le_var, color = "le_var")) +
  geom_point(aes(y = le_ton, color = "le_ton")) +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~latent~heat~(W~m^-2)),
       y = expression(Observed~latent~heat~(W~m^-2))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  scale_color_manual(breaks = c("le_var", "le_ton"),
                     values = c("black", "red")) +
  prestheme.nogridlines  +
  theme(legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/le-comparison.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)


#SH
ggplot(all_wide, aes(x = h_mod)) + 
  geom_point(aes(y = h_var, color = "h_var")) +
  geom_point(aes(y = h_ton, color = "h_ton")) +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~sensible~heat~(W~m^-2)),
       y = expression(Observed~sensible~heat~(W~m^-2))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  scale_color_manual(breaks = c("h_var", "h_ton"),
                     values = c("black", "red")) +
  prestheme.nogridlines +
  theme(legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/sh-comparison.png", width = col2, height = 0.7*col2, 
       units = "cm", dpi = 800)

## LAI
ggplot(all_wide, aes(x = lai_mod)) + 
  geom_point(aes(y = lai_var, color = "lai_var")) +
  geom_point(aes(y = lai_ton, color = "lai_ton")) +
  geom_abline(slope =1, intercept = 0, size = 0.8) +
  labs(x = expression(Model~simulated~LAI~(m^2~m^-2)),
       y = expression(Observed~LAI~(m^2~m^-2))) + 
  scale_x_continuous(labels = comma_format( decimal.mark = ".")) +
  scale_y_continuous(labels = comma_format(decimal.mark = ".")) +
  scale_color_manual(breaks = c("lai_var", "lai_ton"),
                     values = c("black", "red")) +
  prestheme.nogridlines +
  theme(legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/lai-comparison.png", width = col2, height = 0.7*col2,
       units = "cm", dpi = 800)

## seasonal cycle of GPP and NEE for both simulations and obs
#first paste year + month to make a date column
all_long <- ungroup(all_long)
all_long <- all_long %>% mutate(date = (paste(all_long$year, 
                                                 all_long$month, sep = '')))
                                    
all_long$dates <- format(parse_date_time(all_long$date, orders = c("Y/m")), "%Y-%m") 
all_wide <- ungroup(all_wide)
all_wide <- all_wide %>% mutate(date = (paste(all_wide$year, 
                                              all_wide$month, sep = '')))
all_wide$dates <- format(parse_date_time(all_wide$date, orders = c("Y/m")), "%Y-%m")




#as ggplot only works with date class, so convert dates to class date
# with zoo::yearmon()

all_long $date <- as.yearmon(all_long$dates)
all_long$date <- as.Date(all_long$date)
all_long <- select(all_long, -dates)

all_wide $date <- as.yearmon(all_wide$dates)
all_wide$date <- as.Date(all_wide$date)
all_wide <- select(all_wide, -dates)

ggplot(all_long, aes(date, gpp, color = site)) + 
  geom_line() + 
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(GPP~(gC~m^-2~s^-1)))+
  scale_colour_manual(values = c("black", "red", "grey")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
  
ggsave("../results/seasonal-cycle-gpp.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

##time series plot for NEE

ggplot(all_long, aes(date, nee, color = site)) + 
  geom_line() +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(NEE~(gC~m^-2~s^-1)))+
  scale_colour_manual(values = c("black", "red", "grey")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")

ggsave("../results/seasonal-cycle-nee.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

# time series plot for latent heat flux

ggplot(all_long, aes(date, le, color = site)) + 
  geom_line()+
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(Latent~heat~flux~(W~m^-2)))+
  scale_colour_manual(values = c("black", "red", "grey")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-le.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

#SH

ggplot(all_long, aes(date, h, color = site)) + 
  geom_line() +
  scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(Sensible~heat~flux~(W~m^-2)))+
  scale_colour_manual(values = c("black", "red", "grey")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-h.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

##LAI cycle
ggplot(all_wide, aes(date, lai_mod, color = 'lai_mod')) +
  geom_line() +
  geom_line(aes(y = lai_ton, color = 'lai_ton')) +
  geom_line(aes(y = lai_var, color = 'lai_var')) +
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0)) +
  labs(x = 'Time', y = expression(LAI~(m^2~m^-2)))+
  scale_colour_manual(breaks = c('lai_mod', 'lai_ton', 'lai_var'),
    values = c("black", "red", "grey")) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-lai.png", width = 1.7*col2, 
       height = col2, units = 'cm', dpi =800)
  

######## what drives the seasonal cycle of GPP? ########

## contribution of background mortality (plant phenology), 
## carbon starvation mortality, and hydrology failure

## need to convert array of matrix to data frame (wide format)
## as the each row in the matrix is actually one variable and the columns
## are all obs for that variable, so need the transpose of the matrix first 
## so that each column is a variable then convert the matrix to 
## to data table (so in wide format), then combine all mortality to long data

##1:39 are pft (1: pine, 2: cedar, 3: grass)*size (1:13 classes, 0-65cm)
## 1:13~pine between 0-65cm; 14:26~cedar between 1-65; 27:39~grass between 1-65

bgloss_slice <- bgloss[9,13,,1681:1860]
hydloss_slice <- hydloss[9,13,,1681:1860]
carbloss_slice <- carbloss[9,13,,1681:1860]
fireloss_slice <- fireloss[9,13,,1681:1860]
termloss_slice <- termloss[9,13,,1681:1860]
nplants_slice <- nplants[9, 13,, 1681:1860]

name_list <- c(sprintf("pine_scl%s", seq(1:13)), sprintf("cedar_scl%s", seq(1:13)),
               sprintf("grass_scl%s", seq(1:13)))

 
nplants_slice <- as.data.frame(t(nplants_slice))
colnames(nplants_slice) <- name_list

bgloss_slice <- as.data.frame(t(bgloss_slice))
colnames(bgloss_slice) <- name_list
bgloss_slice <- bgloss_slice/(nplants_slice*86400*365)
bgloss_slice $type <- "bg"
bgloss_slice $year <- all_wide$year
bgloss_slice $month <- all_wide$month

hydloss_slice <- as.data.frame(t(hydloss_slice))
colnames(hydloss_slice) <- name_list
hydloss_slice <- hydloss_slice/(nplants_slice*86400*365)
hydloss_slice$type <- "hydro"
hydloss_slice$year <- all_wide$year
hydloss_slice$month <- all_wide$month

carbloss_slice <- as.data.frame(t(carbloss_slice))
colnames(carbloss_slice) <- name_list
carbloss_slice <- carbloss_slice/(nplants_slice*86400*365)
carbloss_slice$type <- "carb"
carbloss_slice$year <- all_wide$year
carbloss_slice$month <- all_wide$month

fireloss_slice <- as.data.frame(t(fireloss_slice))
colnames(fireloss_slice) <- name_list
fireloss_slice <- fireloss_slice/(nplants_slice*86400*365)
fireloss_slice$type <- "fire"
fireloss_slice$year <- all_wide$year
fireloss_slice$month <- all_wide$month

termloss_slice <- as.data.frame(t(termloss_slice))
colnames(termloss_slice) <- name_list
termloss_slice <- termloss_slice/(nplants_slice*86400*365)
termloss_slice$type <- "term"
termloss_slice$year <- all_wide$year
termloss_slice$month <- all_wide$month


#combine all and convert to long data format
death_wide <- rbind(bgloss_slice, hydloss_slice, carbloss_slice, fireloss_slice, termloss_slice)
death_long <- pivot_longer(death_wide, cols = -c(year, month, type),
                           names_to = "pft_scl",
                           values_to = "mortality" )

## the final big data set that includes everything in long format

big_long <- left_join(death_long, all_long, by = c("year", "month"))

## plot ggp seasonal cycle against different mortality sources to see 
## which drives the change in ggp
mod_big <- big_long %>% filter(site == "mod") %>% 
   mutate(pft = str_extract(pft_scl, "[^_]+")) #must be averaged across all size classes 
death_pfttype <- mod_big %>% group_by(date, pft, type) %>% 
  summarise(mortality = mean(mortality, na.rm = TRUE))
death_type <- mod_big %>% group_by(date, type) %>% 
  summarise(mortality = mean(mortality, na.rm=TRUE))
death_typewide <- pivot_wider(death_type, names_from = type, 
                              values_from = mortality)

### seasonal cycle of mortality by source
ggplot(death_type, aes(date, mortality, color = type)) +
  geom_line() +
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0))+ 
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Time', y = expression(Mortality~(Year^-1))) +
  scale_colour_manual(values = schwilkcolors) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")
ggsave("../results/seasonal-cycle-deathbytype.png", width = 1.7*col2,
       height = col2, units = 'cm', dpi = 800)

### seasonal cycle of mortality by pft and source
ggplot(death_pfttype, aes(date, mortality, color = type)) +
  geom_line() +
  facet_grid(pft~., scales = "free_y")+
  #scale_y_continuous(labels = comma_format(decimal.mark = "."))+
  scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               limits = as.Date(c('2000-01-01', '2014-12-31')),
               expand = expansion(0))+ 
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Time', y = expression(Mortality~(Year^-1)))+
  scale_colour_manual(values = schwilkcolors) +
  prestheme.nogridlines +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        legend.title = element_blank(), 
        legend.position = "bottom")

ggsave("../results/seasonal-cycle-deathby-pftype.png", width = 1.7*col2,
       height = 1.4*col2, units = 'cm', dpi = 800)
### it's mainly driven by grasses dying as both hydro failure and carbon starvation???

## plot hydro failure and carbon starvation mortality over gpp
#ggplot(all_wide, aes(date)) +
  #geom_line(aes(y = gpp_mod, color = "gpp_mod"))+
  #geom_line(aes(y = gpp_var, color = "gpp_var")) +
  #geom_line(aes(y = gpp_ton, color = "gpp_ton"))+
  #geom_line(aes(y = death_typewide$hydro, color = "hydro_fail")) +
  #geom_line(aes(y = death_typewide$carb, color = "carb_starv")) +
  #geom_line(aes(y = death_typewide$fire, color = "fire"))+
  #geom_line(aes(y = death_typewide$bg, color = "bg"))+
  #scale_x_date(date_labels= "%Y-%m", date_breaks = '3 months',
               #limits = as.Date(c('2000-01-01', '2014-12-31')),
               #expand = expansion(0)) +
#xlab("Time") + ylab("") + 
 #scale_colour_manual(breaks = c("gpp_mod", "hydro_fail", "carb_starv",
                                #"fire", "bg"), values = schwilkcolors) +
  #prestheme.nogridlines +
  #theme(axis.text.x=element_text(angle=60, hjust=1, size = pressmsz-6),
        #legend.title = element_blank(), 
        #legend.position = "bottom")

#ggsave("../results/all-gpp-death.png", width = 1.7*col2,
       #height = col2, units = 'cm', dpi = 800)

##both hydraulic failure and carbon starvation can be causes, but which one
## is more important?

## build linear mixed effects model to predict gpp
## independent variables: grass_hydro, grass_carb, cedar_hydro, cedar_carb, and
## pine_phen
death_pfttype_wide <- pivot_wider(death_pfttype, 
                                  names_from = c(pft, type),
                                  values_from = mortality)

model_data <- all_wide %>% left_join(death_pfttype_wide, by = "date") %>% 
  left_join(death_typewide, by = "date")

cordt <- model_data %>% select(grass_hydro, grass_carb, cedar_hydro,
                               cedar_carb, cedar_fire, pine_hydro, pine_carb, 
                               pine_fire,hydro, carb, fire)
corr <- cor(cordt)
ggcorrplot::ggcorrplot(corr)
ggsave("../results/correlation.png", dpi = 800)
pairs(cordt, pch = 19, lower.panel = NULL)

zscore <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm = TRUE) 

model_data <- model_data %>% mutate_at(c("grass_hydro", "grass_carb", "cedar_hydro", "cedar_fire",
                                         "cedar_carb", "pine_hydro", "pine_carb", "pine_fire",
                                         "carb", "hydro", "fire", "lai_mod"),
                               list(s = zscore))
model_data <- model_data %>% mutate(gpp_mod_log = log10(gpp_mod))


## plot mortaliry against gpp to explore the relationship then do model
ggplot(model_data, aes(hydro, gpp_mod)) +
  geom_point() +
  scale_y_continuous("Log (GPP)",
    trans = "log10") +
  xlab(expression(Hydraulic~mortality~(Year^-1))) +
  prestheme.nogridlines
ggsave("../results/gpp-phen.png", dpi = 800, width = col2, 
       height = 0.8*col2, units = "cm") #b/c phen is negatively correlated to carb?

## also plot mortality source by pft to each mortality source to see how each pft
## contributes to each mortality


## model
gpp.lm1 <- lm(gpp_mod_log ~ carb_s + hydro_s+fire_s + lai_mod_s, data = model_data)
summary(gpp.lm1)
plot(gpp.lm1)
car::Anova(gpp.lm1)

gpp.lm2 <- lm(gpp_mod_log ~ grass_carb_s + grass_hydro_s + cedar_hydro_s +
                cedar_fire_s + pine_fire_s + pine_carb_s + lai_mod_s, 
              data = model_data)
summary(gpp.lm2)
plot(gpp.lm2)
car::Anova(gpp.lm2)

##close connection to netcdf 
nc_close(ftest7)

## plot function from 
## https://sites.google.com/view/climate-access-cooperative/code?authuser=0

mapCDFtemp <- function(lat,lon,tas){ #model and perc should be a string
  
  titletext <- "title"
  expand.grid(lon, lat) %>%
    rename(lon = Var1, lat = Var2) %>%
    mutate(lon = ifelse(lon > 180, -(360 - lon), lon),
           tas = as.vector(tas)) %>% 
    
    #mutate(tas = convert_temperature(tas, "k", "c")) %>%
    
    ggplot() + 
    geom_point(aes(x = lon, y = lat, color = tas),
               size = 0.8) + 
    borders(database="state", regions = "california", colour="black", fill=NA) + 
    scale_color_viridis(na.value="white",name = "Temperature") + 
    theme(legend.direction="vertical", legend.position="right", legend.key.width=unit(0.4,"cm"), legend.key.heigh=unit(2,"cm")) + 
    coord_quickmap() + 
    ggtitle(titletext) 
}


#clean env.
rm('fluxobs', 'ecoprod', 'gpp', 'nep', 'npp', 'varnames', 
   'start', 'end', 'filter.interval', 'nep_slice', 'bgloss', 'hydloss',
   'carbloss', 'fireloss', 'termloss')

