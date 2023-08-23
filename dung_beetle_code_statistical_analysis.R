#Install packages
install.packages("metafor")
install.packages("tidyverse")
install.packages("Rtools")
install.packages("devtools")
devtools::install_github("daniel1noble/orchaRd", force = TRUE)
install.packages("ggplot2")
install.packages("maps")
install.packages("mapdata")

#Load packages
library(metafor)
library(tidyverse)
library(orchaRd)
library(cowplot)
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)

#Set working directory
setwd()

#Read in the data
study<-read.csv("Beetle_study.csv")
sites<-read.csv("Beetle_sites.csv")
outcomes<-read.csv("Beetle_outcomes.csv")

###################################################
#1 - Data cleaning and imputation #################
###################################################

#Cleaning data
#Remove extra rows
beetle_study<-study%>%
  filter(Person_ID=="Niamh")

beetle_sites<-sites%>%
  filter(Person_ID=="Niamh")

beetle_outcomes<-outcomes%>%
  filter(Person_ID=="Niamh")


#Remove extra columns
col_details<-data.frame(col_name=names(beetle_sites),
                        col_index=seq(1,122))

beetle_sites_clean<-dplyr::select(beetle_sites,-c(32:122))

#There were no extra columns to clean in 'Beetle_study' or 'Beetle_outcomes'

#Combine the data from the different datasets
beetle_all_data<-beetle_outcomes%>%
  left_join(beetle_sites_clean,"Site_ID")%>%
  left_join(beetle_study,"Study_ID")

#calculate the standard deviation
beetle_all_data <- beetle_all_data %>%
  mutate(
    #convert SE to SD
    control_SD=ifelse(var_type=="SE", Control_var_combined*sqrt(control_n), control_var),
    treatment_SD=ifelse(var_type=="SE",Treatment_var_combined*sqrt(treatment_n), treatment_var),
    #Estimate SD
    #Calculate the coefficient of variation (SD/mean)
    control_CV=control_SD/Control_mean_combined,
    treatment_CV=treatment_SD/Treatment_mean_combined,
    #Use the coefficient of variation to estimate the SD for data points where it is missing
    control_SD=ifelse(is.na(control_SD),Control_mean_combined*median(control_CV,na.rm=TRUE),control_SD),
    treatment_SD=ifelse(is.na(treatment_SD),Treatment_mean_combined*median(treatment_CV,na.rm=TRUE),treatment_SD))

#interpolate sample sizes
beetle_all_data <- beetle_all_data %>%
  mutate(
    #replace -999 with NA
    control_n=ifelse(control_n==-999,NA,control_n),
    treatment_n=ifelse(treatment_n==-999,NA,treatment_n),
    #replace missing sample sizes with the median sample sizes for each group
    control_n=ifelse((control_n),median(control_n,na.rm=TRUE),control_n),
    treatment_n=ifelse(is.na(treatment_n),median(treatment_n,na.rm=TRUE),treatment_n))

#convert time since burning to years
beetle_all_data$last_burned_treatment_years<-beetle_all_data$last_burned_treatment/365

#Calculate effect sizes
#Standardized mean difference
dung_beetle_smd<-escalc(data = beetle_all_data,m2i=Control_mean_combined,m1i = Treatment_mean_combined,
                        sd2i=control_SD,sd1i = treatment_SD,
                        n2i=control_n,n1i=treatment_n,measure = "SMD")


#Removing data with zeros in both control and treatment groups
dung_beetle_smd$control_plus_treatment<-dung_beetle_smd$control_mean+dung_beetle_smd$treatment_mean
dung_beetle_smd_filtered<-filter(dung_beetle_smd,control_plus_treatment>0)

#Subset data for each biodiversity metric

#abundance
dung_beetle_smd_abundance<-dung_beetle_smd_filtered%>%
  filter(detailed_outcome=="density of individuals")

#biomass
dung_beetle_smd_biomass<-dung_beetle_smd_filtered%>%
  filter(detailed_outcome=="biomass")

#species richness
dung_beetle_smd_richness<-dung_beetle_smd_filtered%>%
  filter(detailed_outcome=="species richness")

#Shannon diversity
dung_beetle_smd_shannon<-dung_beetle_smd_filtered%>%
  filter(detailed_outcome=="Shannon Index")

#Simpson's diversity
dung_beetle_smd_simpson<-dung_beetle_smd_filtered%>%
  filter(detailed_outcome=="Simpson's Diversity Index")


###################################################
#2 - Analysis #####################################
###################################################

### Overall effect of fire ###

## Abundance ##
#remove rows with effect size equal to na
no_nas_abundance<-filter(dung_beetle_smd_abundance,!is.na(yi))

#overall effect of fire
#the 'overall effect' is the weighted mean of the comparison between control and treatment (unburnt vs burnt)
#yi is the effect size.  Mean of the treatment group – mean of the control group/pooled SD
#vi is the variance of the effect size.  This is used to weight the analysis.  Datapoints that are less variable have more weight in the analysis as more precise estimates are deemed more trustworthy.
null_model_abundance<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#overall effect on each habitat type
habitat_model_abundance<-rma.mv(yi,vi,mods=~Habitat_type_combined2-1,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#calculate R^2 for habitat type
1-(deviance(habitat_model_abundance)/deviance(null_model_abundance))

## Biomass ##
#remove rows with effect size equal to na
no_nas_biomass<-filter(dung_beetle_smd_biomass,!is.na(yi))

#overall effect of fire
null_model_biomass<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_biomass)

## Richness ##
#remove rows with effect size equal to na
no_nas_richness<-filter(dung_beetle_smd_richness,!is.na(yi))

#overall effect of fire
null_model_richness<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#overall effect on each habitat type
habitat_model_richness<-rma.mv(yi,vi,mods=~Habitat_type_combined2-1,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

## Shannon diversity ##
#remove rows with effect size equal to na
no_nas_shannon<-filter(dung_beetle_smd_shannon,!is.na(yi))

#overall effect of fire
null_model_shannon<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_shannon)

## Simpson's diversity ##

#remove rows with effect size equal to na
no_nas_simpson<-filter(dung_beetle_smd_simpson,!is.na(yi))

#overall effect of fire
null_model_simpson<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_simpson)


### Effect of severity ###

## Abundance ##
#overall effect of aboveground severity
above_severity_model_abundance<-rma.mv(yi,vi,mods=~ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#the warning message relates to the fact that there are 146 rows that report abundance but do not have an above ground severity index.  This is also true for below ground severity.
sum(is.na(dung_beetle_smd_abundance$ag_severity_index))

#overall effect of belowground severity
below_severity_model_abundance<-rma.mv(yi,vi,mods=~bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

sum(is.na(dung_beetle_smd_abundance$bg_severity_index))

##calculate R^2 for above and below ground severity
dung_beetle_smd_abundance_ag_severity <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(ag_severity_index))

dung_beetle_smd_abundance_bg_severity <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(bg_severity_index))

null_severity_model_abundance_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_ag_severity)
null_severity_model_abundance_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_bg_severity)

1-(deviance(above_severity_model_abundance)/deviance(null_severity_model_abundance_ag))
1-(deviance(below_severity_model_abundance)/deviance(null_severity_model_abundance_bg))


## Richness ##
#overall effect of aboveground severity
above_severity_model_richness<-rma.mv(yi,vi,mods=~ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#the warning message relates to the fact that there are 10 rows that report richness but do not have an above ground severity index.  This is also true for below ground severity.
sum(is.na(dung_beetle_smd_richness$ag_severity_index))

#overall effect of belowground severity
below_severity_model_richness<-rma.mv(yi,vi,mods=~bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

sum(is.na(dung_beetle_smd_richness$bg_severity_index))

#calculate R^2 for above and below ground severity
dung_beetle_smd_richness_ag_severity <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(ag_severity_index))

dung_beetle_smd_richness_bg_severity <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(bg_severity_index))

null_severity_model_richness_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_ag_severity)
null_severity_model_richness_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_bg_severity)

1-(deviance(above_severity_model_richness)/deviance(null_severity_model_richness_ag))
1-(deviance(below_severity_model_richness)/deviance(null_severity_model_richness_bg))


### Recovery times ###

## Abundance ##
#overall effect of time since burning
time_model_abundance<-rma.mv(yi,vi,mods=~last_burned_treatment_years,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#effect of aboveground severity on the influence of time since burning
time_above_severity_model_abundance<-rma.mv(yi,vi,mods=~last_burned_treatment_years*ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#effect of belowground severity on the influence of time since burning
time_below_severity_model_abundance<-rma.mv(yi,vi,mods=~last_burned_treatment_years*bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#calculate R^2

#time since last fire
dung_beetle_smd_abundance_recovery <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(last_burned_treatment_years))

null_recovery_model_abundance<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_recovery)

1-(deviance(time_model_abundance)/deviance(null_recovery_model_abundance))

#influence of above ground severity on the impact of time since fire
dung_beetle_smd_abundance_recovery_ag <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(last_burned_treatment_years) & !is.na(ag_severity_index))

null_recovery_model_abundance_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_recovery_ag)

1-(deviance(time_above_severity_model_abundance)/deviance(null_recovery_model_abundance_ag))

#influence of below ground severity on the impact of time since fire
dung_beetle_smd_abundance_recovery_bg <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(last_burned_treatment_years) & !is.na(bg_severity_index))

null_recovery_model_abundance_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_recovery_bg)

1-(deviance(time_below_severity_model_abundance)/deviance(null_recovery_model_abundance_bg))

## Sensitivity analysis for abundance (removing studies that have data for over 20 years since last fire)

# Filter data to exclude rows with 'last_burned_treatment_years' over 20
sensitivity_dung_beetle_smd_abundance <- dung_beetle_smd_abundance %>%
  filter(last_burned_treatment_years <= 20)

# Fit the model with filtered data
sensitivity_time_model_abundance <- rma.mv(yi, vi, mods = ~last_burned_treatment_years,
                               random = list(~1|Study_ID/Site_ID, ~1|Site_name),
                               data = sensitivity_dung_beetle_smd_abundance)


## Richness ##
#overall effect of time since burning
time_model_richness<-rma.mv(yi,vi,mods=~last_burned_treatment_years,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#effect of aboveground severity on the influence of time since burning
time_above_severity_model_richness<-rma.mv(yi,vi,mods=~last_burned_treatment_years*ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#effect of belowground severity on the influence of time since burning
time_below_severity_model_richness<-rma.mv(yi,vi,mods=~last_burned_treatment_years*bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#calculate R^2

#time since last fire
dung_beetle_smd_richness_recovery <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(last_burned_treatment_years))

null_recovery_model_richness<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_recovery)

1-(deviance(time_model_richness)/deviance(null_recovery_model_richness))

#influence of above ground severity on the impact of time since fire
dung_beetle_smd_richness_recovery_ag <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(last_burned_treatment_years) & !is.na(ag_severity_index))

null_recovery_model_richness_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_recovery_ag)

1-(deviance(time_above_severity_model_richness)/deviance(null_recovery_model_richness_ag))

#influence of below ground severity on the impact of time since fire
dung_beetle_smd_richness_recovery_bg <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(last_burned_treatment_years) & !is.na(bg_severity_index))

null_recovery_model_richness_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_recovery_bg)

1-(deviance(time_below_severity_model_richness)/deviance(null_recovery_model_richness_bg))


### Influence of historical fire interval ###

## Abundance ##
#overall effect of historical fire interval
historical_fire_abundance<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#historical fires and aboveground severity
historical_fire_ag_severity_abundance<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome*ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#historical fires and belowground severity
historical_fire_bg_severity_abundance<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome*bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#calculate R^2

#historical fire return interval
dung_beetle_smd_abundance_historical <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(historical_fire_frequency_biome))

null_model_abundance_historical<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_historical)

1-(deviance(historical_fire_abundance)/deviance(null_model_abundance_historical))

#influence of above ground severity on the impact of historical fire return interval
dung_beetle_smd_abundance_historical_ag <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(historical_fire_frequency_biome)& !is.na(ag_severity_index))

null_model_abundance_historical_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_historical)

1-(deviance(historical_fire_abundance)/deviance(null_model_abundance_historical_ag))

#influence of below ground severity on the impact of historical fire return interval
dung_beetle_smd_abundance_historical_bg <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "density of individuals" & !is.na(historical_fire_frequency_biome)& !is.na(bg_severity_index))

null_model_abundance_historical_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance_historical)

1-(deviance(historical_fire_bg_severity_abundance)/deviance(null_model_abundance_historical_bg))

## Richness ##
#overall effect of historical fire interval
historical_fire_richness<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#historical fires and aboveground severity
historical_fire_ag_severity_richness<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome*ag_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#historical fires and belowground severity
historical_fire_bg_severity_richness<-rma.mv(yi,vi,mods=~historical_fire_frequency_biome*bg_severity_index,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#calculate R^2

# historical fire return interval
dung_beetle_smd_richness_historical <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(historical_fire_frequency_biome))

null_model_richness_historical<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_historical)

1-(deviance(historical_fire_richness)/deviance(null_model_richness_historical))

#influence of above ground severity on the impact of historical fire return interval
dung_beetle_smd_richness_historical_ag <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(historical_fire_frequency_biome)& !is.na(ag_severity_index))

null_model_richness_historical_ag<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_historical)

1-(deviance(historical_fire_richness)/deviance(null_model_richness_historical_ag))

#influence of below ground severity on the impact of historical fire return interval
dung_beetle_smd_abundance_richness_bg <- dung_beetle_smd_filtered %>%
  filter(detailed_outcome == "species richness" & !is.na(historical_fire_frequency_biome)& !is.na(bg_severity_index))

null_model_richness_historical_bg<-rma.mv(yi,vi,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness_historical)

1-(deviance(historical_fire_bg_severity_richness)/deviance(null_model_richness_historical_bg))


#####################################################
#3 - Effect of studies with low internal validity ###
#####################################################

#Running overall models with only high validity studies
null_model_abundance_validity <- rma.mv(yi, vi, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_abundance, Internal_validity == 'HIGH'))
null_model_biomass_validity <- rma.mv(yi, vi, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_biomass, Internal_validity == 'HIGH'))
null_model_richness_validity <- rma.mv(yi, vi, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_richness, Internal_validity == 'HIGH'))

above_severity_model_abundance_validity <- rma.mv(yi, vi, mods = ~ag_severity_index, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_abundance, Internal_validity == 'HIGH'))
above_severity_model_richness_validity <- rma.mv(yi, vi, mods = ~ag_severity_index, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_richness, Internal_validity == 'HIGH'))

below_severity_model_abundance_validity <- rma.mv(yi, vi, mods = ~bg_severity_index, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_abundance, Internal_validity == 'HIGH'))
below_severity_model_richness_validity <- rma.mv(yi, vi, mods = ~bg_severity_index, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_richness, Internal_validity == 'HIGH'))

time_model_abundance_validity <- rma.mv(yi, vi, mods = ~last_burned_treatment_years, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_abundance, Internal_validity == 'HIGH'))
time_model_richness_validity <- rma.mv(yi, vi, mods = ~last_burned_treatment_years, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_richness, Internal_validity == 'HIGH'))

historical_fire_abundance_validity <- rma.mv(yi, vi, mods = ~historical_fire_frequency_biome, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_abundance, Internal_validity == 'HIGH'))
historical_fire_richness_validity <- rma.mv(yi, vi, mods = ~historical_fire_frequency_biome, random = list(~1|Study_ID/Site_ID, ~1|Site_name), data = subset(dung_beetle_smd_richness, Internal_validity == 'HIGH'))


###################################################
#4 FIGURES ########################################
###################################################

##########
#FIGURE 1. MAP OF STUDY LOCATIONS
##########

#count number of sites within 1 degree radius of each other
#filter data to remove duplicate coordinates
sites_count<-
  beetle_sites_clean%>%
  mutate(lat_1_deg=round(Lat_dec_deg/4)*4,lon_1_deg=round(Lon_dec_deg/4)*4)%>%
  group_by(lat_1_deg,lon_1_deg,Habitat_type_combined2)%>%
  summarise(no_sites=length(lat_1_deg),mean_lat=mean(Lat_dec_deg),mean_lon=mean(Lon_dec_deg))

#Load world map data
world <- map_data("world")

#Create the base map
p <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray40", size = 0.1) +
  coord_quickmap()

subset(beetle_sites_clean,Habitat_type_combined2=="")

# Create the plot
p <- ggplot(data = sites_count, aes(x = lon_1_deg, y = lat_1_deg)) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray40", size = 0.1) +
  geom_point(aes(color = Habitat_type_combined2, size = no_sites), alpha = 0.5) +
  scale_color_viridis_d(name = "", guide = guide_legend(override.aes = list(size = 5))) +
  labs(x = "", y = "", title = "",
       size = "No. of sites") +  # Set size key title
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),  
        legend.box.background = element_blank(),  
        legend.margin = margin(0))  

# Display the plot
print(p)

##############
#FIGURE 2. OVERALL EFFECTS OF FIRE
##############

null_model_results_abundance<-mod_results(null_model_abundance,mod="1",at=NULL,group="Study_ID")
null_model_results_richness<-mod_results(null_model_richness,mod="1",at=NULL,group="Study_ID")
null_model_results_biomass<-mod_results(null_model_biomass,mod="1",at=NULL,group="Study_ID")
null_model_results_shannon<-mod_results(null_model_shannon,mod="1",at=NULL,group="Study_ID")
null_model_results_simpson<-mod_results(null_model_simpson,mod="1",at=NULL,group="Study_ID")

null_model_list<-list(null_model_abundance,null_model_richness,null_model_biomass,
                      null_model_shannon,null_model_simpson)
outcome_list<-c("Abundance","Species richness","Biomass","Shannon Index","Simpson's Index")
null_predictions<-NULL
for(i in 1:length((null_model_list))){
  beta<-null_model_list[[i]]$beta
  LCI<-null_model_list[[i]]$ci.lb
  UCI<-null_model_list[[i]]$ci.ub
  p_val<-null_model_list[[i]]$pval
  k<-null_model_list[[i]]$k.all
  null_predictions_temp<-data.frame(beta,LCI,UCI,p_val,k,outcome=outcome_list[i])
  null_predictions<-rbind(null_predictions,null_predictions_temp)
}

#format names of outcomes for figure
dung_beetle_smd_filtered<-dung_beetle_smd_filtered%>%
  mutate(new_outcome = fct_recode(detailed_outcome, "Species richness" = "species richness",
                                  "Abundance" = "density of individuals",
                                  "Biomass" = "biomass",
                                  "Simpson's Index" = "Simpson's Diversity Index",
                                  "Shannon Index" = "Shannon Index"))

# Set colors for each biodiversity metric
metric_colours <- c("Abundance" = "blue", "Species richness" = "green", "Biomass" = "red", "Shannon Index" = "turquoise", "Simpson's Index" = "orange")

#plot
ggplot(null_predictions, aes(x = beta, xmin = LCI, xmax = UCI, y = outcome)) +
  geom_point(size = 3) +
  geom_errorbarh() +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Standardised mean difference") +
  ylab("")+
  geom_beeswarm(
    data = dung_beetle_smd_filtered,
    aes(x = yi, y = new_outcome, colour = new_outcome),
    alpha = 0.2
  ) +
  scale_color_manual(values = metric_colours, name = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#Plot with re-ordered metrics
plot <- ggplot(null_predictions, aes(x = beta, xmin = LCI, xmax = UCI, y = outcome)) +
  geom_point(size = 3) +
  geom_errorbarh() +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Standardised mean difference") +
  ylab("") +
  geom_beeswarm(
    data = dung_beetle_smd_filtered,
    aes(x = yi, y = new_outcome, colour = new_outcome),
    alpha = 0.2
  ) +
  scale_color_manual(values = metric_colours, name = "", guide = "none") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Reorder the y-axis labels
desired_order <- c("Simpson's Index", "Shannon Index", "Species richness", "Biomass", "Abundance")
plot <- plot + scale_y_discrete(limits = desired_order)

# Display the plot
print(plot)


################
#FIGURE 3. EFFECT OF FIRE IN EACH HABITAT
################

#ABUNDANCE 
habitat_model_results_abundance <- mod_results(habitat_model_abundance, mod = "Habitat_type_combined2", at = NULL, group = "Study_ID")
abundance_habitat_plot <- orchard_plot(habitat_model_results_abundance, xlab = "Standardised mean difference")
abundance_habitat_plot <- abundance_habitat_plot + ggtitle("A")
abundance_habitat_plot <- abundance_habitat_plot + theme(panel.grid = element_blank())
print(abundance_habitat_plot)

#RICHNESS
habitat_model_results_richness <- mod_results(habitat_model_richness, mod = "Habitat_type_combined2", at = NULL, group = "Study_ID")
richness_habitat_plot <- orchard_plot(habitat_model_results_richness, xlab = "Standardised mean difference")
richness_habitat_plot <- richness_habitat_plot + ggtitle("C")
richness_habitat_plot <- richness_habitat_plot + theme(panel.grid = element_blank())
print(richness_habitat_plot)

### GRID ###
create_plot <- function(habitat_model_results, title) {
  plot <- orchard_plot(habitat_model_results, xlab = "Standardised mean difference", legend.pos = "none")
  plot <- plot + ggtitle(title)
  plot <- plot + theme(panel.grid = element_blank(),
                       text = element_text(size = 6),  
                       axis.text.y = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 0.5),  # Adjust y-axis text size, angle, and position
                       axis.text.x = element_text(size = 8),  
                       axis.title.y = element_text(size = 10),  
                       axis.title.x = element_text(size = 10),  
                       axis.ticks.x = element_line(size = 0.2))  
  return(plot)
}

# ABUNDANCE
habitat_model_results_abundance <- mod_results(habitat_model_abundance, mod = "Habitat_type_combined2", at = NULL, group = "Study_ID")
abundance_habitat_plot <- create_plot(habitat_model_results_abundance, "")

# RICHNESS
habitat_model_results_richness <- mod_results(habitat_model_richness, mod = "Habitat_type_combined2", at = NULL, group = "Study_ID")
richness_habitat_plot <- create_plot(habitat_model_results_richness, "")

# Arrange plots in a grid
grid <- plot_grid(abundance_habitat_plot, richness_habitat_plot,
                  ncol = 2, align = "hv", nrow = 1, labels = c("A", "B"))

# Display the grid
print(grid)


############
#FIGURE 4. EFFECT OF SEVERITY
############

### ag abundance
ag_results_abundance<-mod_results(above_severity_model_abundance,mod="ag_severity_index",at=NULL,group="Study_ID",weights = "prop")
ag_results_abundance_preds<-ag_results_abundance$mod_table

#subset data so that only rows with ag severity are present
filtered_ag_abundance<-filter(no_nas_abundance,!is.na(ag_severity_index))

abundance_predictions_ag<-data.frame(filtered_ag_abundance,predict(above_severity_model_abundance))      
ag_results_abundance$mod_table

# Plot the effect of aboveground severity on abundance
ggplot(ag_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) + 
  geom_jitter(data = filtered_ag_abundance, aes(ag_severity_index, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Above ground severity index", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = ""))

### ag richness
ag_results_richness<-mod_results(above_severity_model_richness,mod="ag_severity_index",at=NULL,group="Study_ID",weights = "prop")
ag_results_richness_preds<-ag_results_richness$mod_table

#subset data so that only rows with ag severity are present
filtered_ag_richness<-filter(no_nas_richness,!is.na(ag_severity_index))

richness_predictions_ag<-data.frame(filtered_ag_richness,predict(above_severity_model_richness))      
ag_results_richness$mod_table

# Plot the effect of aboveground severity on richness
ggplot(ag_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +
  geom_jitter(data = filtered_ag_richness, aes(ag_severity_index, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Above ground severity index", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("green")) +
  guides(color = guide_legend(title = ""))


### bg abundance
bg_results_abundance<-mod_results(below_severity_model_abundance,mod="bg_severity_index",at=NULL,group="Study_ID",weights = "prop")
bg_results_abundance_preds<-bg_results_abundance$mod_table

#subset data so that only rows with ag severity are present
filtered_bg_abundance<-filter(no_nas_abundance,!is.na(bg_severity_index))

abundance_predictions_bg<-data.frame(filtered_bg_abundance,predict(below_severity_model_abundance))      
bg_results_abundance$mod_table

# Plot the effect of belowground severity on abundance
ggplot(bg_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +  
  geom_jitter(data = filtered_bg_abundance, aes(bg_severity_index, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Below ground severity index", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = ""))

### bg richness
bg_results_richness<-mod_results(below_severity_model_richness,mod="bg_severity_index",at=NULL,group="Study_ID",weights = "prop")
bg_results_richness_preds<-bg_results_richness$mod_table

#subset data so that only rows with ag severity are present
filtered_bg_richness<-filter(no_nas_richness,!is.na(bg_severity_index))

richness_predictions_bg<-data.frame(filtered_bg_richness,predict(below_severity_model_richness))      
bg_results_richness$mod_table

#Plot the effect of belowground severity on richness
ggplot(bg_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +  
  geom_jitter(data = filtered_bg_richness, aes(bg_severity_index, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Below ground severity index", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("green")) +
  guides(color = guide_legend(title = ""))

### GRID ###
plot_abundance_ag <- ggplot(ag_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +  
  geom_jitter(data = filtered_ag_abundance, aes(ag_severity_index, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Above ground severity index (%)", y = "Standardised mean difference", title = "")+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12),  axis.title.y = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = ""))

plot_richness_ag <- ggplot(ag_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +  
  geom_jitter(data = filtered_ag_richness, aes(ag_severity_index, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Above ground severity index (%)", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12),  axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12)) +
  scale_color_manual(values = c("green")) +
  guides(color = guide_legend(title = ""))

plot_abundance_bg <- ggplot(bg_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "lightblue", alpha = 0.5) +  
  geom_jitter(data = filtered_bg_abundance, aes(bg_severity_index, yi), color = "blue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Below ground severity index (%)", y = "Standardised mean difference", title = "") +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("lightblue")) +
  guides(color = guide_legend(title = ""))

plot_richness_bg <-  ggplot(bg_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "lightgreen", alpha = 0.5) +  
  geom_jitter(data = filtered_bg_richness, aes(bg_severity_index, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Below ground severity index (%)", y = "Standardised mean difference", title = "") +
  theme(axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12),axis.title.x = element_text(size = 12)) +
  scale_color_manual(values = c("lightgreen")) +
  guides(color = guide_legend(title = ""))

# Arrange plots in a grid
grid <- plot_grid(
  plot_abundance_ag,plot_abundance_bg, plot_richness_ag,
  plot_richness_bg,
  nrow = 2, labels = c("A", "B", "C", "D")
)

# Display the grid
print(grid)


###############
#FIGURE 5. EFFECT OF RECOVERY TIME
###############

#abundance
recovery_results_abundance<-mod_results(time_model_abundance,mod="last_burned_treatment_years",at=NULL,group="Study_ID",weights = "prop")
recovery_results_abundance_preds<-recovery_results_abundance$mod_table

#subset data so that only rows with historical fire interval are present
filtered_recovery_abundance<-filter(no_nas_abundance,!is.na(last_burned_treatment_years))

abundance_predictions_recovery<-data.frame(filtered_recovery_abundance,predict(time_model_abundance))      
recovery_results_abundance$mod_table

# Plot the effect of recovery on abundance
ggplot(recovery_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +  # Pale blue color for the confidence interval ribbon
  geom_jitter(data = filtered_recovery_abundance, aes(last_burned_treatment_years, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Time since last fire", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

### richness
recovery_results_richness<-mod_results(time_model_richness,mod="last_burned_treatment_years",at=NULL,group="Study_ID",weights = "prop")
recovery_results_richness_preds<-recovery_results_richness$mod_table

#subset data so that only rows with historical fire interval are present
filtered_recovery_richness<-filter(no_nas_richness,!is.na(last_burned_treatment_years))

richness_predictions_recovery<-data.frame(filtered_recovery_richness,predict(time_model_richness))      
recovery_results_richness$mod_table

#Plot the effect of recovery time on richness
ggplot(recovery_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +  # Pale blue color for the confidence interval ribbon
  geom_jitter(data = filtered_recovery_richness, aes(last_burned_treatment_years, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Time since last fire", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("green")) +
  guides(color = guide_legend(title = "Biome"))

### GRID ###

# Abundance plot
plot_abundance <- ggplot(recovery_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +
  geom_jitter(data = filtered_recovery_abundance, aes(last_burned_treatment_years, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Time since last fire (years)", y = "Standardised mean difference", title = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12),  axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

# Richness plot
plot_richness <- ggplot(recovery_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +
  geom_jitter(data = filtered_recovery_richness, aes(last_burned_treatment_years, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Time since last fire (years)") +
  theme(axis.title.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12),axis.title.x = element_text(size = 10)) +
  scale_color_manual(values = c("green")) +
  guides(color = guide_legend(title = "Biome"))

# Arrange plots in a grid
grid <- plot_grid(
  plot_abundance,plot_richness,
  nrow = 1, align = "hv", labels = c("A", "B")
)

# Display the grid
print(grid)


################
#Figure 6. Effect of severity on recovery time
################

## ag abundance
ag_severity_recovery_ab<-rma.mv(yi,vi,mods=~ag_severity_index*last_burned_treatment_years,
                                random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

ag_severity_recovery_ab

#copy model formula
ag_sev_rec_ab_formula<-(~ag_severity_index*last_burned_treatment_years-1)

#create dataframe with new data for predictions
#lines of best fit for the impact of time since last fire at different severity levels
new_data_ab<-data.frame(expand.grid(
  ag_severity_index = c(25,50,75,100),
  last_burned_treatment_years=seq(0,20,0.1)))

#create a model matrix and remove the intercept
predgrid_ab<-model.matrix(ag_sev_rec_ab_formula,data=new_data_ab)

#predict onto the new model matrix
mypreds_ab<-data.frame(predict.rma(ag_severity_recovery_ab,newmods=predgrid_ab))

#attach predictions to variables for plotting
new_data_ab <- cbind(new_data_ab, mypreds_ab[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#### ag richness ###
ag_severity_recovery_richness<-rma.mv(yi,vi,mods=~ag_severity_index*last_burned_treatment_years,
                                      random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#copy model formula
ag_sev_rec_richness_formula<-(~ag_severity_index*last_burned_treatment_years-1)

#create dataframe with new data for predictions
new_data_richness<-data.frame(expand.grid(
  ag_severity_index = c(25,50,75,100),
  last_burned_treatment_years=seq(0,20,0.1)))

#create a model matrix and remove the intercept
predgrid_richness<-model.matrix(ag_sev_rec_ab_formula,data=new_data_richness)

#predict onto the new model matrix
mypreds_richness<-data.frame(predict.rma(ag_severity_recovery_richness,newmods=predgrid_richness))

#attach predictions to variables for plotting
new_data_richness <- cbind(new_data_richness, mypreds_richness[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

### bg abundance ###
bg_severity_recovery_ab<-rma.mv(yi,vi,mods=~bg_severity_index*last_burned_treatment_years,
                                random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)


#copy model formula
bg_sev_rec_ab_formula<-(~bg_severity_index*last_burned_treatment_years-1)

#create dataframe with new data for predictions
new_data_ab_bg<-data.frame(expand.grid(
  bg_severity_index = c(25,50,75,100),
  last_burned_treatment_years=seq(0,20,0.1)))


#create a model matrix and remove the intercept
predgrid_ab_bg<-model.matrix(bg_sev_rec_ab_formula,data=new_data_ab_bg)

#predict onto the new model matrix
mypreds_ab_bg<-data.frame(predict.rma(bg_severity_recovery_ab,newmods=predgrid_ab_bg))

#attach predictions to variables for plotting
new_data_ab_bg <- cbind(new_data_ab_bg, mypreds_ab_bg[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

### bg richness ###
bg_severity_recovery_richness<-rma.mv(yi,vi,mods=~bg_severity_index*last_burned_treatment_years,
                                      random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)


#copy model formula
bg_sev_rec_richness_formula<-(~bg_severity_index*last_burned_treatment_years-1)

#create dataframe with new data for predictions
new_data_richness_bg<-data.frame(expand.grid(
  bg_severity_index = c(25,50,75,100),
  last_burned_treatment_years=seq(0,20,0.1)))

#create a model matrix and remove the intercept
predgrid_richness_bg<-model.matrix(bg_sev_rec_richness_formula,data=new_data_richness_bg)

#predict onto the new model matrix
mypreds_richness_bg<-data.frame(predict.rma(bg_severity_recovery_richness,newmods=predgrid_richness_bg))

#attach predictions to variables for plotting
new_data_richness_bg <- cbind(new_data_richness_bg, mypreds_richness_bg[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

## 4x4 GRID

##ag abundance
# Define custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common ggplot object with custom colors
plot_ag_abundance_4x4_recovery <- ggplot(new_data_ab, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +   
  theme_cowplot() +
  labs(x = "",
       y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)")+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank()) 


# Separate plots using facet_wrap
plot_ag_abundance_4x4_recovery +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4)  


## ag richness
# Define custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common ggplot object with custom colors
plot_ag_richness_4x4_recovery <- ggplot(new_data_richness, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +   
  theme_cowplot() +
  labs(x = "",
       y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)")+
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank()) 


# Separate plots using facet_wrap
plot_ag_richness_4x4_recovery +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4, strip.position = "bottom") +  
  theme(strip.text.x = element_blank())  


## bg abundance

# Define custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common ggplot object with custom colors
plot_bg_abundance_4x4_recovery <- ggplot(new_data_ab_bg, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +   
  theme_cowplot() +
  labs(x = "",
       y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)")+
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank()) 


# Separate plots using facet_wrap
plot_bg_abundance_4x4_recovery +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +  
  theme(strip.text.x = element_blank())  


## bg richness

# Define custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common ggplot object with custom colors
plot_bg_richness_4x4_recovery <- ggplot(new_data_richness_bg, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +  
  scale_fill_manual(values = custom_colors) +   
  theme_cowplot() +
  labs(x = "Time since last fire (years)",
       y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)")+
  theme(legend.position = "none") 


# Separate plots using facet_wrap
plot_bg_richness_4x4_recovery +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +  
  theme(strip.text.x = element_blank())  


# STACK

# Define the custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common y-axis title
y_axis_title <- "Standardised mean difference"

theme_mod <- theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        axis.title.y.left = element_text(vjust = 2))

# Create the ggplot objects for each plot
plot_ag_abundance_4x4_recovery <- ggplot(new_data_ab, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4)

plot_ag_richness_4x4_recovery <- ggplot(new_data_richness, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_abundance_4x4_recovery <- ggplot(new_data_ab_bg, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_richness_4x4_recovery <- ggplot(new_data_richness_bg, aes(x = last_burned_treatment_years, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "Time since last fire (years)",
       y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none") +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())


# Arrange the plots in a 4x4 grid
recovery_grid<-plot_grid(plot_ag_abundance_4x4_recovery+theme(axis.title.y=element_blank()),
                         plot_ag_richness_4x4_recovery+theme(axis.title.y=element_blank()),
                         plot_bg_abundance_4x4_recovery+theme(axis.title.y=element_blank()),
                         plot_bg_richness_4x4_recovery+theme(axis.title.y=element_blank()),
                         ncol = 1, align = "v",rel_heights = c(1, 1, 1, 1),
                         labels = c("A","B","C","D"))


#create common y labels

y.grob <- textGrob("Standardised mean difference", 
                   gp=gpar(col="black", fontsize=12), rot=90)


recovery_grid_with_y_axis<-grid.arrange(arrangeGrob(recovery_grid, left = y.grob))


##############
#FIGURE 7. EFFECT OF HISTORICAL FIRE INTERVAL
##############

### abundance
historical_results_abundance<-mod_results(historical_fire_abundance,mod="historical_fire_frequency_biome",at=NULL,group="Study_ID",weights = "prop")
historical_results_abundance_preds<-historical_results_abundance$mod_table

#subset data so that only rows with historical fire interval are present
filtered_historical_abundance<-filter(no_nas_abundance,!is.na(historical_fire_frequency_biome))

abundance_predictions_historical<-data.frame(filtered_historical_abundance,predict(historical_fire_abundance))      
historical_results_abundance$mod_table

# Plot the effect of historical fire interval on abundance
ggplot(historical_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +  # Pale blue color for the confidence interval ribbon
  geom_jitter(data = filtered_historical_abundance, aes(historical_fire_frequency_biome, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Historical fire interval", y = "Standardised mean difference", title = "Abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

#richness
historical_results_richness<-mod_results(historical_fire_richness,mod="historical_fire_frequency_biome",at=NULL,group="Study_ID",weights = "prop")
historical_results_richness_preds<-historical_results_richness$mod_table

#subset data so that only rows with historical fire interval are present
filtered_historical_richness<-filter(no_nas_richness,!is.na(historical_fire_frequency_biome))

richness_predictions_historical<-data.frame(filtered_historical_richness,predict(historical_fire_richness))      
historical_results_richness$mod_table

# Plot the effect of historical fire interval on richness
ggplot(historical_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +  # Pale blue color for the confidence interval ribbon
  geom_jitter(data = filtered_historical_richness, aes(historical_fire_frequency_biome, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Historical fire interval", y = "Standardised mean difference", title = "Richness") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), text = element_text(size = 12)) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

### GRID ###
title_size <- 6

plot_abundance <- ggplot(historical_results_abundance_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "blue", alpha = 0.5) +
  geom_jitter(data = filtered_historical_abundance, aes(historical_fire_frequency_biome, yi), color = "darkblue", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Historical fire return interval (years)", y = "Standardised mean difference", title = "") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = title_size)
  ) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

plot_richness <- ggplot(historical_results_richness_preds, aes(moderator, estimate, ymin = lowerCL, ymax = upperCL)) +
  geom_line() +
  geom_ribbon(fill = "green", alpha = 0.5) +  # Pale blue color for the confidence interval ribbon
  geom_jitter(data = filtered_historical_richness, aes(historical_fire_frequency_biome, yi), color = "darkgreen", alpha = 0.5, inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Historical fire return interval (years)", y = "Standardised mean difference", title = "") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = title_size)
  ) +
  scale_color_manual(values = c("blue")) +
  guides(color = guide_legend(title = "Biome"))

# Arrange the plots into a grid
grid <- plot_grid(plot_abundance, plot_richness, ncol = 2, nrow = 1, labels = c("A", "B"))

# Display the grid
grid

#########
#FIGURE 8. EFFECT OF SEVERITY ON THE INFLUENCE OF HISTORICAL FIRE FREQUENCY 
#########

#### ag Abundance
ag_severity_historical_abundance<-rma.mv(yi,vi,mods=~ag_severity_index*historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#copy model formula
ag_sev_historical_abundance_formula<-(~ag_severity_index*historical_fire_frequency_biome-1)

#create dataframe with new data for predictions
#lines of best fit for the impact of historical fire return interval at different severity levels
new_data_abundance_historical<-data.frame(expand.grid(
  ag_severity_index = c(25,50,75,100),
  historical_fire_frequency_biome=seq(12,291,12.1)))

#create a model matrix and remove the intercept
predgrid_abundance_historical<-model.matrix(ag_sev_historical_abundance_formula,data=new_data_abundance_historical)

#predict onto the new model matrix
mypreds_abundance_historical<-data.frame(predict.rma(ag_severity_historical_abundance,newmods=predgrid_abundance_historical))

#attach predictions to variables for plotting
new_data_abundance_historical <- cbind(new_data_abundance_historical, mypreds_abundance_historical[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

##### ag richness
ag_severity_historical_richness<-rma.mv(yi,vi,mods=~ag_severity_index*historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#copy model formula
ag_sev_historical_richness_formula<-(~ag_severity_index*historical_fire_frequency_biome-1)

#create dataframe with new data for predictions
new_data_richness_historical<-data.frame(expand.grid(
  ag_severity_index = c(25,50,75,100),
  historical_fire_frequency_biome=seq(12,291,12.1)))

#create a model matrix and remove the intercept
predgrid_richness_historical<-model.matrix(ag_sev_historical_richness_formula,data=new_data_richness_historical)

#predict onto the new model matrix
mypreds_richness_historical<-data.frame(predict.rma(ag_severity_historical_richness,newmods=predgrid_richness_historical))

#attach predictions to variables for plotting
new_data_richness_historical <- cbind(new_data_richness_historical, mypreds_richness_historical[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

###### bg abundance
bg_severity_historical_abundance<-rma.mv(yi,vi,mods=~bg_severity_index*historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_abundance)

#copy model formula
bg_sev_historical_abundance_formula<-(~bg_severity_index*historical_fire_frequency_biome-1)

#create dataframe with new data for predictions
new_data_abundance_historical_bg<-data.frame(expand.grid(
  bg_severity_index = c(25,50,75,100),
  historical_fire_frequency_biome=seq(12,291,12.1)))

#create a model matrix and remove the intercept
predgrid_abundance_historical_bg<-model.matrix(bg_sev_historical_abundance_formula,data=new_data_abundance_historical_bg)

#predict onto the new model matrix
mypreds_abundance_historical_bg<-data.frame(predict.rma(bg_severity_historical_abundance,newmods=predgrid_abundance_historical_bg))

#attach predictions to variables for plotting
new_data_abundance_historical_bg <- cbind(new_data_abundance_historical_bg, mypreds_abundance_historical_bg[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

##### bg richness ###
bg_severity_historical_richness<-rma.mv(yi,vi,mods=~bg_severity_index*historical_fire_frequency_biome,random=list(~1|Study_ID/Site_ID,~1|Site_name),data=dung_beetle_smd_richness)

#copy model formula
bg_sev_historical_richness_formula<-(~bg_severity_index*historical_fire_frequency_biome-1)

#create dataframe with new data for predictions
new_data_richness_historical_bg<-data.frame(expand.grid(
  bg_severity_index = c(25,50,75,100),
  historical_fire_frequency_biome=seq(12,291,12.1)))

#create a model matrix and remove the intercept
predgrid_richness_historical_bg<-model.matrix(bg_sev_historical_richness_formula,data=new_data_richness_historical_bg)

#predict onto the new model matrix
mypreds_richness_historical_bg<-data.frame(predict.rma(bg_severity_historical_richness,newmods=predgrid_richness_historical_bg))

#attach predictions to variables for plotting
new_data_richness_historical_bg <- cbind(new_data_richness_historical_bg, mypreds_richness_historical_bg[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

# 4x4 grid

# STACK

# Define the custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common y-axis title
y_axis_title <- "Standardised mean difference"

theme_mod <- theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        axis.title.y.left = element_text(vjust = 2))

# Create the ggplot objects for each plot
plot_ag_abundance_4x4 <- ggplot(new_data_abundance_historical, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4)

plot_ag_richness_4x4 <- ggplot(new_data_richness_historical, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_abundance_4x4 <- ggplot(new_data_abundance_historical_bg, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_richness_4x4 <- ggplot(new_data_richness_historical_bg, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "Historical fire return interval (years)",
       y = "Standardised mean difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none") +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())


# Arrange the plots in a 4x4 grid
historical_grid<-plot_grid(plot_ag_abundance_4x4+theme(axis.title.y=element_blank()),
                         plot_ag_richness_4x4+theme(axis.title.y=element_blank()),
                         plot_bg_abundance_4x4+theme(axis.title.y=element_blank()),
                         plot_bg_richness_4x4+theme(axis.title.y=element_blank()),
                         ncol = 1, align = "v",rel_heights = c(1, 1, 1, 1),
                         labels = c("A","B","C","D"))


#create common y labels

y.grob <- textGrob("Standardised mean difference", 
                   gp=gpar(col="black", fontsize=12), rot=90)


historical_grid_with_y_axis<-grid.arrange(arrangeGrob(historical_grid, left = y.grob))






# Define the custom colors
custom_colors <- c("#FF5733", "#FFC300", "#85C1E9", "#27AE60")

# Create a common y-axis title
y_axis_title <- "Standardised mean difference"

theme_mod <- theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        axis.title.y.left = element_text(vjust = 2))

# Create the ggplot objects for each plot
plot_ag_abundance_4x4 <- ggplot(new_data_abundance_historical, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised Mean Difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4)

plot_ag_richness_4x4 <- ggplot(new_data_richness_historical, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(ag_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(ag_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised Mean Difference",
       fill = "Above Ground Severity Index",  
       colour = "Above Ground Severity Index") +  
  scale_fill_discrete(name = "Above Ground Severity Index (%)") +  
  scale_color_discrete(name = "Above Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(ag_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_abundance_4x4 <- ggplot(new_data_abundance_historical_bg, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "", y = "Standardised Mean Difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())

plot_bg_richness_4x4 <- ggplot(new_data_richness_historical_bg, aes(x = historical_fire_frequency_biome, y = pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_line(aes(colour = as.factor(bg_severity_index))) +
  geom_ribbon(alpha = 0.1, aes(fill = as.factor(bg_severity_index))) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 7),  
        axis.title.y = element_text(size = 8))+
  labs(x = "Historical fire interval (years)",
       y = "Standardised Mean Difference",
       fill = "Below Ground Severity Index",  
       colour = "Below Ground Severity Index") +  
  scale_fill_discrete(name = "Below Ground Severity Index (%)") +  
  scale_color_discrete(name = "Below Ground Severity Index (%)") +
  theme(legend.position = "none") +
  facet_wrap(~ as.factor(bg_severity_index), ncol = 4, strip.position = "bottom") +
  theme(strip.text.x = element_blank())


# Arrange the plots in a 4x4 grid
plot_grid(plot_ag_abundance_4x4, plot_ag_richness_4x4, plot_bg_abundance_4x4, plot_bg_richness_4x4,
          ncol = 1, align = "v",rel_heights = c(1, 1, 1, 1))+
  ylab(y_axis_title) 
