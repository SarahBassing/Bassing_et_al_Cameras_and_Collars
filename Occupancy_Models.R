  #'  R code associated with: 
  #'  
  #'  Bassing, S. B., M. Devivo, T. R. Ganz, B. N. Kertson, L. R. Prugh, T. Roussin, 
  #'  L. Satterfield, R. M. Windell, A. J. Wirsing, & B. Gardner 
  #'  
  #'  Are we telling the same story? Comparing inferences made from camera trap 
  #'  and telemetry data for wildlife monitoring 
  #'  --------------------------------------------------------------------------
  #'  Script to create unmarked data frames and run single-species, single-season
  #'  occupancy models for deer, elk, cougars, wolves, coyotes, and bobcats for
  #'  summer (July 1 - Sept 29, 2018 & 2019) and winter (Dec 1 - Mar 1, 2018 - 2019 
  #'  & 2019 - 2020), respectively. Each single-season occupancy model includes 
  #'  13 7-day sampling occasions comprising the warmest months in summer and
  #'  coldest months in winter with the most consistent snow. Result tables 
  #'  created at end of script.
  #'  
  #'  MUST download Data folder from Dryad repository for script to work.
  #'  
  #'  Camera data collected as part of the Washington Predator-Prey Project:
  #'  https://predatorpreyproject.weebly.com/
  #'  Covariate data included in occupancy models were collected at each
  #'  camera site or extracted from remotely sensed data.
  #'  --------------------------------------------------------------------------

  #'  Load libraries
  library(unmarked)
  library(tidyverse)
  library(data.table)

  #'  Read in detection histories
  DH_files <- list.files("./Data/Detection Histories/", pattern=".csv", full.names = TRUE)
  DH_files # this will be the order of the DH_list!
  DH_list <- lapply(DH_files, read.csv)
  
  #'  Read in covariate data & scale
  #'  Elevation, Slope, Road Density, % Forest, Grass, Shrub => site-level occupancy covs
  #'  Distance, Height, monitoring => survey/site-level detection covs
  stations <- read.csv("./Data/Covariate Data/Camera_Station_Covariates.csv") %>% 
    mutate(
      CameraLocation = as.factor(CameraLocation),
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      Distance = scale(Distance),
      Height = scale(Height),
      Trail = as.factor(Trail),
      PercForestMix = scale(PercForestMix),
      PercXericShrub = scale(PercXericShrub),
      PercXericGrass = scale(PercXericGrass),
      Elev = scale(Elev), 
      Slope = scale(Slope), 
      RoadDensity = scale(RoadDensity),
    ) %>%
    arrange(Year) # necessary to match camera location order in detection histories

  #'  Adjust reference category for Monitoring factors
  order_trail <- c("Trail", "Dirt road", "Decommissioned road")
  stations <- stations %>%
    mutate(
      Trail = fct_relevel(Trail, order_trail)
    )
    
  #'  Missing camera height & distance values for two camera sites 
  #'  Replacing NAs with mean camera height and distance (0 when scaled)
  stations$Distance[is.na(stations$Distance),] <- 0
  stations$Height[is.na(stations$Height),] <- 0
  
  #'  Save study-area specific covariates for deer & elk models
  stations_NE <- filter(stations, Study_Area == "NE")
  stations_OK <- filter(stations, Study_Area == "OK")
  
  #'  Check for correlation among continuous covariates
  #'  Using r = |0.7| as cutoff 
  cov_data <- stations[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity")]
  (cov_matrix <- cor(cov_data, use = "complete.obs")) 
  cov_data_NE <- stations_NE[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity")]
  (cov_matrix_NE <- cor(cov_data_NE, use = "complete.obs"))
  cov_data_OK <- stations_OK[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity")]
  (cov_matrix_OK <- cor(cov_data_OK, use = "complete.obs")) 
  
  #'  Survey-level mean weekly temperature per camera site
  temp_smr <- read.csv("./Data/Covariate Data/Sampling_Occasion_Summer_Temperatures.csv") %>%
    dplyr::select(-X)
  Temp_smr <- temp_smr %>%
    dplyr::select(-c(CameraLocation, Study_Area)) 
  temp_wtr <- read.csv("./Data/Covariate Data/Sampling_Occasion_Winter_Temperatures.csv") %>%
    dplyr::select(-X)
  Temp_wtr <- temp_wtr %>%
    dplyr::select(-c(CameraLocation, Study_Area)) 

  #'  Function to scale survey-level covariates
  scale_srvy_cov <- function(srvycov) {
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(srvycov), na.rm = TRUE)
    sd <- sd(as.matrix(srvycov), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((srvycov - mu) / sd)
    
    return(scaled)
  }
  #'  Scale each survey-level covariate
  Temp_smr_scaled <- scale_srvy_cov(Temp_smr)
  Temp_wtr_scaled <- scale_srvy_cov(Temp_wtr)
  
  #'  Create survey-level covariate matrix
  #'  Requires unique column for each sampling occasion and covariate
  nrows <- nrow(stations)
  ncols <- 13
  
  srvy_covs <- list(
    Height = matrix(c(Hgt1 = stations$Height, Hgt2 = stations$Height,
                      Hgt3 = stations$Height, Hgt4 = stations$Height,
                      Hgt5 = stations$Height, Hgt6 = stations$Height,
                      Hgt7 = stations$Height, Hgt8 = stations$Height,
                      Hgt9 = stations$Height, Hgt10 = stations$Height,
                      Hgt11 = stations$Height, Hgt12 = stations$Height,
                      Hgt13 = stations$Height),
                    nrow = nrows, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = stations$Distance, Dist2 = stations$Distance,
                        Dist3 = stations$Distance, Dist4 = stations$Distance,
                        Dist5 = stations$Distance, Dist6 = stations$Distance,
                        Dist7 = stations$Distance, Dist8 = stations$Distance,
                        Dist9 = stations$Distance, Dist10 = stations$Distance,
                        Dist11 = stations$Distance, Dist12 = stations$Distance,
                        Dist13 = stations$Distance),
                      nrow = nrows, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr_scaled,
    Temp_wtr = Temp_wtr_scaled
    )

  #'  Create study-area specific survey covariates & call scaling function
  Hgt_NE <- stations$Height[stations$Study_Area == "NE"]
  Hgt_OK <- stations$Height[stations$Study_Area == "OK"]
  Dist_NE <- stations$Distance[stations$Study_Area == "NE"]
  Dist_OK <- stations$Distance[stations$Study_Area == "OK"]
  
  Temp_smr_NE <- temp_smr[temp_smr$Study_Area == "NE",] %>%
    dplyr::select(-c(CameraLocation, Study_Area)) 
  Temp_smr_OK <- temp_smr[temp_smr$Study_Area == "OK",] %>%
    dplyr::select(-c(CameraLocation, Study_Area)) 
  Temp_wtr_NE <- temp_wtr[temp_wtr$Study_Area == "NE",] %>%
    dplyr::select(-c(CameraLocation, Study_Area)) 
  Temp_wtr_OK <- temp_wtr[temp_wtr$Study_Area == "OK",] %>%
    dplyr::select(-c(CameraLocation, Study_Area))
  
  Temp_smr_NE_scaled <- scale_srvy_cov(Temp_smr_NE)
  Temp_smr_OK_scaled <- scale_srvy_cov(Temp_smr_OK)
  Temp_wtr_NE_scaled <- scale_srvy_cov(Temp_wtr_NE)
  Temp_wtr_OK_scaled <- scale_srvy_cov(Temp_wtr_OK)
  
  nrows_NE <- nrow(stations_NE)
  nrows_OK <- nrow(stations_OK)
  
  #'  Format survey-level covariate data for study area-specific occupancy models
  #'  (important for deer & elk models)
  srvy_covs_NE <- list(
    Height = matrix(c(Hgt1 = Hgt_NE, Hgt2 = Hgt_NE, Hgt3 = Hgt_NE, Hgt4 = Hgt_NE,
                      Hgt5 = Hgt_NE, Hgt6 = Hgt_NE, Hgt7 = Hgt_NE, Hgt8 = Hgt_NE,
                      Hgt9 = Hgt_NE, Hgt10 = Hgt_NE, Hgt11 = Hgt_NE, Hgt12 = Hgt_NE,
                      Hgt13 = Hgt_NE),
                    nrow = nrows_NE, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = Dist_NE, Dist2 = Dist_NE, Dist3 = Dist_NE, 
                        Dist4 = Dist_NE, Dist5 = Dist_NE, Dist6 = Dist_NE,
                        Dist7 = Dist_NE, Dist8 = Dist_NE, Dist9 = Dist_NE, 
                        Dist10 = Dist_NE, Dist11 = Dist_NE, Dist12 = Dist_NE,
                        Dist13 = Dist_NE),
                      nrow = nrows_NE, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr_NE_scaled,
    Temp_wtr = Temp_wtr_NE_scaled
  )
  
  srvy_covs_OK <- list(
    Height = matrix(c(Hgt1 = Hgt_OK, Hgt2 = Hgt_OK, Hgt3 = Hgt_OK, Hgt4 = Hgt_OK,
                      Hgt5 = Hgt_OK, Hgt6 = Hgt_OK, Hgt7 = Hgt_OK, Hgt8 = Hgt_OK,
                      Hgt9 = Hgt_OK, Hgt10 = Hgt_OK, Hgt11 = Hgt_OK, Hgt12 = Hgt_OK,
                      Hgt13 = Hgt_OK),
                    nrow = nrows_OK, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = Dist_OK, Dist2 = Dist_OK, Dist3 = Dist_OK, 
                        Dist4 = Dist_OK, Dist5 = Dist_OK, Dist6 = Dist_OK,
                        Dist7 = Dist_OK, Dist8 = Dist_OK, Dist9 = Dist_OK, 
                        Dist10 = Dist_OK, Dist11 = Dist_OK, Dist12 = Dist_OK,
                        Dist13 = Dist_OK),
                      nrow = nrows_OK, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr_OK_scaled,
    Temp_wtr = Temp_wtr_OK_scaled
  )

  
  #'  Create unmarked data frames
  #'  --------------------------------------------------------------------------
  #'  Format detection histories and covariate data for unmarked
  #'  --------------------------------------------------------------------------
  ####  BOBCAT UMFs  ####
  bob_s1819_UMF <- unmarkedFrameOccu(DH_list[[1]][,-1], #' [,-1] drops CameraLocation column
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         PercForMix = stations$PercForestMix,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         RoadDensity = stations$RoadDensity),
                                   obsCovs = srvy_covs)

  bob_w1820_UMF <- unmarkedFrameOccu(DH_list[[2]][,-1],
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           PercForMix = stations$PercForestMix,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           RoadDensity = stations$RoadDensity),
                                     obsCovs = srvy_covs)
  summary(bob_w1820_UMF)
  
  ####  COUGAR UMFs  ####
  coug_s1819_UMF <- unmarkedFrameOccu(DH_list[[3]][,-1],
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity),
                                    obsCovs = srvy_covs)
  
  coug_w1820_UMF <- unmarkedFrameOccu(DH_list[[4]][,-1],
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity),
                                    obsCovs = srvy_covs)
  
  
  ####  COYOTE UMFs  ####
  coy_s1819_UMF <- unmarkedFrameOccu(DH_list[[5]][,-1],
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity),
                                    obsCovs = srvy_covs)
  
  coy_w1820_UMF <- unmarkedFrameOccu(DH_list[[6]][,-1],
                                      siteCovs = data.frame(Year = stations$Year,
                                                            Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            PercForMix = stations$PercForestMix,
                                                            PercXGrass = stations$PercXericGrass,
                                                            PercXShrub = stations$PercXericShrub,
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            RoadDensity = stations$RoadDensity),
                                      obsCovs = srvy_covs)
  
  ####  ELK UMFs  ####
  #'  Data from NE study area only to mirror GPS collar data
  elk_s1819_UMF <- unmarkedFrameOccu(DH_list[[7]][,-1],
                                     siteCovs = data.frame(Year = stations_NE$Year,
                                                           Trail = stations_NE$Trail,
                                                           PercForMix = stations_NE$PercForestMix,
                                                           PercXGrass = stations_NE$PercXericGrass,
                                                           PercXShrub = stations_NE$PercXericShrub,
                                                           Elev = stations_NE$Elev,
                                                           Slope = stations_NE$Slope,
                                                           RoadDensity = stations_NE$RoadDensity),
                                     obsCovs = srvy_covs_NE)
  
  elk_w1820_UMF <- unmarkedFrameOccu(DH_list[[8]][,-1],
                                     siteCovs = data.frame(Year = stations_NE$Year,
                                                           Trail = stations_NE$Trail,
                                                           PercForMix = stations_NE$PercForestMix,
                                                           PercXGrass = stations_NE$PercXericGrass,
                                                           PercXShrub = stations_NE$PercXericShrub,
                                                           Elev = stations_NE$Elev,
                                                           Slope = stations_NE$Slope,
                                                           RoadDensity = stations_NE$RoadDensity),
                                     obsCovs = srvy_covs_NE)
  
  ####  MULE DEER UMFs  ####
  #'  Data from OK study area only to mirror GPS collar data
  md_s1819_UMF <- unmarkedFrameOccu(DH_list[[9]][,-1],
                                    siteCovs = data.frame(Year = stations_OK$Year,
                                                          Trail = stations_OK$Trail,
                                                          PercForMix = stations_OK$PercForestMix,
                                                          PercXGrass = stations_OK$PercXericGrass,
                                                          PercXShrub = stations_OK$PercXericShrub,
                                                          Elev = stations_OK$Elev,
                                                          Slope = stations_OK$Slope,
                                                          RoadDensity = stations_OK$RoadDensity),
                                    obsCovs = srvy_covs_OK)
  
  md_w1820_UMF <- unmarkedFrameOccu(DH_list[[10]][,-1],
                                    siteCovs = data.frame(Year = stations_OK$Year,
                                                          Trail = stations_OK$Trail,
                                                          PercForMix = stations_OK$PercForestMix,
                                                          PercXGrass = stations_OK$PercXericGrass,
                                                          PercXShrub = stations_OK$PercXericShrub,
                                                          Elev = stations_OK$Elev,
                                                          Slope = stations_OK$Slope,
                                                          RoadDensity = stations_OK$RoadDensity),
                                    obsCovs = srvy_covs_OK)
  
  ####  WHITE-TAILED DEER UMFs  ####
  #'  Data from NE study area only to mirror GPS collar data
  wtd_s1819_UMF <- unmarkedFrameOccu(DH_list[[11]][,-1],
                                     siteCovs = data.frame(Year = stations_NE$Year,
                                                           Trail = stations_NE$Trail,
                                                           PercForMix = stations_NE$PercForestMix,
                                                           PercXGrass = stations_NE$PercXericGrass,
                                                           PercXShrub = stations_NE$PercXericShrub,
                                                           Elev = stations_NE$Elev,
                                                           Slope = stations_NE$Slope,
                                                           RoadDensity = stations_NE$RoadDensity),
                                     obsCovs = srvy_covs_NE)
  
  wtd_w1820_UMF <- unmarkedFrameOccu(DH_list[[12]][,-1],
                                     siteCovs = data.frame(Year = stations_NE$Year,
                                                           Trail = stations_NE$Trail,
                                                           PercForMix = stations_NE$PercForestMix,
                                                           PercXGrass = stations_NE$PercXericGrass,
                                                           PercXShrub = stations_NE$PercXericShrub,
                                                           Elev = stations_NE$Elev,
                                                           Slope = stations_NE$Slope,
                                                           RoadDensity = stations_NE$RoadDensity),
                                     obsCovs = srvy_covs_NE)
  
  
  ####  WOLF UMFs  ####
  wolf_s1819_UMF <- unmarkedFrameOccu(DH_list[[13]][,-1],
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         PercForMix = stations$PercForestMix,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         RoadDensity = stations$RoadDensity),
                                   obsCovs = srvy_covs)
  
  wolf_w1820_UMF <- unmarkedFrameOccu(DH_list[[14]][,-1],
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           PercForMix = stations$PercForestMix,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           RoadDensity = stations$RoadDensity),
                                     obsCovs = srvy_covs)
  
  
  ####  Occupancy models  ####
  #'  --------------------------------------------------------------------------
  #'  unmarked formula: ~detection ~occupancy
  #'  
  #'  Run occupancy model for all species and seasons. Removed Study Area covariate 
  #'  if data only collected in one area. Additional habitat covariates removed 
  #'  only if failed to converge or were highly correlated with another covariate 
  #'  (usually with elevation).
  #'  --------------------------------------------------------------------------

  ####  BOBCAT MODELS  ####                   
  #'  SUMMERS 2018 & 2019
  (bob_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height
                         ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area,
                         bob_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(bob_s1819_occ, type = "state")
  unmarked::vif(bob_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020
  (bob_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height
                         ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + Area,
                         bob_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(bob_w1820_occ, type = "state")
  unmarked::vif(bob_w1820_occ, type = "det")
  
  
  ####  COUGAR MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coug_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Height*Distance 
                          ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area, 
                          coug_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coug_s1819_occ, type = "state")
  unmarked::vif(coug_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020     
  (coug_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height
                          ~Elev + Slope + PercForMix + PercXShrub + RoadDensity + Area, 
                          coug_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coug_w1820_occ, type = "state") 
  unmarked::vif(coug_w1820_occ, type = "det")
  
  ####  COYOTE MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coy_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height
                         ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area, 
                         coy_s1819_UMF)) 
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coy_s1819_occ, type = "state") 
  unmarked::vif(coy_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020    
  (coy_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                         ~Elev + Slope + PercForMix + PercXShrub + RoadDensity + Area,
                         coy_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coy_w1820_occ, type = "state") 
  unmarked::vif(coy_w1820_occ, type = "det")

  
  ####  ELK MODELS ####                          
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  (elk_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                         ~Elev + Slope + PercForMix + RoadDensity, 
                         elk_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(elk_s1819_occ, type = "state") 
  unmarked::vif(elk_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  NE study area only so no Area effect 
  (elk_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                         ~Elev + Slope + PercForMix + RoadDensity,
                         elk_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(elk_w1820_occ, type = "state") 
  unmarked::vif(elk_w1820_occ, type = "det")
  
  
  ####  MULE DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  OK study area only so no Area effect
  (md_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height
                        ~Elev + Slope + PercForMix + PercXGrass + RoadDensity,
                        md_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(md_s1819_occ, type = "state") 
  unmarked::vif(md_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  OK study area only so no Area effect                    
  (md_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                        ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity, 
                        md_w1820_UMF))  
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(md_w1820_occ, type = "state") 
  unmarked::vif(md_w1820_occ, type = "det")
  
  
  ####  WHITE-TAILED DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  (wtd_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                         ~Elev + Slope + PercForMix +  RoadDensity, 
                         wtd_s1819_UMF)) 
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wtd_s1819_occ, type = "state") 
  unmarked::vif(wtd_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  NE study area only so no Area effect     
  (wtd_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height
                         ~Elev + Slope + PercForMix + RoadDensity,
                         wtd_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wtd_w1820_occ, type = "state") 
  unmarked::vif(wtd_w1820_occ, type = "det")
  
  
  ####  WOLF MODELS  ####
  #'  SUMMERS 2018 & 2019    
  (wolf_s1819_occ <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height
                          ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + Area, 
                          wolf_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wolf_s1819_occ, type = "state") 
  unmarked::vif(wolf_s1819_occ, type = "det")
  
  #'  WINTERS 2018-2019 & 2019-2020 
  (wolf_w1820_occ <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                          ~Elev + Slope + PercForMix + RoadDensity + Area, 
                          wolf_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wolf_w1820_occ, type = "state") 
  unmarked::vif(wolf_w1820_occ, type = "det")

  
  ####  Summary tables  ####
  #'  --------------------------------------------------------------------------
  #'  Function to extract coefficient estimates & p-values from occupancy sub-model
  occ_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    return(out)
  }
  
  #'  Run each model through function
  bob_s1819_occ_tbl <- occ_out(bob_s1819_occ, "Bobcat", "Summer") 
  bob_w1820_occ_tbl <- occ_out(bob_w1820_occ, "Bobcat", "Winter")
  coug_s1819_occ_tbl <- occ_out(coug_s1819_occ, "Cougar", "Summer")
  coug_w1820_occ_tbl <- occ_out(coug_w1820_occ, "Cougar", "Winter")
  coy_s1819_occ_tbl <- occ_out(coy_s1819_occ, "Coyote", "Summer")
  coy_w1820_occ_tbl <- occ_out(coy_w1820_occ, "Coyote", "Winter")
  wolf_s1819_occ_tbl <- occ_out(wolf_s1819_occ, "Wolf", "Summer")
  wolf_w1820_occ_tbl <- occ_out(wolf_w1820_occ, "Wolf", "Winter")
  elk_s1819_occ_tbl <- occ_out(elk_s1819_occ, "Elk", "Summer")
  elk_w1820_occ_tbl <- occ_out(elk_w1820_occ, "Elk", "Winter")
  md_s1819_occ_tbl <- occ_out(md_s1819_occ, "Mule Deer", "Summer")
  md_w1820_occ_tbl <- occ_out(md_w1820_occ, "Mule Deer", "Winter")
  wtd_s1819_occ_tbl <- occ_out(wtd_s1819_occ, "White-tailed Deer", "Summer")
  wtd_w1820_occ_tbl <- occ_out(wtd_w1820_occ, "White-tailed Deer", "Winter")
  
  #'  Merge into larger data frames for easy comparison
  summer_occ_tbl <- rbind(bob_s1819_occ_tbl, coug_s1819_occ_tbl, coy_s1819_occ_tbl, 
                          elk_s1819_occ_tbl, md_s1819_occ_tbl, wtd_s1819_occ_tbl,
                          wolf_s1819_occ_tbl)
  winter_occ_tbl <- rbind(bob_w1820_occ_tbl, coug_w1820_occ_tbl, coy_w1820_occ_tbl, 
                          elk_w1820_occ_tbl, md_w1820_occ_tbl, wtd_w1820_occ_tbl,
                          wolf_w1820_occ_tbl)
  #'  One big table
  occ_results <- rbind(summer_occ_tbl, winter_occ_tbl) %>%
    arrange(Species)
  colnames(occ_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")  
  
  #'  Round so estimates are easier to read
  rounddig <- 2
  results_psi <- occ_results %>%
    mutate(
      Estimate = round(Estimate, rounddig),
      SE = round(SE, rounddig),
      z = round(z, rounddig),
      Pval = round(Pval, rounddig)
    )
  
  #'  Save for predicting probability of occupancy across study area
  write.csv(results_psi, "./Data/Results tables/OccMod_OccProb_Results.csv")  
  
  
  ####  Covariate effects on detection probability  ####
  #'  --------------------------------------------------------------------------
  #'  Function to extract coefficient estimates & p-values from detection sub-model
  det_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    return(out)
  }
  
  #'  Run each model through detection function
  bob_s1819_det <- det_out(bob_s1819_occ, "Bobcat", "Summer") 
  bob_w1820_det <- det_out(bob_w1820_occ, "Bobcat", "Winter")
  coug_s1819_det <- det_out(coug_s1819_occ, "Cougar", "Summer")
  coug_w1820_det <- det_out(coug_w1820_occ, "Cougar", "Winter")
  coy_s1819_det <- det_out(coy_s1819_occ, "Coyote", "Summer")
  coy_w1820_det <- det_out(coy_w1820_occ, "Coyote", "Winter")
  wolf_s1819_det <- det_out(wolf_s1819_occ, "Wolf", "Summer")
  wolf_w1820_det <- det_out(wolf_w1820_occ, "Wolf", "Winter")
  elk_s1819_det <- det_out(elk_s1819_occ, "Elk", "Summer")
  elk_w1820_det <- det_out(elk_w1820_occ, "Elk", "Winter")
  md_s1819_det <- det_out(md_s1819_occ, "Mule Deer", "Summer")
  md_w1820_det <- det_out(md_w1820_occ, "Mule Deer", "Winter")
  wtd_s1819_det <- det_out(wtd_s1819_occ, "White-tailed Deer", "Summer")
  wtd_w1820_det <- det_out(wtd_w1820_occ, "White-tailed Deer", "Winter")

  #'  Merge into larger data frames for easier comparison
  summer_det <- rbind(bob_s1819_det, coug_s1819_det, coy_s1819_det, wolf_s1819_det,
                      elk_s1819_det, md_s1819_det, wtd_s1819_det)
  winter_det <- rbind(bob_w1820_det, coug_w1820_det, coy_w1820_det, wolf_w1820_det,
                      elk_w1820_det, md_w1820_det, wtd_w1820_det)
  #' One big table
  det_results <- rbind(summer_det, winter_det) %>%
    arrange(Species)
  colnames(det_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  
  #'  Round so values are a little easier to interpret
  results_det <- det_results %>%
    mutate(
      Estimate = round(Estimate, 2),
      SE = round(SE, 2),
      z = round(z, 2),
      Pval = round(Pval, 2)
    )
 
  
  ####  Mean occupancy and detection probabilities  ####
  #'  --------------------------------------------------------------------------
  mu_occ <- function(mod, species, season) {
    #'  Predict occupancy probability for all camera sties
    occu_mean <- predict(object = mod, type = "state") %>% 
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE)
    #'  Predict occupancy probability for all camera sties
    det_mean <- predict(object = mod, type = "det") %>%
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE) 
    predicted <- as.data.frame(rbind(occu_mean, det_mean))
    colnames(predicted) <- c("Mean", "SE")
    Parameter <- c("Occupancy", "Detection")
    Species <- species
    Season <- season
    predicted <- cbind(predicted, Parameter)
    predicted <- cbind(predicted, Species)
    predicted <- cbind(predicted, Season)
    return(predicted)
  }
  #'  Estimate mean probability of occupancy and detection per species and season
  md_predict_smr <- mu_occ(md_s1819_occ, "Mule Deer", "Summer")
  md_predict_wtr <- mu_occ(md_w1820_occ, "Mule Deer", "Winter")
  elk_predict_smr <- mu_occ(elk_s1819_occ, "Elk", "Summer")
  elk_predict_wtr <- mu_occ(elk_w1820_occ, "Elk", "Winter")
  wtd_predict_smr <- mu_occ(wtd_s1819_occ, "White-tailed Deer", "Summer")
  wtd_predict_wtr <- mu_occ(wtd_w1820_occ, "White-tailed Deer", "Winter")
  coug_predict_smr <- mu_occ(coug_s1819_occ, "Cougar", "Summer")
  coug_predict_wtr <- mu_occ(coug_w1820_occ, "Cougar", "Winter")
  wolf_predict_smr <- mu_occ(wolf_s1819_occ, "Wolf", "Summer")
  wolf_predict_wtr <- mu_occ(wolf_w1820_occ, "Wolf", "Winter")
  bob_predict_smr <- mu_occ(bob_s1819_occ, "Bobcat", "Summer")
  bob_predict_wtr <- mu_occ(bob_w1820_occ, "Bobcat", "Winter")
  coy_predict_smr <- mu_occ(coy_s1819_occ, "Coyote", "Summer")
  coy_predict_wtr <- mu_occ(coy_w1820_occ, "Coyote", "Winter")

  #'  Make a pretty table
  Mean_tbl <- bind_rows(md_predict_smr, md_predict_wtr, elk_predict_smr, elk_predict_wtr, wtd_predict_smr, 
              wtd_predict_wtr, coug_predict_smr, coug_predict_wtr, wolf_predict_smr, 
              wolf_predict_wtr, bob_predict_smr, bob_predict_wtr, coy_predict_smr, 
              coy_predict_wtr) %>%
    relocate(Species, .before = Mean) %>%
    relocate(Season, .after = Species) %>%
    relocate(Parameter, .after = Season) %>%
    arrange(Parameter, Mean, Species)
  
  #'  End
  