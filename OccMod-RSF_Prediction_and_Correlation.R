  #'  R code associated with: 
  #'  
  #'  Bassing, S. B., M. Devivo, T. R. Ganz, B. N. Kertson, L. R. Prugh, T. Roussin, 
  #'  L. Satterfield, R. M. Windell, A. J. Wirsing, & B. Gardner 
  #'  
  #'  Are we telling the same story? Comparing inferences made from camera trap 
  #'  and telemetry data for wildlife monitoring 
  #'  --------------------------------------------------------------------------
  #'  Script to predict probability of use and relative probability of selection
  #'  for each species across the relevant study area(s) by season, based on 
  #'  occupancy model and RSF results. Script scales covariate data across both
  #'  study areas based on mean & SD values from each species & season analysis. 
  #'  Script then predicts probability of use and relative selection across
  #'  study areas and calculates pixel by pixel correlation between results. 
  #'  
  #'  MUST download data from Dryad repository and create "Data" folder for script to work.
  #'  --------------------------------------------------------------------------

  #'  Load libraries
  library(tidyverse)
  
  ####  Read in data  ####
  #'  Occupancy model output
  occ_out <- read.csv("./Data/OccMod_OccProb_Results.csv") %>%
    #'  Calculate 90% confidence intervals
    mutate(
      l90 = (Estimate - (1.64 * SE)),  
      u90 = (Estimate + (1.64 * SE))   
    ) 
  
  #'  RSF results output
  rsf_out <- read.csv("./Data/RSF_Results.csv") %>% 
    #'  Calculate 95% confidence intervals
    mutate(
      l95 = (Estimate - (1.96 * SE)),  
      u95 = (Estimate + (1.96 * SE))
    ) 
  
  #'  Original covariate data from occupancy models
  stations <- read.csv("./Data/Camera_Station_Covariates.csv") %>%
    transmute(
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      CameraLocation = as.factor(CameraLocation),
      PercForMix = PercForestMix, 
      PercXShrub = PercXericShrub,
      PercXGrass = PercXericGrass,
      Elev = as.numeric(Elev), 
      Slope = Slope, 
      RoadDen = RoadDensity
    ) %>%
    arrange(Year)
  
  #'  Original covariate data from RSFs
  load("./Data/RSF_used_avail_pts_noXY.RData")
  
  #'  Covariate values for each pixel in study area for predicting across study area
  #'  Using 1km grid cells for ease of computation but values for 30m grid cells
  #'  used in full analysis are provided as well
  NE_covs_1km <- read.csv("./Data/StudyAreaWide_NE_Covariates_1km.csv") %>% 
    #'  Add binary study area indicator for occupancy predictions
    mutate(Area = 0)
  OK_covs_1km <- read.csv("./Data/StudyAreaWide_OK_Covariates_1km.csv") %>%
    #'  Add binary study area indicator for occupancy predictions
    mutate(Area = 1)
  all_covs_1km <- as.data.frame(rbind(NE_covs_1km, OK_covs_1km))
  
  #'  Exclude elevation values in cells >2150 m elevation for camera covariates 
  #'  because not sampled by camera traps
  all_covs_adj_1km <- all_covs_1km %>%
    mutate(Elev = ifelse(Elev > 2150, NA, Elev))

  
  ####  Center & Scale Covariate Data  ####
  #'  --------------------------------------------------------------------------
  #'  Function to find mean and standard deviation for each covariate
  #'  Used when center & scaling covariates for original models and will be used 
  #'  to standardize covariate data when predicting across study areas
  cov_summary <- function(covs) {
    mu.cov <- covs %>%
      summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
    sd.cov <- covs %>%
      summarise(across(where(is.numeric), ~sd(.x, na.rm = TRUE)))
    mu.sd.cov <- rbind(mu.cov, sd.cov)
    parameter <- as.data.frame(c("Mean", "SD"))
    colnames(parameter) <- "Parameter"
    cov_summary <- cbind(parameter, mu.sd.cov)
    return(cov_summary)
  }
  #'  Find mean & SD across all camera locations
  summary_occ_covs <- cov_summary(stations)
  #'  Find mean & SD across all used & available collar data, but separately by season & species
  summary_bob_smr <- cov_summary(rsf_pts[[1]][rsf_pts[[1]]$Season == "Summer18" | rsf_pts[[1]]$Season == "Summer19",])
  summary_bob_wtr <- cov_summary(rsf_pts[[1]][rsf_pts[[1]]$Season == "Winter1819" | rsf_pts[[1]]$Season == "Winter1920",])
  summary_coug_smr <- cov_summary(rsf_pts[[2]][rsf_pts[[2]]$Season == "Summer18" | rsf_pts[[2]]$Season == "Summer19",])
  summary_coug_wtr <- cov_summary(rsf_pts[[2]][rsf_pts[[2]]$Season == "Winter1819" | rsf_pts[[2]]$Season == "Winter1920",])
  summary_coy_smr <- cov_summary(rsf_pts[[3]][rsf_pts[[3]]$Season == "Summer18" | rsf_pts[[3]]$Season == "Summer19",])
  summary_coy_wtr <- cov_summary(rsf_pts[[3]][rsf_pts[[3]]$Season == "Winter1819" | rsf_pts[[3]]$Season == "Winter1920",])
  summary_elk_smr <- cov_summary(rsf_pts[[4]][rsf_pts[[4]]$Season == "Summer18" | rsf_pts[[4]]$Season == "Summer19",])
  summary_elk_wtr <- cov_summary(rsf_pts[[4]][rsf_pts[[4]]$Season == "Winter1819" | rsf_pts[[4]]$Season == "Winter1920",])
  summary_md_smr <- cov_summary(rsf_pts[[5]][rsf_pts[[5]]$Season == "Summer18" | rsf_pts[[5]]$Season == "Summer19",])
  summary_md_wtr <- cov_summary(rsf_pts[[5]][rsf_pts[[5]]$Season == "Winter1819" | rsf_pts[[5]]$Season == "Winter1920",])
  summary_wtd_smr <- cov_summary(rsf_pts[[6]][rsf_pts[[6]]$Season == "Summer18" | rsf_pts[[6]]$Season == "Summer19",])
  summary_wtd_wtr <- cov_summary(rsf_pts[[6]][rsf_pts[[6]]$Season == "Winter1819" | rsf_pts[[6]]$Season == "Winter1920",])
  summary_wolf_smr <- cov_summary(rsf_pts[[7]][rsf_pts[[7]]$Season == "Summer18" | rsf_pts[[7]]$Season == "Summer19",])
  summary_wolf_wtr <- cov_summary(rsf_pts[[7]][rsf_pts[[7]]$Season == "Winter1819" | rsf_pts[[7]]$Season == "Winter1920",])
    
  #'  Function to scale the individual covariates based on the covariate-specific 
  #'  mean and SD values used to standardize the covariate data in the original
  #'  models. This will differ by species and season for the RSF predictions since
  #'  the used and available locations differed with each model. This will be the
  #'  same for all occupancy model predictions since the same camera sites were
  #'  included in all analyses.
  scaling_covs <- function(covs, mu.sd) {
    scaling_covs <- covs %>%
      transmute(
        obs = obs,
        Elev = (Elev - mu.sd$Elev[1]) / mu.sd$Elev[2],
        Slope = (Slope - mu.sd$Slope[1]) / mu.sd$Slope[2],
        RoadDen = (RoadDen - mu.sd$RoadDen[1]) / mu.sd$RoadDen[2],
        PercForMix = (PercForestMix - mu.sd$PercForMix[1]) / mu.sd$PercForMix[2],
        PercXGrass = (PercXericGrass - mu.sd$PercXGrass[1]) / mu.sd$PercXGrass[2],
        PercXShrub = (PercXericShrub - mu.sd$PercXShrub[1]) / mu.sd$PercXShrub[2],
        x = x,
        y = y,
        Area = Area)
  }
  #'  Standardize covariate data based on cameras
  cam_zcovs <- scaling_covs(all_covs_adj_1km, summary_occ_covs)  
  #'  Standardize covariate data based on used/available locations
  bob_smr_zcovs <- scaling_covs(all_covs_1km, summary_bob_smr)
  bob_wtr_zcovs <- scaling_covs(all_covs_1km, summary_bob_wtr)
  coug_smr_zcovs <- scaling_covs(all_covs_1km, summary_coug_smr)
  coug_wtr_zcovs <- scaling_covs(all_covs_1km, summary_coug_wtr)
  coy_smr_zcovs <- scaling_covs(all_covs_1km, summary_coy_smr)
  coy_wtr_zcovs <- scaling_covs(all_covs_1km, summary_coy_wtr)
  elk_smr_zcovs <- scaling_covs(all_covs_1km, summary_elk_smr)
  elk_wtr_zcovs <- scaling_covs(all_covs_1km, summary_elk_wtr)
  md_smr_zcovs <- scaling_covs(all_covs_1km, summary_md_smr)  
  md_wtr_zcovs <- scaling_covs(all_covs_1km, summary_md_wtr)
  wtd_smr_zcovs <- scaling_covs(all_covs_1km, summary_wtd_smr)
  wtd_wtr_zcovs <- scaling_covs(all_covs_1km, summary_wtd_wtr)
  wolf_smr_zcovs <- scaling_covs(all_covs_1km, summary_wolf_smr)
  wolf_wtr_zcovs <- scaling_covs(all_covs_1km, summary_wolf_wtr)
  
  zcovs <- list(cam_zcovs, bob_smr_zcovs, bob_wtr_zcovs, coug_smr_zcovs, coug_wtr_zcovs,
                coy_smr_zcovs, coy_wtr_zcovs, elk_smr_zcovs, elk_wtr_zcovs, 
                md_smr_zcovs, md_wtr_zcovs, wtd_smr_zcovs, wtd_wtr_zcovs, wolf_smr_zcovs,
                wolf_wtr_zcovs)

  
  ####  Predict probability of use across study areas  ####
  #'  --------------------------------------------------------------------------
  #'  Manipulate occupancy result table
  occ_coefs_signif <- occ_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDensity,
      B.area = AreaOK) %>% 
    #'  Change NAs to 0 (no effect) for coefficients not included in species-specific models
    mutate(
      B.elev = ifelse(is.na(B.elev), 0, B.elev),
      B.slope = ifelse(is.na(B.slope), 0, B.slope),
      B.for = ifelse(is.na(B.for), 0, B.for),
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.rd = ifelse(is.na(B.rd), 0, B.rd),
      B.area = ifelse(is.na(B.area), 0, B.area)
    )
  
  #'  Function to predict across all grid cells based on occupancy model results
  #'  Should end up with 1 predicted value per grid cell
  predict_occ <- function(cov, coef) {
    predict_odds <- c()
    predict_prob <- c()
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$alpha + coef$B.area*cov$Area[i] + coef$B.elev*cov$Elev[i] + 
                               coef$B.slope*cov$Slope[i]+ coef$B.for*cov$PercForMix[i] + 
                               coef$B.grass*cov$PercXGrass[i] + coef$B.shrub*cov$PercXShrub[i] + 
                               coef$B.rd*cov$RoadDen[i])
      predict_prob[i] <- predict_odds[i] / (1 + predict_odds[i])
    }
    predict_prob <- as.data.frame(predict_prob) %>%
      transmute(
        Predicted_Occ = predict_prob
      )
    return(predict_prob)
  }
  #'  Run estimates from occupancy sub-model through function to predict probability of use
  #'  Reminder: zcovs[[1]] is standardized based on camera sites and applies to all occ mods
  #'  Area == 1 if only predicting across Okanogan study area, Area == 0 if only predicting
  #'  across Northeast study area
  md_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Summer",])
  md_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Winter",])
  elk_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Summer",])
  elk_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Winter",])
  wtd_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Winter",])
  coug_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Summer",])
  coug_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Winter",])
  wolf_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Winter",])
  bob_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Summer",])
  bob_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Winter",])
  coy_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Summer",])
  coy_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Winter",])
  
  #'  Combine into a large data frame
  #'  Start with predator data that spans both study areas
  Predicted_occ <- as.data.frame(zcovs[[1]]) %>%  
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(bob_smr_predict_occ_sgnf, bob_wtr_predict_occ_sgnf,
          coug_smr_predict_occ_sgnf, coug_wtr_predict_occ_sgnf, 
          coy_smr_predict_occ_sgnf, coy_wtr_predict_occ_sgnf,
          wolf_smr_predict_occ_sgnf, wolf_wtr_predict_occ_sgnf)
  #'  Make sure you have the order right when you change the names!
  colnames(Predicted_occ) <- c("obs", "Area", "x", "y",  
                               "BOB_smr_occ", "BOB_wtr_occ", 
                               "COUG_smr_occ", "COUG_wtr_occ", 
                               "COY_smr_occ", "COY_wtr_occ",
                               "WOLF_smr_occ", "WOLF_wtr_occ")
  
  #'  Okanogan-only data (mule deer)
  OK_occ <- as.data.frame(cbind(md_smr_predict_occ_sgnf, md_wtr_predict_occ_sgnf)) 
  colnames(OK_occ) <- c("MD_smr_occ", "MD_wtr_occ")
  
  #'  Northeast-only data (elk & white-tailed deer)
  NE_occ <- as.data.frame(cbind(elk_smr_predict_occ_sgnf, elk_wtr_predict_occ_sgnf, 
                                wtd_smr_predict_occ_sgnf, wtd_wtr_predict_occ_sgnf))  
  colnames(NE_occ) <- c("ELK_smr_occ", "ELK_wtr_occ", "WTD_smr_occ", "WTD_wtr_occ") 
  
  #'  Merge ungulate & predator data by study area
  Predicted_occ_OK <- Predicted_occ[Predicted_occ$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area data frame
    mutate(
      ELK_smr_occ = NA,
      ELK_wtr_occ = NA,
      WTD_smr_occ = NA,
      WTD_wtr_occ = NA
    ) %>%
    cbind(OK_occ) 
  Predicted_occ_NE <- Predicted_occ[Predicted_occ$Area == "Northeast",] %>%
    cbind(NE_occ) %>%
    #'  Need to account for columns that are present in other study area data frame
    mutate(
      MD_smr_occ = NA,
      MD_wtr_occ = NA
    )
  
  #'  Merge NE and OK predictions together
  Predicted_occ <- as.data.frame(rbind(Predicted_occ_NE, Predicted_occ_OK)) 

  
  ####  Predict relative probability of selection across study areas  ####
  #'  --------------------------------------------------------------------------
  #'  Manipulate RSF result tables
  rsf_coefs_signif <- rsf_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept, 
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDen) %>% 
    mutate(
      B.elev = ifelse(is.na(B.elev), 0, B.elev),
      B.slope = ifelse(is.na(B.slope), 0, B.slope),
      B.for = ifelse(is.na(B.for), 0, B.for),
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.rd = ifelse(is.na(B.rd), 0, B.rd)
    )
  
  #'  Function to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  #'  NOTE: Because these are RSFs, not using a logit transformation like with a 
  #'  normal logistic regression. Instead, dropping the intercept from the model 
  #'  and just exponentiating the coeffs*covs (Fieberg et al. 2020)
  predict_rsf <- function(cov, coef) {
    predict_odds <- c()
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$B.elev*cov$Elev[i] + coef$B.slope*cov$Slope[i] + 
                               coef$B.for*cov$PercForMix[i] + coef$B.grass*cov$PercXGrass[i] + 
                               coef$B.shrub*cov$PercXShrub[i] + coef$B.rd*cov$RoadDen[i]) 
    }
    predict_rsf <- as.data.frame(predict_odds) 
    return(predict_rsf)
  }
  #'  Run estimated coefficients from RSFs through function to predict relative probability of selection
  #'  Use correct z-transformed covariates for each species- & season-specific RSF
  #'  e.g., zcovs[[2]] = covs standardized based on summer bobcat mean & SDs
  #'  e.g., zcovs[[9]] = covs standardized based on winter elk mean & SDs
  bob_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[2]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Summer",])
  bob_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[3]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Winter",])
  coug_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[4]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Summer",])
  coug_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[5]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Winter",])
  coy_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[6]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Summer",])
  coy_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[7]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Winter",])
  elk_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[8]][zcovs[[8]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Summer",])
  elk_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[9]][zcovs[[9]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Winter",])
  md_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[10]][zcovs[[10]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Summer",])
  md_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[11]][zcovs[[11]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Winter",])
  wtd_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[12]][zcovs[[12]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[13]][zcovs[[13]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Winter",])
  wolf_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[14]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[15]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predators
  Predicted_rsf <- as.data.frame(all_covs_1km) %>%
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(bob_smr_predict_rsf_sgnf, bob_wtr_predict_rsf_sgnf,
          coug_smr_predict_rsf_sgnf, coug_wtr_predict_rsf_sgnf, 
          coy_smr_predict_rsf_sgnf, coy_wtr_predict_rsf_sgnf,
          wolf_smr_predict_rsf_sgnf, wolf_wtr_predict_rsf_sgnf)
  #'  Make sure you have the order right when you change the names!
  colnames(Predicted_rsf) <- c("obs", "Area", "x", "y",  
                               "BOB_smr_rsf", "BOB_wtr_rsf", 
                               "COUG_smr_rsf", "COUG_wtr_rsf", 
                               "COY_smr_rsf", "COY_wtr_rsf",
                               "WOLF_smr_rsf", "WOLF_wtr_rsf")
  
  #'  Okanogan-only data (mule deer)
  OK_rsf <- as.data.frame(cbind(md_smr_predict_rsf_sgnf, md_wtr_predict_rsf_sgnf))  
  colnames(OK_rsf) <- c("MD_smr_rsf", "MD_wtr_rsf")  
  
  #'  Northeast-only data (elk & white-tailed deer)
  NE_rsf <- as.data.frame(cbind(elk_smr_predict_rsf_sgnf, elk_wtr_predict_rsf_sgnf,  
                                wtd_smr_predict_rsf_sgnf, wtd_wtr_predict_rsf_sgnf)) 
  colnames(NE_rsf) <- c("ELK_smr_rsf", "ELK_wtr_rsf", "WTD_smr_rsf", "WTD_wtr_rsf")  
  
  #'  Merge ungulate & predator data by study area
  Predicted_rsf_OK <- Predicted_rsf[Predicted_rsf$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area data frame
    mutate(
      ELK_smr_rsf = NA,
      ELK_wtr_rsf = NA,
      WTD_smr_rsf = NA,
      WTD_wtr_rsf = NA
    ) %>%
    cbind(OK_rsf) 
  Predicted_rsf_NE <- Predicted_rsf[Predicted_rsf$Area == "Northeast",] %>%
    cbind(NE_rsf) %>%
    #'  Need to account for columns that are present in other study area data frame
    mutate(
      MD_smr_rsf = NA,
      MD_wtr_rsf = NA
    )
  
  #'  Merge NE and OK predictions together
  Predicted_rsf <- as.data.frame(rbind(Predicted_rsf_NE, Predicted_rsf_OK)) 
  
  #'  Function to identify any outliers in predicted RSF values and set to value
  #'  at 99th percentile
  outliers <- function(predicted, title) { 
    #'  Summarize predicted values
    hist(predicted, breaks = 100, main = title)
    boxplot(predicted, main = title)
    #'  What value represents the 99th percentile in the predicted RSF values
    quant <- quantile(predicted, c(0.99), na.rm = TRUE)
    #'  Print that value and maximum prediction
    print(quant); print(max(predicted, na.rm = TRUE))
    #'  Identify the 1% most extreme values and set to 99th percentile value
    predicted <- as.data.frame(predicted) %>%
      mutate(outlier = ifelse(predicted > quant, "outlier", "not_outlier"),
             adjusted_rsf = ifelse(outlier == "outlier", quant, predicted))
    #'  How many predicted values are above the 99th percentile?
    outlier <- predicted[predicted$outlier == "outlier",]
    outlier <- filter(outlier, !is.na(outlier))
    print(nrow(outlier))
    #'  What percentage of pixels are being forced to 99th percentile value
    print(nrow(outlier)/nrow(predicted))

    #'  Return these adjusted RSF predictions
    adjusted_rsf <- predicted$adjusted_rsf
    return(adjusted_rsf)
  }
  #'  Identify outlier predictions
  Predicted_rsf$BOB_smr_rsf2 <- outliers(Predicted_rsf$BOB_smr_rsf, "Bobcat Summer RSF Predictions")
  Predicted_rsf$BOB_wtr_rsf2 <- outliers(Predicted_rsf$BOB_wtr_rsf, "Bobcat Winter RSF Predictions")
  Predicted_rsf$COUG_smr_rsf2 <- outliers(Predicted_rsf$COUG_smr_rsf, "Cougar Summer RSF Predictions")
  Predicted_rsf$COUG_wtr_rsf2 <- outliers(Predicted_rsf$COUG_wtr_rsf, "Cougar Winter RSF Predictions")
  Predicted_rsf$COY_smr_rsf2 <- outliers(Predicted_rsf$COY_smr_rsf, "Coyote Summer RSF Predictions")
  Predicted_rsf$COY_wtr_rsf2 <- outliers(Predicted_rsf$COY_wtr_rsf, "Coyote Winter RSF Predictions")
  Predicted_rsf$ELK_smr_rsf2 <- outliers(Predicted_rsf$ELK_smr_rsf, "Elk Summer RSF Predictions")
  Predicted_rsf$ELK_wtr_rsf2 <- outliers(Predicted_rsf$ELK_wtr_rsf, "Elk Winter RSF Predictions")
  Predicted_rsf$MD_smr_rsf2 <- outliers(Predicted_rsf$MD_smr_rsf, "Mule Deer Summer RSF Predictions")
  Predicted_rsf$MD_wtr_rsf2 <- outliers(Predicted_rsf$MD_wtr_rsf, "Mule Deer Winter RSF Predictions")
  Predicted_rsf$WTD_smr_rsf2 <- outliers(Predicted_rsf$WTD_smr_rsf, "White-tailed Deer Summer RSF Predictions")
  Predicted_rsf$WTD_wtr_rsf2 <- outliers(Predicted_rsf$WTD_wtr_rsf, "White-tailed Deer Winter RSF Predictions")
  Predicted_rsf$WOLF_smr_rsf2 <- outliers(Predicted_rsf$WOLF_smr_rsf, "Wolf Summer RSF Predictions")
  Predicted_rsf$WOLF_wtr_rsf2 <- outliers(Predicted_rsf$WOLF_wtr_rsf, "Wolf Winter RSF Predictions")

  
  ####  Calculate Correlations between OccMod & RSF Predictions ####
  #'  --------------------------------------------------------------------------
  #'  Merge predicted occupancy with predicted resource selection
  Predicted_Occ_RSF <- Predicted_occ %>%
    full_join(Predicted_rsf, by = c("obs", "Area", "x", "y"))
  
  #'  Evaluate correlation between predicted space use for each paired set of models
  predict_corr <- function(predict_occu, predict_rsfs) {
    #'  Identify the maximum value in the predicted RSFs
    m <- max(predict_rsfs, na.rm = TRUE)
    #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
    stand_rsf <- as.data.frame(predict_rsfs) %>%
      mutate(scaled_rsf = predict_rsfs/m)
    #'  Combine occupancy and re-scaled RSF values
    predicted <- as.data.frame(cbind(predict_occu, stand_rsf$scaled_rsf)) %>%
      mutate(CellID = seq(1:nrow(.)))
    colnames(predicted) <- c("Occ_predictions", "RSF_rs_predictions", "CellID")
    #'  Calculated correlation between occupancy and RSF predictions
    pred_corr <- cor(predicted$Occ_predictions, predicted$RSF_rs_predictions, use = "complete.obs")
    
    return(pred_corr)
  }
  #'  Run each set of occupancy and RSF predictions through function
  bob_smr_corr <- predict_corr(Predicted_occ$BOB_smr_occ, Predicted_rsf$BOB_smr_rsf2)  
  bob_wtr_corr <- predict_corr(Predicted_occ$BOB_wtr_occ, Predicted_rsf$BOB_wtr_rsf2)  
  coug_smr_corr <- predict_corr(Predicted_occ$COUG_smr_occ, Predicted_rsf$COUG_smr_rsf2)
  coug_wtr_corr <- predict_corr(Predicted_occ$COUG_wtr_occ, Predicted_rsf$COUG_wtr_rsf2)  
  coy_smr_corr <- predict_corr(Predicted_occ$COY_smr_occ, Predicted_rsf$COY_smr_rsf2)
  coy_wtr_corr <- predict_corr(Predicted_occ$COY_wtr_occ, Predicted_rsf$COY_wtr_rsf2)  
  elk_smr_corr <- predict_corr(Predicted_occ$ELK_smr_occ, Predicted_rsf$ELK_smr_rsf2)  
  elk_wtr_corr <- predict_corr(Predicted_occ$ELK_wtr_occ, Predicted_rsf$ELK_wtr_rsf2)  
  md_smr_corr <- predict_corr(Predicted_occ$MD_smr_occ, Predicted_rsf$MD_smr_rsf2)  
  md_wtr_corr <- predict_corr(Predicted_occ$MD_wtr_occ, Predicted_rsf$MD_wtr_rsf2)  
  wtd_smr_corr <- predict_corr(Predicted_occ$WTD_smr_occ, Predicted_rsf$WTD_smr_rsf2)  
  wtd_wtr_corr <- predict_corr(Predicted_occ$WTD_wtr_occ, Predicted_rsf$WTD_wtr_rsf2)  
  wolf_smr_corr <- predict_corr(Predicted_occ$WOLF_smr_occ, Predicted_rsf$WOLF_smr_rsf2)  
  wolf_wtr_corr <- predict_corr(Predicted_occ$WOLF_wtr_occ, Predicted_rsf$WOLF_wtr_rsf2)  
    
  #'  Wrangle results into a table
  spp <- rep(c("Bobcat", "Cougar", "Coyote", "Elk", "Mule Deer", "White-tailed Deer", "Wolf"), each = 2)
  season <- rep(c("Summer", "Winter"), 7)
  corr <- c(bob_smr_corr, bob_wtr_corr, coug_smr_corr, coug_wtr_corr, coy_smr_corr, 
            coy_wtr_corr, elk_smr_corr, elk_wtr_corr, md_smr_corr, md_wtr_corr, 
            wtd_smr_corr, wtd_wtr_corr, wolf_smr_corr, wolf_wtr_corr)
  corr <- as.data.frame(corr)
  corr_results <- as.data.frame(cbind(spp, season, corr)) %>%
    transmute(
      Species = spp,
      Season = season,
      Correlation = round(corr, digits = 2)
    ) %>%
    arrange(Species)
  
  ####  NOTE: These correlation coefficients differ from published values due to 
  ####  larger pixels used here (1km^2 pixels). Switch to 30m^2 pixels to reproduce
  ####  exact published results (requires massive computation time and power).
  
  #'  End!
  