  #'  R code associated with: 
  #'  
  #'  Bassing, S. B., M. Devivo, T. R. Ganz, B. N. Kertson, L. R. Prugh, T. Roussin, 
  #'  L. Satterfield, R. M. Windell, A. J. Wirsing, & B. Gardner 
  #'  
  #'  Are we telling the same story? Comparing inferences made from camera trap 
  #'  and telemetry data for wildlife monitoring 
  #'  --------------------------------------------------------------------------
  #'  Script to run resource selection function (RSF) analyses for deer, elk, 
  #'  cougars, wolves, coyotes, and bobcats for summer (July 1 - Sept 29, 2018 & 
  #'  2019) and winter (Dec 1 - Mar 1, 2018 - 2019 & 2019 - 2020), respectively.
  #'  Specific date periods mirror 13-week sampling window used in corresponding
  #'  occupancy models. Result tables created at end of script.
  #'  
  #'  MUST download Data folder from Dryad repository for script to work.
  #'  
  #'  GPS collar data collected as part of the Washington Predator-Prey Project:
  #'  https://predatorpreyproject.weebly.com/
  #'  Covariate data extracted from remotely sensed data.
  #'  --------------------------------------------------------------------------

  #'  Load packages for selecting available points
  library(tidyverse)
  library(car)
  library(lme4)
  
  #'  Load 1:20 ratio of used and available points for all species
  #'  Available points drawn for 2nd order resource selection
  load("./Data/Use-Available Data/RSF_used_avail_pts_noXY.RData")
  
  #'  Center & scale covariates 
  #'  Note: standardizing across all IDs & years, but separately by season & species
  spp_dataPrep <- function(locs){
    #'  Make categorical variables factors
    locs$ID <- as.factor(locs$ID)
    locs$Used <- as.factor(locs$Used)
    locs$Area <- as.factor(locs$Area)
    locs$Year <- as.factor(locs$Year)
    locs$Season <- as.factor(locs$Season)
    #'  Standardize continuous variables
    locs$Elev <- scale(locs$Elev)
    locs$Slope <- scale(locs$Slope)
    locs$RoadDen <- scale(locs$RoadDen)
    locs$PercForMix <- scale(locs$PercForMix)
    locs$PercXGrass <- scale(locs$PercXGrass)
    locs$PercXShrub <- scale(locs$PercXShrub)
    #'  Leave weights as is
    locs$w <- locs$w
    
    locs <- as.data.frame(locs)
  
    return(locs)
  }
  #'  Format season & species-specific data
  bobData_smr <- spp_dataPrep(rsf_pts[[1]][rsf_pts[[1]]$Season == "Summer18" | rsf_pts[[1]]$Season == "Summer19",])
  bobData_wtr <- spp_dataPrep(rsf_pts[[1]][rsf_pts[[1]]$Season == "Winter1819" | rsf_pts[[1]]$Season == "Winter1920",])
  cougData_smr <- spp_dataPrep(rsf_pts[[2]][rsf_pts[[2]]$Season == "Summer18" | rsf_pts[[2]]$Season == "Summer19",])
  cougData_wtr <- spp_dataPrep(rsf_pts[[2]][rsf_pts[[2]]$Season == "Winter1819" | rsf_pts[[2]]$Season == "Winter1920",])
  coyData_smr <- spp_dataPrep(rsf_pts[[3]][rsf_pts[[3]]$Season == "Summer18" | rsf_pts[[3]]$Season == "Summer19",])
  coyData_wtr <- spp_dataPrep(rsf_pts[[3]][rsf_pts[[3]]$Season == "Winter1819" | rsf_pts[[3]]$Season == "Winter1920",])
  elkData_smr <- spp_dataPrep(rsf_pts[[4]][rsf_pts[[4]]$Season == "Summer18" | rsf_pts[[4]]$Season == "Summer19",])
  elkData_wtr <- spp_dataPrep(rsf_pts[[4]][rsf_pts[[4]]$Season == "Winter1819" | rsf_pts[[4]]$Season == "Winter1920",])
  mdData_smr <- spp_dataPrep(rsf_pts[[5]][rsf_pts[[5]]$Season == "Summer18" | rsf_pts[[5]]$Season == "Summer19",])
  mdData_wtr <- spp_dataPrep(rsf_pts[[5]][rsf_pts[[5]]$Season == "Winter1819" | rsf_pts[[5]]$Season == "Winter1920",])
  wtdData_smr <- spp_dataPrep(rsf_pts[[6]][rsf_pts[[6]]$Season == "Summer18" | rsf_pts[[6]]$Season == "Summer19",])
  wtdData_wtr <- spp_dataPrep(rsf_pts[[6]][rsf_pts[[6]]$Season == "Winter1819" | rsf_pts[[6]]$Season == "Winter1920",])
  wolfData_smr <- spp_dataPrep(rsf_pts[[7]][rsf_pts[[7]]$Season == "Summer18" | rsf_pts[[7]]$Season == "Summer19",])
  wolfData_wtr <- spp_dataPrep(rsf_pts[[7]][rsf_pts[[7]]$Season == "Winter1819" | rsf_pts[[7]]$Season == "Winter1920",])
  
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    used <- dat[dat$Used == 1,]
    covs <- used[,c("Elev", "Slope", "PercForMix", "PercXGrass", "PercXShrub", "RoadDen")]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (bob_smr_corr <- cov_correlation(bobData_smr))
  (bob_wtr_corr <- cov_correlation(bobData_wtr)) 
  (coug_smr_corr <- cov_correlation(cougData_smr)) 
  (coug_wtr_corr <- cov_correlation(cougData_wtr)) 
  (coy_smr_corr <- cov_correlation(coyData_smr)) 
  (coy_wtr_corr <- cov_correlation(coyData_wtr)) 
  (elk_smr_corr <- cov_correlation(elkData_smr)) 
  (elk_wtr_corr <- cov_correlation(elkData_wtr)) 
  (md_smr_corr <- cov_correlation(mdData_smr)) 
  (md_wtr_corr <- cov_correlation(mdData_wtr))
  (wtd_smr_corr <- cov_correlation(wtdData_smr))
  (wtd_wtr_corr <- cov_correlation(wtdData_wtr))
  (wolf_smr_corr <- cov_correlation(wolfData_smr))
  (wolf_wtr_corr <- cov_correlation(wolfData_wtr)) 

  
  #'  Resource Selection Functions (RSFs)
  #'  --------------------------------------------------------------------------
  #'  Function to run logistic mixed effects models that include a random effect 
  #'  for individual. Specific habitat covariates excluded from species- and 
  #'  season-specific models if there was high correlation (r = |0.7|) or they 
  #'  failed to converge in the corresponding occupancy model.
  #'  
  #'  NOTE: Deer and elk RSFs take a long time to run on a standard computer
  #'  --------------------------------------------------------------------------
  
  ####  Bobcat RSF  ####
  #'  SUMMERS 2018 & 2019
  bob_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                         data = bobData_smr, weights = w, family = binomial(link = "logit")) 
  summary(bob_s1819_rsf)
  car::vif(bob_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  bob_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                         data = bobData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(bob_w1820_rsf)
  car::vif(bob_w1820_rsf)
  
  
  ####  Cougar RSF  ####
  #'  SUMMERS 2018 & 2019
  coug_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                          data = cougData_smr, weights = w, family = binomial(link = "logit")) 
  summary(coug_s1819_rsf)
  car::vif(coug_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020
  coug_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                          data = cougData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(coug_w1820_rsf)
  car::vif(coug_w1820_rsf)
  
  
  ####  Coy RSF  ####
  #'  SUMMERS 2018 & 2019
  coy_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                         data = coyData_smr, weights = w, family = binomial(link = "logit")) 
  summary(coy_s1819_rsf)
  car::vif(coy_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  coy_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                         data = coyData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(coy_w1820_rsf)
  car::vif(coy_w1820_rsf)
 
  
  ####  Elk RSF  ####
  #'  SUMMERS 2018 & 2019
  #'  DROPPING RANDOM EFFECT (AND USING GLM) due to singularity issue 
  elk_s1819_rsf <- glm(Used ~ 1 + Elev + Slope + PercForMix + RoadDen, 
                       data = elkData_smr, weights = w, family = binomial(link = "logit")) 
  summary(elk_s1819_rsf)
  car::vif(elk_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020 
  elk_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),
                         data = elkData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(elk_w1820_rsf)
  car::vif(elk_w1820_rsf)
  
  
  ####  Mule Deer RSF  ####
  #'  SUMMERS 2018 & 2019
  md_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),
                        data = mdData_smr, weights = w, family = binomial(link = "logit")) 
  summary(md_s1819_rsf)
  car::vif(md_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020
  md_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID),
                        data = mdData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(md_w1820_rsf)
  car::vif(md_w1820_rsf)
  
  ####  White-tailed Deer RSF  ####
  #'  SUMMERS 2018 & 2019
  wtd_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),
                         data = wtdData_smr, weights = w, family = binomial(link = "logit"))
  summary(wtd_s1819_rsf)
  car::vif(wtd_s1819_rsf)

  #'  WINTERS 2018-2019 & 2019-2020
  wtd_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),
                         data = wtdData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(wtd_w1820_rsf)
  car::vif(wtd_w1820_rsf)

  
  ####  Wolf RSF  ####
  #'  SUMMERS 2018 & 2019
  wolf_s1819_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),
                          data = wolfData_smr, weights = w, family = binomial(link = "logit"))
  summary(wolf_s1819_rsf)
  car::vif(wolf_s1819_rsf)
  
  #'  WINTERS 2018-2019 & 2019-2020
  wolf_w1820_rsf <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),
                          data = wolfData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(wolf_w1820_rsf)
  car::vif(wolf_w1820_rsf)

  
  #'  List RSF outputs together
  rsf_outputs <- list(bob_s1819_rsf, bob_w1820_rsf, coug_s1819_rsf, coug_w1820_rsf, 
                      coy_s1819_rsf, coy_w1820_rsf, elk_s1819_rsf, elk_w1820_rsf, 
                      md_s1819_rsf, md_w1820_rsf, wtd_s1819_rsf, wtd_w1820_rsf, 
                      wolf_s1819_rsf, wolf_w1820_rsf)
  
  #'  Save RSF outputs because they take forever to run!
  save(rsf_outputs, "./Data/Use-Available Data/RSF_outputs.RData")

  
  ####  Summary tables  ####
  #'  --------------------------------------------------------------------------
  #'  Load RSF outputs data if needed
  load("./Data/Use-Available Data/RSF_outputs.RData")
  
  #'  Function to extract parameter estimates & p-values from rsf outputs
  rounddig <- 2
  rsf_out <- function(mod, spp, season){
    betas <- summary(mod)$coef[,1]
    se <- summary(mod)$coef[,2]
    z <- summary(mod)$coef[,3]
    pval <- summary(mod)$coef[,4]
    out <- as.data.frame(cbind(betas, se, pval)) %>%
      transmute(
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Parameter = row.names(.),
        Estimate = round(betas, rounddig),
        SE = round(se, rounddig),
        Z = round(z, rounddig),
        Pval = round(pval, rounddig)) 
    rownames(out) <- NULL
    return(out)
  }
  #'  Create species- and season-specific results tables
  bob_smr_rsf_tbl <- rsf_out(rsf_outputs[[1]], "Bobcat", "Summer")
  bob_wtr_rsf_tbl <- rsf_out(rsf_outputs[[2]], "Bobcat", "Winter")
  coug_smr_rsf_tbl <- rsf_out(rsf_outputs[[3]], "Cougar", "Summer")
  coug_wtr_rsf_tbl <- rsf_out(rsf_outputs[[4]], "Cougar", "Winter")
  coy_smr_rsf_tbl <- rsf_out(rsf_outputs[[5]], "Coyote", "Summer")
  coy_wtr_rsf_tbl <- rsf_out(rsf_outputs[[6]], "Coyote", "Winter")
  elk_smr_rsf_tbl <- rsf_out(rsf_outputs[[7]], "Elk", "Summer")
  elk_wtr_rsf_tbl <- rsf_out(rsf_outputs[[8]], "Elk", "Winter")
  md_smr_rsf_tbl <- rsf_out(rsf_outputs[[9]], "Mule Deer", "Summer")
  md_wtr_rsf_tbl <- rsf_out(rsf_outputs[[10]], "Mule Deer", "Winter")
  wtd_smr_rsf_tbl <- rsf_out(rsf_outputs[[11]], "White-tailed Deer", "Summer")
  wtd_wtr_rsf_tbl <- rsf_out(rsf_outputs[[12]], "White-tailed Deer", "Winter")
  wolf_smr_rsf_tbl <- rsf_out(rsf_outputs[[13]], "Wolf", "Summer")
  wolf_wtr_rsf_tbl <- rsf_out(rsf_outputs[[14]], "Wolf", "Winter")
  
  #'  Merge into larger data frames for easy comparison
  summer_rsf_tbl <- rbind(bob_smr_rsf_tbl, coug_smr_rsf_tbl, coy_smr_rsf_tbl, 
                          elk_smr_rsf_tbl, md_smr_rsf_tbl, wtd_smr_rsf_tbl, 
                          wolf_smr_rsf_tbl) 
  winter_rsf_tbl <- rbind(bob_wtr_rsf_tbl, coug_wtr_rsf_tbl, coy_wtr_rsf_tbl, 
                          elk_wtr_rsf_tbl, md_wtr_rsf_tbl, wtd_wtr_rsf_tbl, 
                          wolf_wtr_rsf_tbl) 
  #'  One big table
  rsf_results <- rbind(summer_rsf_tbl, winter_rsf_tbl) %>%
    arrange(Species)
  colnames(rsf_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  
  #'  Save for predicting relative prob. of resource selection across study area
  write.csv(rsf_results, "./Data/Results tables/RSF_Results.csv")  
  
  #'  End
  