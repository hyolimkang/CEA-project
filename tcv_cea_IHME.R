# TCV cost-effectiveness

# things to do 
# update to latin hypercube sampling 
# ggplot for CEA Acceptability curve
# ggplot for tornado plot
# scenario analysis (how many scenarios?)
# correction factor
# any adjustments?
# https://github.com/pierucci/heemod/blob/master/R/acceptability_curve.R


# load libraries
library (dampack)
library (data.table)
library (ggplot2)
library (tictoc)
library (caret)
library (rriskDistributions)
library (BCEA)
library (writexl)
library (tidyverse)
library (lhs)
library (truncnorm)

# clear workspace
rm (list = ls())
tic ()

# ------------------------------------------------------------------------------
# function
# create psa samples for input parameters
# ------------------------------------------------------------------------------
create_psa_sample <- function (sample_n) {
  
  # parameters
  my_params <-c ("v_cost",               # vaccine cost (including syringes): 2018 USD
                 "delivery_cost",        # delivery cost: 2018 USD
                 "tot_pop",              # total population in Phase 1 
                 "vacc_pop",             # total number of vaccinated people
                 "postvacc_p1_ipd",      # post vaccination cases of IPD in Phase1 
                 "postvacc_p1_opd",      # post vaccination cases of OPD in Phase1 
                 "postvacc_p1",          # post vaccination cases (all)
                 "facility_cost",        # health facility cost of OPD : 2021 USD
                 "dmc_ipd",              # direct medical cost of IPD : 2021 USD
                 "dmc_opd",              # direct medical cost of OPD : 2021 USD
                 "dmc",                  # average dmc (for unvaccinated)
                 "dnmc_ipd",             # direct non-medical cost of IPD : 2021 USD
                 "dnmc_opd",             # direct non-medical csot of OPD : 2021 USD
                 "dnmc",                 # average dnmc (for unvaccinated)
                 "indirect_ipd",         # indirect cost of IPD : 2021 USD
                 "indirect_opd",         # indirect cost of OPD : 2021 USD
                 "indirect",             # average indirect cost (for unvaccinated)
                 "ve",                   # vaccine-effectiveness
                 "dw_ipd",               # disability weight of IPD
                 "dw_opd",               # disability weight of OPD (severe case)
                 "illness_duration_ipd", # illness duration of IPD: days
                 "illness_duration_opd", # illness duration of OPD: days 
                 "cfr_ipd",              # case fatality rate of IPD
                 "cfr_opd",              # case fatality rate of OPD
                 "age_death",            # average age at death
                 "life_exp"   )          # reference life expectancy
  
  
  # constant value
  # "postvacc"        # post-vaccination case
  # "vacc_pop"
  # "postvacc_p1"
  
  # trickle down uncertainty
  # "prevacc",         # pre-vaccination case
  # "totcost_vacc",    # total cost of vaccinated group
  # "totcost_unvacc")  # total cost of unvaccinated group
  
  
  # distributions for parameters
  my_dists <- c("constant",
                "gamma",
                "constant",
                "constant",
                "constant",
                "constant",
                "constant",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "gamma",
                "truncated-normal",
                "truncated-normal",
                "truncated-normal",
                "gamma",
                "gamma",
                "truncated-normal",
                "truncated-normal",
                "normal",
                "normal")
  
  
  # parameter types
  my_parameterization_types <-c ("val",
                                 "mean, sd",
                                 "val",
                                 "val",
                                 "val",
                                 "val",
                                 "val",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul",
                                 "mean, sd",
                                 "mean, sd")
  
  
  # values for parameters
  my_dists_params <- list( c (2.96),
                           c (1.49, 1.38),
                           c (159831),
                           c (113420), 
                           c (200),
                           c (459) ,
                           c (659),
                           c (97.33,  85.33),
                           c (234.77, 265.99),
                           c (115.38, 149.92),
                           c (183.34, 229.85),
                           c (48.59,  37.40),
                           c (21.49,  25.90),
                           c (36.91,  35.39),
                           c (77.54,  29.51),
                           c (74.29,  26.90),
                           c (76.13,  28.25),
                           c (.81,  2.05, .54, .926),
                           c (.21,  1.37, .19, .23),
                           c (.052,  1.37, .031, .079),
                           c (.040, .0153),
                           c (.046, .025),
                           c (.028, .211, .02, .036),
                           c (.019, .72, 0.01, 0.045),
                           c (7.56, 6.53),
                           c (72.68, 3.07))
  
  # input samples for PSA
  my_psa <- gen_psa_samp (params                 = my_params,
                          dists                  = my_dists,
                          parameterization_types = my_parameterization_types,
                          dists_params           = my_dists_params,
                          n                      = sample_n)
  
  return (my_psa)
  
}

# ------------------------------------------------------------------------------
# Latin Hypercube Sampling
# ------------------------------------------------------------------------------

# https://stats.stackexchange.com/questions/495215/standard-error-standard-deviation-and-variance-confusion

set.seed (3)
runs <- 4000

# a design with n samples from k parameters
A <- randomLHS (n = runs, 
                k = 30) 
# It is common to transform the margins of the design (the columns) 
# into other distributions
# for gamma dist: a= shape (m^2/sigma^2), b= scale = 1/ (sigma^2/mean) 
# (m is mean and sigma is SD)
# for beta dist function: 

beta <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
# beta parameterization for VE
beta(0.816,  0.01982924)
# beta parameterization for DW
beta(0.21, 0.0001)
beta(0.052, 0.0001)
# beta parameterization for CFR
beta(0.028, 0.00007972)
beta(0.019, 0.00001666)


lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A))
lhs_sample [,1]  <- 2.93
lhs_sample [,2]  <- qgamma (p = A[,2], shape = ((1.52/0.520408163)^2), rate = (1.52/(0.520408163)^2),
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,3]  <- 0.3218391* 159831
lhs_sample [,4]  <- 0.3295019*159831
lhs_sample [,5]  <- 0.348659 * 159831
lhs_sample [,6]  <- 20405 *(113420/79836)
lhs_sample [,7]  <- 28911*(113420/79836)
lhs_sample [,8]  <- 30520*(113420/79836)
lhs_sample [,9]  <- qunif (p = A[,9], min =128.705, max = 588.66, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,10]  <- qunif (p = A[,10], min = 429.7415, max = 1772.195, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,11]  <- qunif (p = A[,11], min = 438.77, max = 1766.42, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,12]  <- qgamma (p = A[,12], shape = ((97.33/12.33)^2), rate = (97.33/(12.33)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,13]  <- qgamma (p = A[,13], shape = ((234.7688/265.9934)^2), rate = (234.7688/(265.9934)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,14] <- qgamma (p = A[,14], shape = ((115.3841/149.9214)^2), rate = (115.3841/(149.9214)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,15] <- qgamma (p = A[,15], shape = ((183.3415/229.8478)^2), rate = (183.3415/(229.8478)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,16] <- qgamma (p = A[,16], shape = ((48.58656/37.39749)^2), rate = (48.58656/(37.39749)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,17] <- qgamma (p = A[,17],  shape = ((21.48999/25.8971)^2), rate = (21.48999/(25.8971)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,18] <- qgamma (p = A[,18],  shape = ((36.91419/35.39041)^2), rate = (36.91419/(35.39041)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,19] <- qgamma (p = A[,19],  shape = ((77.53769/29.51196)^2), rate = (77.53769/(29.51196)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,20] <- qgamma (p = A[,20],   shape = ((74.28838/26.90198)^2), rate = (74.28838/(26.90198)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,21] <- qgamma (p = A[,21],  shape = ((76.13798/28.24641)^2),rate = (76.13798/(28.24641)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,22] <- qbeta  (p = A[,22], shape1 = 5.362628, shape2 = 1.20922, ncp=0,lower.tail = TRUE, log.p = FALSE)
lhs_sample [,23] <- qbeta (p = A[,23], shape1 = 348.18, shape2 = 1309.82, ncp=0, lower.tail = TRUE, log.p = FALSE) 
lhs_sample [,24] <- qbeta (p = A[,24], shape1 = 25.58192, shape2 = 466.3781, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,25] <- qgamma (p = A[,25], shape = ((0.0404/0.01532)^2), rate = (0.0404/(0.01532)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,26] <- qgamma (p = A[,26], shape = ((0.04559688/ 0.02500143)^2), rate = (0.04559688/(0.02500143)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,27] <- qbeta (p = A[,27], shape1 = 9.531057, shape2 = 330.8638, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,28] <- qbeta (p = A[,28], shape1 =  21.23796, shape2 = 1096.55, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,29] <- qgamma (p = A[,29], shape = ((7.56/0.527)^2), rate = (7.56/(0.527)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,30] <- qnorm (p = A[,30], mean = 72.68, sd = 1.73979, lower.tail = TRUE, log.p = FALSE)

# change the name of the lhs_sample
cols <- c(  "v_cost", 
            "delivery_cost",    
            "tot_pop_age1",
            "tot_pop_age2",
            "tot_pop_age3",
            "vacc_pop_age1",  
            "vacc_pop_age2",
            "vacc_pop_age3",
            "prevacc_age1",      
            "prevacc_age2",
            "prevacc_age3",
            "facility_cost",        
            "dmc_ipd",              
            "dmc_opd",         
            "dmc",               
            "dnmc_ipd",          
            "dnmc_opd",             
            "dnmc",                 
            "indirect_ipd",      
            "indirect_opd",        
            "indirect",             
            "ve",                   
            "dw_ipd",               
            "dw_opd",               
            "illness_duration_ipd", 
            "illness_duration_opd", 
            "cfr_ipd",             
            "cfr_opd",              
            "age_death",            
            "life_exp" ) 
colnames (lhs_sample) <- cols
# change the matrix into the data table. (R doesn't recognize atomic values inside matrix)
lhs_sample <- as.data.table (lhs_sample)

# function - computer icer for a given input parameter sample

compute_icer <- function (v_cost            = 2.93,
                          delivery_cost,
                          tot_pop_age1      = 0.3218391* 159831,
                          tot_pop_age2      = 0.3295019*159831,
                          tot_pop_age3      = 0.348659 * 159831,
                          vacc_pop_age1     = 20405 *(113420/79836),
                          vacc_pop_age2     = 28911*(113420/79836),
                          vacc_pop_age3     = 30520*(113420/79836),
                          prevacc_age1, 
                          prevacc_age2,
                          prevacc_age3,
                          facility_cost,
                          dmc_ipd,
                          dmc_opd,
                          dmc,
                          dnmc_ipd,
                          dnmc_opd,
                          dnmc,
                          indirect_ipd,
                          indirect_opd,
                          indirect,
                          ve,
                          dw_ipd,
                          dw_opd,
                          illness_duration_ipd,
                          illness_duration_opd,
                          cfr_ipd,
                          cfr_opd,
                          age_death,
                          life_exp) {

  # prevacc age-group1 (IPD)
  
  prevacc_age1_ipd <- prevacc_age1 * 0.011
  
  # prevacc age-group1 (OPD)
  
  prevacc_age1_opd <- prevacc_age1 - prevacc_age1_ipd
  
  # prevacc age-group2 (IPD)
  
  prevacc_age2_ipd <- prevacc_age2 * 0.011
  
  # prevacc age-group2 (OPD)
  
  prevacc_age2_opd <- prevacc_age2 - prevacc_age2_ipd
  
  # prevacc age-group3 (IPD)
  
  prevacc_age3_ipd <- prevacc_age3 * 0.011
  
  # prevacc age-group3 (OPD)
  
  prevacc_age3_opd <- prevacc_age3 - prevacc_age3_ipd
  
  # post-vacc age-group1 
  
  postvacc_age1 <- (prevacc_age1)*(1 - ve*(vacc_pop_age1/tot_pop_age1))
  
  # post-vacc age-group2 
  
  postvacc_age2 <- (prevacc_age2)*(1 - ve*(vacc_pop_age2/tot_pop_age2))
  
  # post-vacc age-group3 
  
  postvacc_age3 <- (prevacc_age3)*(1 - ve*(vacc_pop_age3/tot_pop_age3))
  
  # post-vacc age-group1 (IPD)
  
  postvacc_age1_ipd <- (prevacc_age1_ipd)*(1 - ve*(vacc_pop_age1/tot_pop_age1))
  
  # post-vacc age-group1 (OPD)
  
  postvacc_age1_opd <- (prevacc_age1_opd)*(1 - ve*(vacc_pop_age1/tot_pop_age1))
  
  # post-vacc age-group2 (IPD)
  
  postvacc_age2_ipd <- (prevacc_age2_ipd)*(1 - ve*(vacc_pop_age2/tot_pop_age2))
  
  # post-vacc age-group2 (OPD)
  
  postvacc_age2_opd <- (prevacc_age2_opd)*(1 - ve*(vacc_pop_age2/tot_pop_age2))
  
  # post-vacc age-group3 (IPD)
  
  postvacc_age3_ipd <- (prevacc_age3_ipd)*(1 - ve*(vacc_pop_age3/tot_pop_age3))
  
  # post-vacc age-group3 (OPD)
  
  postvacc_age3_opd <- (prevacc_age3_opd)*(1 - ve*(vacc_pop_age3/tot_pop_age3))
  
  # case-averted age-group1 (per 100 000)
  
  case_avert_age1 <- (prevacc_age1 - postvacc_age1)
  
  # case-averted age-group2 (per 100 000)
  
  case_avert_age2 <- (prevacc_age2 - postvacc_age2)
  
  # case-averted age-group3 (per 100 000)
  
  case_avert_age3 <- (prevacc_age3 - postvacc_age3)
  
  # case-averted age-group1 IPD (per 100 000)
  
  case_avert_age1_ipd <- (case_avert_age1) * 0.011
  
  # case-averted age-group1 OPD (per 100 000)
  
  case_avert_age1_opd <- (case_avert_age1) - (case_avert_age1_ipd)
  
  # case-averted age-group2 IPD (per 100 000)
  
  case_avert_age2_ipd <- (case_avert_age2) * 0.011
  
  # case-averted age-group2 OPD (per 100 000)
  
  case_avert_age2_opd <- (case_avert_age2) - (case_avert_age2_ipd)
  
  # case-averted age-group3 IPD (per 100 000)
  
  case_avert_age3_ipd <- (case_avert_age3) * 0.011
  
  # case-averted age-group3 OPD (per 100 000)
  
  case_avert_age3_opd <- (case_avert_age3) - (case_avert_age3_ipd)
  
  # incremental cost of inpatients 
  
  incremental_cost_ipd   <- (v_cost + delivery_cost) * (vacc_pop_age1 + vacc_pop_age2 + vacc_pop_age3) - (case_avert_age1_ipd + case_avert_age2_ipd + case_avert_age3_ipd)*(facility_cost + dmc_ipd + dnmc_ipd)
  
  # incremental cost of outpatients
  
  incremental_cost_opd   <- (v_cost + delivery_cost) * (vacc_pop_age1 + vacc_pop_age2 + vacc_pop_age3) - (case_avert_age1_opd + case_avert_age2_opd + case_avert_age3_opd)*(facility_cost + dmc_opd + dnmc_opd)
  
  # incremental cost total (used for the final value as a delta C)
  
  incremental_cost_total <-  (case_avert_age1 + case_avert_age2 + case_avert_age3)*(facility_cost + dmc +dnmc + indirect) - (v_cost + delivery_cost)* (vacc_pop_age1 + vacc_pop_age2 + vacc_pop_age3)
  
  # total averted cases (all age group)
  
  case_avert_total      <- (case_avert_age1 + case_avert_age2 + case_avert_age3)
  
  # YLL calculation component: amount of years of life loss between life-expectancy and average age at death
  
  distance               <- (life_exp) - (age_death)
  
  # number of deaths in no-vaccination situation (inpatients)
  
  pre_death_ipd          <- (cfr_ipd)*(prevacc_age1_ipd + prevacc_age2_ipd + prevacc_age3_ipd)
  
  pre_death_age1_ipd     <- (cfr_ipd)*(prevacc_age1_ipd)
  
  pre_death_age2_ipd     <- (cfr_ipd)*(prevacc_age2_ipd)
  
  pre_death_age3_ipd     <- (cfr_ipd)*(prevacc_age3_ipd) 
  
  # number of deaths in no-vaccination situation (outpatients)
  
  pre_death_opd          <- (cfr_opd)*(prevacc_age1_opd + prevacc_age2_opd + prevacc_age3_opd)
  
  pre_death_age1_opd     <- (cfr_opd)*(prevacc_age1_opd)
  
  pre_death_age2_opd     <- (cfr_opd)*(prevacc_age2_opd)
  
  pre_death_age3_opd     <- (cfr_opd)*(prevacc_age3_opd)
  
  # number of deaths in no-vaccination situation (total)
  
  pre_death_total        <- (pre_death_ipd) + (pre_death_opd)
  
  pre_death_age1         <- (pre_death_age1_ipd) + (pre_death_age1_opd)
  
  pre_death_age2         <- (pre_death_age2_ipd) + (pre_death_age2_opd)
  
  pre_death_age3         <- (pre_death_age3_ipd) + (pre_death_age3_opd)
  
  # number of deaths post-vaccination situation (inpatients)
  
  post_death_ipd         <- (cfr_ipd)*(postvacc_age1_ipd + postvacc_age2_ipd + postvacc_age3_ipd)
  
  post_death_age1_ipd    <- (cfr_ipd)*(postvacc_age1_ipd)
  
  post_death_age2_ipd    <- (cfr_ipd)*(postvacc_age2_ipd)

  post_death_age3_ipd    <- (cfr_ipd)*(postvacc_age3_ipd)
  
  # number of deaths post-vaccination situation (outpatients)
  
  post_death_opd         <- (cfr_opd)*(postvacc_age1_opd + postvacc_age2_opd + postvacc_age3_opd)
  
  post_death_age1_opd    <- (cfr_opd)*(postvacc_age1_opd)
  
  post_death_age2_opd    <- (cfr_opd)*(postvacc_age2_opd)

  post_death_age3_opd    <- (cfr_opd)*(postvacc_age3_opd)
  
  # number of deaths post-vaccination situation (total)
  
  post_death_total       <-  (post_death_ipd) + (post_death_opd) 
  
  post_death_age1        <-  (post_death_age1_ipd) + (post_death_age1_opd)
  
  post_death_age2        <-  (post_death_age2_ipd) + (post_death_age2_opd)
  
  post_death_age3        <-  (post_death_age3_ipd) + (post_death_age3_opd)
  
  # number of deaths averted (total)
  
  death_averted          <-  (pre_death_total) - (post_death_total)
  
  death_averted_age1     <-  ((pre_death_age1_ipd) - (post_death_age1_ipd)) + ((pre_death_age1_opd) - (post_death_age1_opd))
  
  death_averted_age2     <-  ((pre_death_age2_ipd) - (post_death_age2_ipd)) + ((pre_death_age2_opd) - (post_death_age2_opd))
  
  death_averted_age3     <-  ((pre_death_age3_ipd) - (post_death_age3_ipd)) + ((pre_death_age3_opd) - (post_death_age3_opd))
  
  # number of deaths averted (inpatients)
  
  death_averted_ipd      <-  (pre_death_ipd) - (post_death_ipd)
  
  # number of deaths averted (outpatients)
  
  death_averted_opd      <-  (pre_death_opd) - (post_death_opd)
  
  # YLL of inpatients in the pre-vaccination situation (total number of inpatient deaths * amount of years loss per person)
  
  yll_pre_ipd            <- (pre_death_ipd)*(distance)
  
  # YLL of outpatients in the pre-vaccination situation (total number of outpatient deaths * amount of years loss per person)
  
  yll_pre_opd            <- (pre_death_opd)*(distance)
  
  # YLL total in pre-vaccination 
  
  yll_pre_total          <- (yll_pre_ipd)  + (yll_pre_opd)
  
  # amount of years loss for inpatients in post-vaccination
  
  yll_post_ipd           <- (post_death_ipd)*(distance)
  
  # amount of years loss for outpatients in post-vaccination
  
  yll_post_opd           <- (post_death_opd)*(distance)
  
  # total amount of years loss in post-vaccination
  
  yll_post_total         <- (yll_post_ipd) + (yll_post_opd)
  
  # years of life with disease for inpatients in pre-vacc situation (number of inpatient cases * disability weight * total illness duration presented as years)
  
  yld_pre_ipd            <- (prevacc_age1_ipd + prevacc_age2_ipd + prevacc_age3_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  # years of life with disease for outpatients in pre-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
  
  yld_pre_opd            <- (prevacc_age1_opd + prevacc_age2_opd + prevacc_age3_opd)*(dw_opd)*(illness_duration_opd)
  
  # total YLD of pre-vaccination typhoid cases
  
  yld_pre_total          <- (yld_pre_ipd) + (yld_pre_opd)
  
  # total DALYs for pre-vaccination typhoid cases (YLL total + YLD total)
  
  dalys_pre_total        <- (yll_pre_total) + (yld_pre_total)
  
  # years of life with disease for inpatients in post-vacc situation (number of inpatient cases * disability weight * total illness duration presented as years)
  
  yld_post_ipd           <- (postvacc_age1_ipd + postvacc_age2_ipd + postvacc_age3_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  # years of life with disease for outpatients in post-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
  
  yld_post_opd           <- (postvacc_age1_opd + postvacc_age2_opd + postvacc_age3_opd)*(dw_opd)*(illness_duration_opd)
  
  # total YLD for post-vaccination typhoid cases (YLL total + YLD total)
  
  yld_post_total         <- (yld_post_ipd) + (yld_post_opd)
  
  # total DALYs for post-vaccination typhoid cases (YLL total + YLD total)
  
  dalys_post_total       <- (yll_post_total) + (yld_post_total)
  
  # YLD averted for inpatients
  
  yld_averted_ipd        <- (case_avert_age1_ipd + case_avert_age2_ipd + case_avert_age3_ipd)* (dw_ipd) *(illness_duration_ipd)
  
  # YLD averted for outpatients
  
  yld_averted_opd        <- (case_avert_age1_opd + case_avert_age2_opd + case_avert_age3_opd)* (dw_opd) *(illness_duration_opd)
  
  # total YLD averted
  
  yld_averted_total      <- (yld_averted_ipd) + (yld_averted_opd) 
  
  # YLL averted  for inpatients
  
  yll_averted_ipd        <- (yll_pre_ipd)  -(yll_post_ipd)
  
  # YLL averted for outpatients
  
  yll_averted_opd        <- (yll_pre_opd)  -(yll_post_opd)
  
  # total YLL averted
  
  yll_averted_total      <- (yll_averted_ipd) + (yll_averted_opd)
  
  # DALYs averted for inpatients
  
  incremental_daly_ipd   <- (yld_averted_ipd) + (yll_averted_ipd)
  
  # DALYs averted for outpatients
  
  incremental_daly_opd   <- (yld_averted_opd) + (yll_averted_opd)
  
  # total DALYs averted (used for final ICER value in the paper: delta E as a denominator)
  
  incremental_daly_total <- (dalys_pre_total) - (dalys_post_total) 
  
  # ICER for inpatients presented as cost-per-DALY averted
  
  icer_daly_ipd          <- incremental_cost_ipd / incremental_daly_ipd
  
  # ICER for outpatients presented as cost-per-DALY averted
  
  icer_daly_opd          <- incremental_cost_opd / incremental_daly_opd
  
  # ICER presented as cost-per-DALY averted (used in the paper as a final value)
  
  icer_daly_total        <- incremental_cost_total / incremental_daly_total
  
  
  
  return (list (prevacc_age1_ipd              = prevacc_age1_ipd,
                prevacc_age1_opd              = prevacc_age1_opd,
                prevacc_age2_ipd              = prevacc_age2_ipd,
                prevacc_age2_opd              = prevacc_age2_opd,
                prevacc_age3_ipd              = prevacc_age3_ipd,
                prevacc_age3_opd              = prevacc_age3_opd,
                postvacc_age1                 = postvacc_age1,
                postvacc_age2                 = postvacc_age2,
                postvacc_age3                 = postvacc_age3,
                postvacc_age1_ipd             = postvacc_age1_ipd,
                postvacc_age1_opd             = postvacc_age1_opd,
                postvacc_age2_ipd             = postvacc_age2_ipd,
                postvacc_age2_opd             = postvacc_age2_opd,
                postvacc_age3_ipd             = postvacc_age3_ipd,
                postvacc_age3_opd             = postvacc_age3_opd,
                case_avert_age1               = case_avert_age1,
                case_avert_age2               = case_avert_age2,
                case_avert_age3               = case_avert_age3,
                case_avert_age1_ipd           = case_avert_age1_ipd,
                case_avert_age1_opd           = case_avert_age1_opd,
                case_avert_age2_ipd           = case_avert_age2_ipd,
                case_avert_age2_opd           = case_avert_age2_opd,
                case_avert_age3_ipd           = case_avert_age3_ipd,
                case_avert_age3_opd           = case_avert_age3_opd,
                incremental_cost_ipd          = incremental_cost_ipd,
                incremental_cost_opd          = incremental_cost_opd,
                incremental_cost_total        = incremental_cost_total,
                case_avert_total              = case_avert_total,
                pre_death_ipd                 = pre_death_ipd,
                pre_death_opd                 = pre_death_opd,
                pre_death_total               = pre_death_total,
                pre_death_age1                = pre_death_age1,
                pre_death_age2                = pre_death_age2,
                pre_death_age3                = pre_death_age3,
                post_death_ipd                = post_death_ipd,
                post_death_opd                = post_death_opd,
                post_death_total              = post_death_total,
                post_death_age1               = post_death_age1,
                post_death_age2               = post_death_age2,
                post_death_age3               = post_death_age3,
                death_averted                 = death_averted,
                death_averted_ipd             = death_averted_ipd,
                death_averted_opd             = death_averted_opd,
                death_averted_age1            = death_averted_age1,
                death_averted_age2            = death_averted_age2,
                death_averted_age3            = death_averted_age3,
                dalys_pre_total               = dalys_pre_total,
                dalys_post_total              = dalys_post_total,
                incremental_daly_total        = incremental_daly_total,
                icer_daly_ipd                 = icer_daly_ipd,
                icer_daly_opd                 = icer_daly_opd,
                icer_daly_total               = icer_daly_total)) 
}

# start time
tic ()
print (Sys.time ())

# create empty icer data table 
icer_dt <- data.table (run_id = 1:runs)

# combine data tables -- run_id & lhs sample
icer_dt <- cbind (icer_dt, lhs_sample)

# loop through psa sample to generate icers
for (i in 1:runs) {
  
  icer_sample <- compute_icer (lhs_sample$v_cost                    [i],
                               lhs_sample$delivery_cost             [i],
                               lhs_sample$tot_pop_age1              [i],
                               lhs_sample$tot_pop_age2              [i],
                               lhs_sample$tot_pop_age3              [i],
                               lhs_sample$vacc_pop_age1             [i],
                               lhs_sample$vacc_pop_age2             [i],
                               lhs_sample$vacc_pop_age3             [i],
                               lhs_sample$prevacc_age1              [i],
                               lhs_sample$prevacc_age2              [i],
                               lhs_sample$prevacc_age3              [i],
                               lhs_sample$facility_cost             [i],   
                               lhs_sample$dmc_ipd                   [i],
                               lhs_sample$dmc_opd                   [i],
                               lhs_sample$dmc                       [i],
                               lhs_sample$dnmc_ipd                  [i],
                               lhs_sample$dnmc_opd                  [i],
                               lhs_sample$dnmc                      [i],
                               lhs_sample$indirect_ipd              [i],
                               lhs_sample$indirect_opd              [i],
                               lhs_sample$indirect                  [i],
                               lhs_sample$ve                        [i],
                               lhs_sample$dw_ipd                    [i],
                               lhs_sample$dw_opd                    [i],
                               lhs_sample$illness_duration_ipd      [i],
                               lhs_sample$illness_duration_opd      [i],
                               lhs_sample$cfr_ipd                   [i],
                               lhs_sample$cfr_opd                   [i],
                               lhs_sample$age_death                 [i],
                               lhs_sample$life_exp                  [i])
  
  icer_dt [i, `:=` ( incremental_cost_ipd            = icer_sample$incremental_cost_ipd,
                     incremental_cost_opd            = icer_sample$incremental_cost_opd,
                     incremental_cost_total          = icer_sample$incremental_cost_total,
                     prevacc_age1_ipd                = icer_sample$prevacc_age1_ipd,
                     prevacc_age2_ipd                = icer_sample$prevacc_age2_ipd,
                     prevacc_age3_ipd                = icer_sample$prevacc_age3_ipd,
                     prevacc_age1_opd                = icer_sample$prevacc_age1_opd,
                     prevacc_age2_opd                = icer_sample$prevacc_age2_opd,
                     prevacc_age3_opd                = icer_sample$prevacc_age3_opd,
                     postvacc_age1                   = icer_sample$postvacc_age1,
                     postvacc_age2                   = icer_sample$postvacc_age2, 
                     postvacc_age3                   = icer_sample$postvacc_age3,
                     postvacc_age1_ipd               = icer_sample$postvacc_age1_ipd,
                     postvacc_age2_ipd               = icer_sample$postvacc_age2_ipd,
                     postvacc_age3_ipd               = icer_sample$postvacc_age3_ipd,
                     postvacc_age1_opd               = icer_sample$postvacc_age1_opd,
                     postvacc_age2_opd               = icer_sample$postvacc_age2_opd,
                     postvacc_age3_opd               = icer_sample$postvacc_age3_opd,
                     case_avert_age1                 = icer_sample$case_avert_age1,
                     case_avert_age2                 = icer_sample$case_avert_age2,
                     case_avert_age3                 = icer_sample$case_avert_age3,
                     case_avert_age1_ipd             = icer_sample$case_avert_age1_ipd,
                     case_avert_age2_ipd             = icer_sample$case_avert_age2_ipd,
                     case_avert_age3_ipd             = icer_sample$case_avert_age3_ipd,
                     case_avert_age1_opd             = icer_sample$case_avert_age1_opd,
                     case_avert_age2_opd             = icer_sample$case_avert_age2_opd,
                     case_avert_age3_opd             = icer_sample$case_avert_age3_opd,
                     pre_death_ipd                   = icer_sample$pre_death_ipd,
                     pre_death_opd                   = icer_sample$pre_death_opd,
                     pre_death_total                 = icer_sample$pre_death_total,
                     pre_death_age1                  = icer_sample$pre_death_age1,
                     pre_death_age2                  = icer_sample$pre_death_age2,
                     pre_death_age3                  = icer_sample$pre_death_age3,
                     post_death_ipd                  = icer_sample$post_death_ipd,
                     post_death_opd                  = icer_sample$post_death_opd,
                     post_death_total                = icer_sample$post_death_total,
                     post_death_age1                 = icer_sample$post_death_age1,
                     post_death_age2                 = icer_sample$post_death_age2,
                     post_death_age3                 = icer_sample$post_death_age3,
                     death_averted                   = icer_sample$death_averted,
                     death_averted_ipd               = icer_sample$death_averted_ipd,
                     death_averted_opd               = icer_sample$death_averted_opd,
                     death_averted_age1              = icer_sample$death_averted_age1,
                     death_averted_age2              = icer_sample$death_averted_age2,
                     death_averted_age3              = icer_sample$death_averted_age3,
                     dalys_pre_total                 = icer_sample$dalys_pre_total,
                     dalys_post_total                = icer_sample$dalys_post_total,
                     incremental_daly_total          = icer_sample$incremental_daly_total,
                     icer_daly_ipd                   = icer_sample$icer_daly_ipd,
                     icer_daly_opd                   = icer_sample$icer_daly_opd,
                     icer_daly_total                 = icer_sample$icer_daly_total 
                     
                     
  )]
}

# icer plot 
options(scipen=9999)
plot ( x = icer_dt$incremental_daly_total, y =icer_dt$incremental_cost_total)

# cea plane
cea_plane <- ggplot(data = icer_dt, aes(x=incremental_daly_total,
                                        y=incremental_cost_total))+
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("DALYs averted") +
  ylab("Incremental cost") +
  theme_bw()
cea_plane

# ceac
print(summary(icer_dt$icer_daly_total))
print(quantile(icer_dt$icer_daly_total, c(.5 , .025, .975 )))

# Pre-vaccination burden of disease (case by age-group)
prevacc_age1 <- quantile(icer_dt$prevacc_age1, c(0.025, 0.5, 0.975))
prevacc_age1 <- as.data.table(prevacc_age1)
prevacc_age1$ui_interval <- c(2.5, 50, 97.5)

prevacc_age2 <- quantile(icer_dt$prevacc_age2, c(0.025, 0.5, 0.975))
prevacc_age2 <- as.data.table(prevacc_age2)
prevacc_age2$ui_interval <- c(2.5, 50, 97.5)

prevacc_age3 <- quantile(icer_dt$prevacc_age3, c(0.025, 0.5, 0.975))
prevacc_age3 <- as.data.table(prevacc_age3)
prevacc_age3$ui_interval <- c(2.5, 50, 97.5)

# Pre-vaccination burden of disease (hospitalization by age-group)
prevacc_age1_ipd <- quantile(icer_dt$prevacc_age1_ipd, c(0.025, 0.5, 0.975))
prevacc_age1_ipd <- as.data.table(prevacc_age1_ipd)
prevacc_age1_ipd$ui_interval <- c(2.5, 50, 97.5)

prevacc_age2_ipd <- quantile(icer_dt$prevacc_age2_ipd, c(0.025, 0.5, 0.975))
prevacc_age2_ipd <- as.data.table(prevacc_age2_ipd)
prevacc_age2_ipd$ui_interval <- c(2.5, 50, 97.5)

prevacc_age2_ipd <- quantile(icer_dt$prevacc_age2_ipd, c(0.025, 0.5, 0.975))
prevacc_age2_ipd <- as.data.table(prevacc_age2_ipd)
prevacc_age2_ipd$ui_interval <- c(2.5, 50, 97.5)

prevacc_age3_ipd <- quantile(icer_dt$prevacc_age3_ipd, c(0.025, 0.5, 0.975))
prevacc_age3_ipd <- as.data.table(prevacc_age3_ipd)
prevacc_age3_ipd$ui_interval <- c(2.5, 50, 97.5)

# Pre-vaccination burden of disease (death by age-group)
pre_death_age1    <- quantile(icer_dt$pre_death_age1, c(0.025, 0.5, 0.975))
pre_death_age1    <- as.data.table(pre_death_age1)
pre_death_age1$ui_interval <- c(2.5, 50, 97.5)

pre_death_age2    <- quantile(icer_dt$pre_death_age2, c(0.025, 0.5, 0.975))
pre_death_age2    <- as.data.table(pre_death_age2)
pre_death_age2$ui_interval <- c(2.5, 50, 97.5)

pre_death_age3    <- quantile(icer_dt$pre_death_age3, c(0.025, 0.5, 0.975))
pre_death_age3    <- as.data.table(pre_death_age3)
pre_death_age3$ui_interval <- c(2.5, 50, 97.5)

# vaccine impact
# case avert by age gorup
case_avert_age1 <- quantile(icer_dt$case_avert_age1, c(0.025, 0.5, 0.975))
case_avert_age1 <- as.data.table(case_avert_age1)
case_avert_age1$ui_interval <- c(2.5, 50, 97.5)

case_avert_age2 <- quantile(icer_dt$case_avert_age2, c(0.025, 0.5, 0.975))
case_avert_age2 <- as.data.table(case_avert_age2)
case_avert_age2$ui_interval <- c(2.5, 50, 97.5)

case_avert_age3 <- quantile(icer_dt$case_avert_age3, c(0.025, 0.5, 0.975))
case_avert_age3 <- as.data.table(case_avert_age3)
case_avert_age3$ui_interval <- c(2.5, 50, 97.5)

#hospitalization avert by age group
ipd_avert_age1 <- quantile(icer_dt$case_avert_age1_ipd, c(0.025, 0.5, 0.975))
ipd_avert_age1 <- as.data.table(ipd_avert_age1)
ipd_avert_age1$ui_interval <- c(2.5, 50, 97.5)

ipd_avert_age2 <- quantile(icer_dt$case_avert_age2_ipd, c(0.025, 0.5, 0.975))
ipd_avert_age2 <- as.data.table(ipd_avert_age2)
ipd_avert_age2$ui_interval <- c(2.5, 50, 97.5)

ipd_avert_age3 <- quantile(icer_dt$case_avert_age3_ipd, c(0.025, 0.5, 0.975))
ipd_avert_age3 <- as.data.table(ipd_avert_age3)
ipd_avert_age3$ui_interval <- c(2.5, 50, 97.5)


# bar graph (Case averted by age group)
age_group <- c(rep("0-4yrs",2), rep("5-9yrs",2), rep("10-14yrs",2))
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

case <- c(prevacc_age1  = quantile(icer_dt$prevacc_age1, c(0.5)),
          postvacc_age1 = quantile(icer_dt$postvacc_age1, c(0.5)),
          prevacc_age2  = quantile(icer_dt$prevacc_age2, c(0.5)),
          postvacc_age2 = quantile(icer_dt$postvacc_age2, c(0.5)),
          prevacc_age3  = quantile(icer_dt$prevacc_age3, c(0.5)),
          postvacc_age3 = quantile(icer_dt$postvacc_age3, c(0.5)))

ui_low    <- c(ui_pre_age1_low      = quantile(icer_dt$prevacc_age1, c(0.025)),
               ui_post_age1_low     = quantile(icer_dt$postvacc_age1, c(0.025)),
               ui_pre_age2_low      = quantile(icer_dt$prevacc_age2, c(0.025)),
               ui_post_age2_low     = quantile(icer_dt$postvacc_age2, c(0.025)),
               ui_pre_age3_low      = quantile(icer_dt$prevacc_age3, c(0.025)),
               ui_post_age3_low     = quantile(icer_dt$postvacc_age3, c(0.025)))

ui_high    <- c(ui_pre_age1_high      = quantile(icer_dt$prevacc_age1, c(0.975)),
                ui_post_age1_high     = quantile(icer_dt$postvacc_age1, c(0.97)),
                ui_pre_age2_high      = quantile(icer_dt$prevacc_age2, c(0.975)),
                ui_post_age2_high     = quantile(icer_dt$postvacc_age2, c(0.975)),
                ui_pre_age3_high      = quantile(icer_dt$prevacc_age3, c(0.975)),
                ui_post_age3_high     = quantile(icer_dt$postvacc_age3, c(0.975)))

data <- data.frame(age_group, case, prepost, ui_low, ui_high)

data <- data.frame(age_group = c("0-4yrs", "0-4yrs","5-9yrs", "5-9yrs","10-14yrs", "10-14yrs"),
                   prepost   = factor(c("pre-vaccination","post-vaccination","pre-vaccination",
                                        "post-vaccination","pre-vaccination","post-vaccination"),
                                      levels = c("pre-vaccination", "post-vaccination")))
data$age_group <- factor(data$age_group, levels = c("0-4yrs", "5-9yrs", "10-14yrs"))

  ggplot(data = data, aes(fill = prepost, y =case, x = age_group))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = age_group, ymin = ui_low, ymax = ui_high),
                width = 0.2, position=position_dodge(.9)) +
    labs(x = "Age group", y = "case per 100,000 persons", fill = "vaccination status") +
    theme_bw()

# bar graph (deaths averted by age group)
age_group <- c(rep("0-4yrs",2), rep("5-9yrs",2), rep("10-14yrs",2))
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

death <-c(pre_death_age1  = quantile(icer_dt$pre_death_age1, c(0.5)), 
          post_death_age1 = quantile(icer_dt$post_death_age1, c(0.5)),
          pre_death_age2  = quantile(icer_dt$pre_death_age2, c(0.5)),
          post_death_age2 = quantile(icer_dt$post_death_age2, c(0.5)),
          pre_death_age3  = quantile(icer_dt$pre_death_age3, c(0.5)),
          post_death_age3 = quantile(icer_dt$post_death_age3, c(0.5)))

ui_low    <- c(ui_pre_age1_low      = quantile(icer_dt$pre_death_age1, c(0.025)),
               ui_post_age1_low     = quantile(icer_dt$post_death_age1, c(0.025)),
               ui_pre_age2_low      = quantile(icer_dt$pre_death_age2, c(0.025)),
               ui_post_age2_low     = quantile(icer_dt$post_death_age2, c(0.025)),
               ui_pre_age3_low      = quantile(icer_dt$pre_death_age3, c(0.025)),
               ui_post_age3_low     = quantile(icer_dt$post_death_age3, c(0.025)))

ui_high    <- c(ui_pre_age1_high      = quantile(icer_dt$pre_death_age1, c(0.975)),
                ui_post_age1_high     = quantile(icer_dt$post_death_age1, c(0.97)),
                ui_pre_age2_high      = quantile(icer_dt$pre_death_age2, c(0.975)),
                ui_post_age2_high     = quantile(icer_dt$post_death_age2, c(0.975)),
                ui_pre_age3_high      = quantile(icer_dt$pre_death_age3, c(0.975)),
                ui_post_age3_high     = quantile(icer_dt$post_death_age3, c(0.975)))


data <- data.frame(age_group, death, prepost, ui_low, ui_high)
data <- data.frame(age_group = c("0-4yrs", "0-4yrs","5-9yrs", "5-9yrs","10-14yrs", "10-14yrs"),
                   prepost   = factor(c("pre-vaccination","post-vaccination","pre-vaccination",
                                        "post-vaccination","pre-vaccination","post-vaccination"),
                                      levels = c("pre-vaccination", "post-vaccination")))
data$age_group <- factor(data$age_group, levels = c("0-4yrs", "5-9yrs", "10-14yrs"))

ggplot(data = data, aes(fill = prepost, y = death, x = age_group))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = age_group, ymin = ui_low, ymax = ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Age group", y = "Death per 100,000 persons", fill = "vaccination status") +
  theme_bw()

# bar graph (hospitalisation averted by age group)
age_group <- c(rep("0-4yrs",2), rep("5-9yrs",2), rep("10-14yrs",2))
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

ipd <-  c(prevacc_age1_ipd  = quantile(icer_dt$prevacc_age1_ipd, c(0.5)), 
          postvacc_age1_ipd = quantile(icer_dt$postvacc_age1_ipd, c(0.5)),
          prevacc_age2_ipd  = quantile(icer_dt$prevacc_age2_ipd, c(0.5)),
          postvacc_age2_ipd = quantile(icer_dt$postvacc_age2_ipd, c(0.5)),
          prevacc_age3_ipd  = quantile(icer_dt$prevacc_age3_ipd, c(0.5)),
          postvacc_age3_ipd = quantile(icer_dt$postvacc_age3_ipd, c(0.5)))

ui_low    <- c(ui_pre_age1_low      = quantile(icer_dt$prevacc_age1_ipd, c(0.025)),
               ui_post_age1_low     = quantile(icer_dt$postvacc_age1_ipd, c(0.025)),
               ui_pre_age2_low      = quantile(icer_dt$prevacc_age2_ipd, c(0.025)),
               ui_post_age2_low     = quantile(icer_dt$postvacc_age2_ipd, c(0.025)),
               ui_pre_age3_low      = quantile(icer_dt$prevacc_age3_ipd, c(0.025)),
               ui_post_age3_low     = quantile(icer_dt$postvacc_age3_ipd, c(0.025)))

ui_high    <- c(ui_pre_age1_high      = quantile(icer_dt$prevacc_age1_ipd, c(0.975)),
                ui_post_age1_high     = quantile(icer_dt$postvacc_age1_ipd, c(0.97)),
                ui_pre_age2_high      = quantile(icer_dt$prevacc_age2_ipd, c(0.975)),
                ui_post_age2_high     = quantile(icer_dt$postvacc_age2_ipd, c(0.975)),
                ui_pre_age3_high      = quantile(icer_dt$prevacc_age3_ipd, c(0.975)),
                ui_post_age3_high     = quantile(icer_dt$postvacc_age3_ipd, c(0.975)))

data <- data.frame(age_group, ipd, prepost, ui_low, ui_high)
data <- data.frame(age_group = c("0-4yrs", "0-4yrs","5-9yrs", "5-9yrs","10-14yrs", "10-14yrs"),
                   prepost   = factor(c("pre-vaccination","post-vaccination","pre-vaccination",
                                        "post-vaccination","pre-vaccination","post-vaccination"),
                                      levels = c("pre-vaccination", "post-vaccination")))
data$age_group <- factor(data$age_group, levels = c("0-4yrs", "5-9yrs", "10-14yrs"))

ggplot(data = data, aes(fill = prepost, y = ipd, x = age_group))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = age_group, ymin = ui_low, ymax = ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Age group", y = "Hospitalisation per 100,000 persons", fill = "vaccination status") +
  theme_bw()

# probability 
probability_cea <- seq(from=0, to= 1, by= .01)
wtp<- quantile(icer_dt$icer_daly_total, probability_cea)
wtp_dt <- as.data.table(wtp)

# wtp table with probs
wtp_prob <- cbind (wtp_dt, probability_cea)

# keep only positive wtp values (negative values are dominant values, no need to show)
wtp_post <- wtp_prob %>% filter(wtp > 0)

# ceac ggplot
ceac <- ggplot(data = wtp_post, aes(x = wtp, y = probability_cea)) +
  geom_line(color = 'black') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Willingenss to pay") +
  ylab("Probability of cost-effectiveness") +
  theme_bw()
ceac

# opportunity costs lines 

vertDf <- data.frame(wtp = c(166, 279, 2191), labels = c("lower-bound", "upper-bound", "GDP per capita 2021"))

# add two graphs 

ceac + geom_vline(aes(xintercept = wtp, color = labels), data = vertDf, show.legend=T) +
  scale_colour_manual("WTP", values = c("lower-bound" = "blue", "upper-bound" = "red", "GDP per capita 2021" = "green")) +
  theme_bw() +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"))

# ------------------------------------------------------------------------------
# function - icer calculation
# icer calculation - add three columns 
# (i)   icer_mid for mid value of each parameter (should be same value for all parameters)
# (ii)  icer_low for low value of each parameter while other parameters have mid values
# (iii) icer_high for high value of each parameter while other parameters have mid values
# ------------------------------------------------------------------------------
icer_calc <- function (param_dt) {
  
  # add 3 empty columns for icer - mid, low, high (columns 5, 6, 7)
  param_dt [, icer_mid  := 0]
  param_dt [, icer_low  := 0]
  param_dt [, icer_high := 0]
  
  # loop through parameters
  for (i in 1:nrow (param_dt)) {
    
    # loop through mid, low, and high values (which are in columns 2, 3, 4)
    for (j in 2:4) {
      
      # select the mid value for each parameter
      v_cost               <- param_dt [parameter == "v_cost",               mid ]
      delivery_cost        <- param_dt [parameter == "delivery_cost",        mid ]
      tot_pop              <- param_dt [parameter == "tot_pop",              mid ]
      vacc_pop             <- param_dt [parameter == "vacc_pop",             mid ]
      postvacc_p1_ipd      <- param_dt [parameter == "postvacc_p1_ipd",      mid ]
      postvacc_p1_opd      <- param_dt [parameter == "postvacc_p1_opd",      mid ]
      postvacc_p1          <- param_dt [parameter == "postvacc_p1",          mid ]
      facility_cost        <- param_dt [parameter == "facility_cost",        mid ]
      dmc_ipd              <- param_dt [parameter == "dmc_ipd",              mid ]
      dmc_opd              <- param_dt [parameter == "dmc_opd",              mid ]
      dmc                  <- param_dt [parameter == "dmc",                  mid ]
      dnmc_ipd             <- param_dt [parameter == "dnmc_ipd",             mid ]
      dnmc_opd             <- param_dt [parameter == "dnmc_opd",             mid ]
      dnmc                 <- param_dt [parameter == "dnmc",                 mid ]
      indirect_ipd         <- param_dt [parameter == "indirect_ipd",         mid ]
      indirect_opd         <- param_dt [parameter == "indirect_opd",         mid ]
      indirect             <- param_dt [parameter == "indirect",             mid ]
      ve                   <- param_dt [parameter == "ve",                   mid ]
      dw_ipd               <- param_dt [parameter == "dw_ipd",               mid ]
      dw_opd               <- param_dt [parameter == "dw_opd",               mid ]
      dw                   <- param_dt [parameter == "dw",                   mid ]
      illness_duration_ipd <- param_dt [parameter == "illness_duration_ipd", mid ]
      illness_duration_opd <- param_dt [parameter == "illness_duration_opd", mid ]
      cfr_ipd              <- param_dt [parameter == "cfr_ipd",              mid ]
      cfr_opd              <- param_dt [parameter == "cfr_opd",              mid ]
      age_death            <- param_dt [parameter == "age_death",            mid ]
      life_exp             <- param_dt [parameter == "life_exp",             mid ]
      
      
      
      
      # ------------------------------------------------------------------------
      # assign value for a specific parameter to mid/low/high values
      assign (x     = param_dt [i, parameter], 
              value = param_dt [i, j, with = FALSE] )
      # ------------------------------------------------------------------------
      
      # ------------------------------------------------------------------------
      # formulas for icer
      
      totcost_vacc_ipd       <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_ipd) * (facility_cost + dmc_ipd + dnmc_ipd + indirect_ipd)
      
      totcost_vacc_opd       <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_opd) * (facility_cost + dmc_opd + dnmc_opd + indirect_opd)
      
      totcost_vacc_tot       <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_ipd) * (facility_cost + dmc_ipd + dnmc_ipd + indirect_ipd) + (postvacc_p1_opd) * (facility_cost + dmc_opd + dnmc_opd + indirect_opd)
      
      prevacc_ipd            <- (postvacc_p1_ipd)/(1-ve*(vacc_pop/tot_pop))
      
      prevacc_opd            <- (postvacc_p1_opd)/(1-ve*(vacc_pop/tot_pop))
      
      prevacc                <- (prevacc_ipd)  + (prevacc_opd)
      
      totcost_unvacc         <- (prevacc)*(facility_cost + dmc + dnmc + indirect)
      
      incremental_cost_ipd   <- (totcost_vacc_ipd) - (totcost_unvacc)
      
      incremental_cost_opd   <- (totcost_vacc_opd) - (totcost_unvacc)
      
      incremental_cost_total <- (totcost_vacc_tot) - (totcost_unvacc)
      
      case_averted_ipd       <- (prevacc_ipd)      - (postvacc_p1_ipd)
      
      case_averted_opd       <- (prevacc_opd)      - (postvacc_p1_opd)
      
      case_averted_total     <- (case_averted_ipd) + (case_averted_opd)
      
      icer_ipd               <- incremental_cost_ipd / case_averted_ipd
      
      icer_opd               <- incremental_cost_opd / case_averted_opd
      
      icer_tot               <- incremental_cost_total / case_averted_total
      
      distance               <- (life_exp) - (age_death)
      
      pre_death_ipd          <- (cfr_ipd)*(prevacc_ipd)
      
      pre_death_opd          <- (cfr_opd)*(prevacc_opd)
      
      pre_death_total        <- (pre_death_ipd) + (pre_death_opd)
      
      post_death_ipd         <- (cfr_ipd)*(postvacc_p1_ipd)
      
      post_death_opd         <- (cfr_opd)*(postvacc_p1_opd)
      
      post_death_total       <-  (post_death_ipd) + (post_death_opd) 
      
      death_averted          <-  (pre_death_total) - (post_death_total)
      
      death_averted_ipd      <-  (pre_death_ipd) - (post_death_ipd)
      
      death_averted_opd      <-  (pre_death_opd) - (post_death_opd)
      
      yll_pre_ipd            <- (pre_death_ipd)*(distance)
      
      yll_pre_opd            <- (pre_death_opd)*(distance)
      
      yll_pre_total          <- (yll_pre_ipd)  + (yll_pre_opd)
      
      yll_post_ipd           <- (post_death_ipd)*(distance)
      
      yll_post_opd           <- (post_death_opd)*(distance)
      
      yll_post_total         <- (yll_post_ipd) + (yll_post_opd)
      
      yld_pre_ipd            <- (prevacc_ipd)*(dw_ipd)*(illness_duration_ipd)
      
      yld_pre_opd            <- (prevacc_opd)*(dw_opd)*(illness_duration_opd)
      
      yld_pre_total          <- (yld_pre_ipd) + (yld_pre_opd)
      
      dalys_pre_total        <- (yll_pre_total) + (yld_pre_total)
      
      yld_post_ipd           <- (postvacc_p1_ipd)*(dw_ipd)*(illness_duration_ipd)
      
      yld_post_opd           <- (postvacc_p1_opd)*(dw_opd)*(illness_duration_opd)
      
      yld_post_total         <- (yld_post_ipd) + (yld_post_opd)
      
      dalys_post_total       <- (yll_post_total) + (yld_post_total)
      
      yld_averted_ipd        <- (prevacc_ipd - postvacc_p1_ipd)* (dw_ipd) *(illness_duration_ipd)
      
      yld_averted_opd        <- (prevacc_opd - postvacc_p1_opd)* (dw_opd) *(illness_duration_opd)
      
      yld_averted_total      <- (yld_averted_ipd) + (yld_averted_opd) 
      
      yll_averted_ipd        <- (yll_pre_ipd)  -(yll_post_ipd)
      
      yll_averted_opd        <- (yll_pre_opd)  -(yll_post_opd)
      
      yll_averted_total      <- (yll_averted_ipd) + (yll_averted_opd)
      
      incremental_daly_ipd   <- (yld_averted_ipd) + (yll_averted_ipd)
      
      incremental_daly_opd   <- (yld_averted_opd) + (yll_averted_opd)
      
      incremental_daly_total <- (yll_averted_total) + (yld_averted_total)
      
      icer_daly_ipd          <- incremental_cost_ipd / incremental_daly_ipd
      
      icer_daly_opd          <- incremental_cost_opd / incremental_daly_opd
      
      icer_daly_total        <- incremental_cost_total / incremental_daly_total
      
      
      # ----------------------------------------------------------------------
      
      # ------------------------------------------------------------------------
      # save icer_daly value to corresponding row and column
      # note: icer (mid, low, high) columns are 3 columns after parameter (mid, low, high) columns
      param_dt [i, j + 3] <- icer_daly_total
      
      # ----------------------------------------------------------------------
      
    } # end of -- loop through mid, low, and high values
    
    
  } # end of -- loop through parameters
  
  # return data table with 3 additional columns for icer - mid, low, high
  return (param_dt)
  
} # end of function - icer calculation
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# main program (start)
# ------------------------------------------------------------------------------


tornado_param <- sapply (lhs_sample, function(x) quantile(x, probs = seq(0, 1, 0.025)))
tornado_param <- as.data.table(tornado_param)  
tornado_param$ui_interval <- seq(0, 100, 2.5)
tornado_param_low  <- tornado_param %>% filter(ui_interval == 2.5)
tornado_param_mid  <- tornado_param %>% filter(ui_interval == 50)
tornado_param_high <- tornado_param %>% filter(ui_interval == 97.5)
tornado_param_all  <- rbind (tornado_param_low, tornado_param_mid, tornado_param_high)

# initialise parameter data table
param <- data.frame (v_cost               = c (mid = 2.96,   low = 2.96,   high = 2.96), 
                     delivery_cost        = c (mid = 1.489151,   low = 1.372637,  high = 1.611899),
                     tot_pop              = c (mid = 159831, low = 159831, high = 159831),
                     vacc_pop             = c (mid = 113420, low = 113420, high = 113420),
                     postvacc_p1_ipd      = c (mid = 200,     low = 200,     high = 200),
                     postvacc_p1_opd      = c (mid = 459,     low = 459,     high = 459),
                     postvacc_p1          = c (mid = 659,     low = 659,     high = 659),
                     facility_cost        = c (mid = 93.72781,  low = 72.05578, high = 119.42376 ),
                     dmc_ipd              = c (mid = 145.085922, low = 2.443228, high = 963.940431),
                     dmc_opd              = c (mid = 60.1547144, low = 0.3206837, high = 535.8029029),
                     dmc                  = c (mid = 100.4446615, low = 0.7537101, high = 822.2790403 ),
                     dnmc_ipd             = c (mid = 39.421173,  low = 4.440873,  high = 144.388621), 
                     dnmc_opd             = c (mid = 12.3908994,  low = 0.1278219,  high = 93.308175), 
                     dnmc                 = c (mid = 26.416384,  low = 1.219572,  high = 130.697268), 
                     indirect_ipd         = c (mid = 73.83249,  low = 31.05202, high = 144.91829),
                     indirect_opd         = c (mid = 71.07460,  low = 31.30255, high = 135.47818),
                     indirect             = c (mid = 72.67934,  low = 31.25501, high = 140.55149),
                     ve                   = c (mid = 0.8285294,   low = 0.4651203,   high = 0.9841985), 
                     dw_ipd               = c (mid = 0.2099958,  low =0.1907379,  high = 0.2298767),
                     dw_opd               = c (mid = 0.05139319,  low =0.03475859,  high = 0.07320066),
                     illness_duration_ipd = c (mid = 1.362,  low = 1.30890, high = 4.3047),
                     illness_duration_opd = c (mid = 0.04154,  low = 0.010775, high = 0.106362),
                     cfr_ipd              = c (mid = 0.02780828,  low = 0.02057399,  high = 0.03651082),
                     cfr_opd              = c (mid = 0.01863717,  low = 0.01107713,  high = 0.02896699),
                     age_death            = c (mid = 7.547850,    low = 6.562944,  high = 8.624685),
                     life_exp             = c (mid = 72.67957,  low = 69.27136, high = 76.08917) )


# transpose the table
param_table <- data.frame (t(param))

# change to data table
param_table <- setDT (param_table)

# give name to each column
names (param_table)   <- rownames (param)

# add parameter column
param_table$parameter <- colnames (param)

# set parameter as first column (column reordering)
setcolorder (x        = param_table, 
             neworder = "parameter")

# icer calculation - add three columns 
# (i)   icer_mid for mid value of each parameter (should be same value for all parameters)
# (ii)  icer_low for low value of each parameter while other parameters have mid values
# (iii) icer_high for high value of each parameter while other parameters have mid values
param_table <- icer_calc (param_dt = param_table)

# Note: 
# When a parameter takes a low/high value in comparison to mid value, icer value 
# can decrease/increase depending on the relationship in the same/inverse direction.

# ------------------------------------------------------------------------------
# basic tornado plot 


param_table <- param_table[-c(1,3,4,5,6,7)]

# length of the ICER diff (high - low)  

param_table <- param_table %>%
  mutate(length = (icer_high) - (icer_low)) 

param_table <- arrange(param_table, length)

# absolute value of length

param_table$length_abs <- abs(param_table$length)

# arrange the table again by the order of absolute value 

param_table <- param_table %>% 
  arrange(desc(length_abs))
# tornado plot 

tornado_plot <- ggplot (param_table,
                        aes (reorder(parameter, length_abs), 
                             ymin = icer_low, 
                             ymax = icer_high, 
                             color = parameter)) +
  geom_linerange (size = 5) +
  geom_hline (yintercept = param_table$icer_mid [1], 
              linetype   = "dotted") +
  xlab ("parameter") +
  ylab ("ICER (cost per DALY averted)") + 
  coord_flip() +
  theme_bw () 

print (tornado_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
toc ()
# main program (end)
# ------------------------------------------------------------------------------