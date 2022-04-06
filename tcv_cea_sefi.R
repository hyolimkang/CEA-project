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
library (ggpubr)
library (gridExtra)
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
                k = 24) 
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
lhs_sample [,1]  <- 3.274093
lhs_sample [,2]  <- qgamma (p = A[,2], shape = ((1.695151/0.579475561)^2), rate = (1.695151/(0.579475561)^2),
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,3]  <- 159831
lhs_sample [,4]  <- 113420
lhs_sample [,5]  <- qlnorm (p = A[,5], meanlog = 6.244167, sdlog =  0.1587704, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,6]  <- qgamma (p = A[,6], shape = ((97.33/12.33)^2), rate = (97.33/(12.33)^2), 
                             lower.tail = TRUE, log.p = FALSE)
lhs_sample [,7]  <- qgamma (p = A[,7], shape = ((234.7688/265.9934)^2), rate = (234.7688/(265.9934)^2), 
                             lower.tail = TRUE, log.p = FALSE)
lhs_sample [,8] <- qgamma (p = A[,8], shape = ((115.3841/149.9214)^2), rate = (115.3841/(149.9214)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,9] <- qgamma (p = A[,9], shape = ((183.3415/229.8478)^2), rate = (183.3415/(229.8478)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,10] <- qgamma (p = A[,10], shape = ((48.58656/37.39749)^2), rate = (48.58656/(37.39749)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,11] <- qgamma (p = A[,11],  shape = ((21.48999/25.8971)^2), rate = (21.48999/(25.8971)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,12] <- qgamma (p = A[,12],  shape = ((36.91419/35.39041)^2), rate = (36.91419/(35.39041)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,13] <- qgamma (p = A[,13],  shape = ((77.53769/29.51196)^2), rate = (77.53769/(29.51196)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,14] <- qgamma (p = A[,14],   shape = ((74.28838/26.90198)^2), rate = (74.28838/(26.90198)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,15] <- qgamma (p = A[,15],  shape = ((76.13798/28.24641)^2),rate = (76.13798/(28.24641)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,16] <- qbeta  (p = A[,16], shape1 = 5.362628, shape2 = 1.20922, ncp=0,lower.tail = TRUE, log.p = FALSE)
lhs_sample [,17] <- qbeta (p = A[,17], shape1 = 348.18, shape2 = 1309.82, ncp=0, lower.tail = TRUE, log.p = FALSE) 
lhs_sample [,18] <- qbeta (p = A[,18], shape1 = 25.58192, shape2 = 466.3781, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,19] <- qgamma (p = A[,19], shape = ((0.0404/0.01532)^2), rate = (0.0404/(0.01532)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,20] <- qgamma (p = A[,20], shape = ((0.04559688/ 0.02500143)^2), rate = (0.04559688/(0.02500143)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,21] <- qbeta (p = A[,21], shape1 = 9.531057, shape2 = 330.8638, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,22] <- qbeta (p = A[,22], shape1 =  21.23796, shape2 = 1096.55, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,23] <- qgamma (p = A[,23], shape = ((7.56/0.527)^2), rate = (7.56/(0.527)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample [,24] <- qnorm (p = A[,24], mean = 72.68, sd = 1.73979, lower.tail = TRUE, log.p = FALSE)

# change the name of the lhs_sample
cols <- c(  "v_cost", 
            "delivery_cost",    
            "tot_pop",
            "vacc_pop",  
            "prevacc",      
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

compute_icer <- function (v_cost            = 3.274093,
                          delivery_cost,
                          tot_pop           = 159831,
                          vacc_pop          = 113420,
                          prevacc, 
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
  
  # prevacc (IPD)
  
  prevacc_ipd <- prevacc * 0.011
  
  # prevacc (OPD)
  
  prevacc_opd <- prevacc - prevacc_ipd
  
  # post-vacc
  
  postvacc <-  (prevacc)*(1 - ve*(vacc_pop/tot_pop))
  
  # post-vacc (IPD)
  
  postvacc_ipd <- (prevacc_ipd)*(1 - ve*(vacc_pop/tot_pop))
  
  # post-vacc (OPD)
  
  postvacc_opd <- postvacc - postvacc_ipd
  
  # case-averted (per 100 000)
  
  case_avert <- (prevacc - postvacc)
  
  # case-averted IPD (per 100 000)
  
  case_avert_ipd <- (case_avert) * 0.011
  
  # case-averted OPD (per 100 000)
  
  case_avert_opd <- (case_avert) - (case_avert_ipd)
  
  # incremental cost of inpatients 
  
  incremental_cost_ipd   <- (v_cost + delivery_cost) * (vacc_pop) - (case_avert_ipd)*(facility_cost + dmc_ipd + dnmc_ipd)
  
  # incremental cost of outpatients
  
  incremental_cost_opd   <- (v_cost + delivery_cost) * (vacc_pop) - (case_avert_opd)*(facility_cost + dmc_opd + dnmc_opd)
  
  # incremental cost total (used for the final value as a delta C)
  
  incremental_cost_total <-  (v_cost + delivery_cost)* (vacc_pop) - (case_avert)*(facility_cost + dmc +dnmc + indirect) 
  
  # YLL calculation component: amount of years of life loss between life-expectancy and average age at death
  
  distance               <- (life_exp) - (age_death)
  
  # number of deaths in no-vaccination situation (inpatients)
  
  pre_death_ipd          <- (cfr_ipd)*(prevacc_ipd)
  
  # number of deaths in no-vaccination situation (outpatients)
  
  pre_death_opd          <- (cfr_opd)*(prevacc_opd)
  
  # number of deaths in no-vaccination situation (total)
  
  pre_death_total        <- (pre_death_ipd) + (pre_death_opd)
  
  # number of deaths post-vaccination situation (inpatients)
  
  post_death_ipd         <- (cfr_ipd)*(postvacc_ipd)

  # number of deaths post-vaccination situation (outpatients)
  
  post_death_opd         <- (cfr_opd)*(postvacc_opd)
  
  # number of deaths post-vaccination situation (total)
  
  post_death_total       <-  (post_death_ipd) + (post_death_opd) 
  
  # number of deaths averted (total)
  
  death_averted          <-  (pre_death_total) - (post_death_total)
  
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
  
  yld_pre_ipd            <- (prevacc_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  # years of life with disease for outpatients in pre-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
  
  yld_pre_opd            <- (prevacc_opd)*(dw_opd)*(illness_duration_opd)
  
  # total YLD of pre-vaccination typhoid cases
  
  yld_pre_total          <- (yld_pre_ipd) + (yld_pre_opd)
  
  # total DALYs for pre-vaccination typhoid cases (YLL total + YLD total)
  
  dalys_pre_total        <- (yll_pre_total) + (yld_pre_total)
  
  # years of life with disease for inpatients in post-vacc situation (number of inpatient cases * disability weight * total illness duration presented as years)
  
  yld_post_ipd           <- (postvacc_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  # years of life with disease for outpatients in post-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
  
  yld_post_opd           <- (postvacc_opd)*(dw_opd)*(illness_duration_opd)
  
  # total YLD for post-vaccination typhoid cases (YLL total + YLD total)
  
  yld_post_total         <- (yld_post_ipd) + (yld_post_opd)
  
  # total DALYs for post-vaccination typhoid cases (YLL total + YLD total)
  
  dalys_post_total       <- (yll_post_total) + (yld_post_total)
  
  # YLD averted for inpatients
  
  yld_averted_ipd        <- (case_avert_ipd)* (dw_ipd) *(illness_duration_ipd)
  
  # YLD averted for outpatients
  
  yld_averted_opd        <- (case_avert_opd)* (dw_opd) *(illness_duration_opd)
  
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
  
  
  
  return (list (prevacc                          = prevacc,
                prevacc_ipd                      = prevacc_ipd,
                prevacc_opd                      = prevacc_opd,
                postvacc                         = postvacc,
                postvacc_ipd                     = postvacc_ipd,
                postvacc_opd                     = postvacc_opd,
                case_avert                       = case_avert,
                case_avert_ipd                   = case_avert_ipd,
                case_avert_opd                   = case_avert_opd,
                incremental_cost_ipd             = incremental_cost_ipd,
                incremental_cost_opd             = incremental_cost_opd,
                incremental_cost_total           = incremental_cost_total,
                pre_death_ipd                    = pre_death_ipd,
                pre_death_opd                    = pre_death_opd,
                pre_death_total                  = pre_death_total,
                post_death_ipd                   = post_death_ipd,
                post_death_opd                   = post_death_opd,
                post_death_total                 = post_death_total,
                death_averted                    = death_averted,
                death_averted_ipd                = death_averted_ipd,
                death_averted_opd                = death_averted_opd,
                yll_pre_ipd                      = yll_pre_ipd,
                yll_pre_opd                      = yll_pre_opd,
                yll_pre_total                    = yll_pre_total,
                yll_post_ipd                     = yll_post_ipd,
                yll_post_opd                     = yll_post_opd,
                yll_post_total                   = yll_post_total,
                yld_pre_ipd                      = yld_pre_ipd,
                yld_pre_opd                      = yld_pre_opd,
                yld_pre_total                    = yld_pre_total,
                dalys_pre_total                  = dalys_pre_total,
                yld_post_ipd                     = yld_post_ipd,
                yld_post_opd                     = yld_post_opd,
                yld_post_total                   = yld_post_total,
                dalys_post_total                 = dalys_post_total,
                yld_averted_ipd                  = yld_averted_ipd,
                yld_averted_opd                  = yld_averted_opd,
                yld_averted_total                = yld_averted_total,
                yll_averted_ipd                  = yll_averted_ipd,
                yll_averted_opd                  = yll_averted_opd,
                yll_averted_total                = yll_averted_total,
                incremental_daly_ipd             = incremental_daly_ipd,
                incremental_daly_opd             = incremental_daly_opd,
                incremental_daly_total           = incremental_daly_total,
                icer_daly_ipd                    = icer_daly_ipd,
                icer_daly_opd                    = icer_daly_opd,
                icer_daly_total                  = icer_daly_total)) 
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
                               lhs_sample$tot_pop                   [i],
                               lhs_sample$vacc_pop                  [i],
                               lhs_sample$prevacc                   [i],
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
  
  icer_dt [i, `:=` ( prevacc                         = icer_sample$prevacc,
                     prevacc_ipd                     = icer_sample$prevacc_ipd,
                     prevacc_opd                     = icer_sample$prevacc_opd,
                     postvacc                        = icer_sample$postvacc,
                     postvacc_ipd                    = icer_sample$postvacc_ipd,
                     postvacc_opd                    = icer_sample$postvacc_opd,
                     case_avert                      = icer_sample$case_avert,
                     case_avert_ipd                  = icer_sample$case_avert_ipd,
                     case_avert_opd                  = icer_sample$case_avert_opd,
                     incremental_cost_ipd            = icer_sample$incremental_cost_ipd,
                     incremental_cost_opd            = icer_sample$incremental_cost_opd,
                     incremental_cost_total          = icer_sample$incremental_cost_total,
                     pre_death_ipd                   = icer_sample$pre_death_ipd,
                     pre_death_opd                   = icer_sample$pre_death_opd,
                     pre_death_total                 = icer_sample$pre_death_total,
                     post_death_ipd                  = icer_sample$post_death_ipd,
                     post_death_opd                  = icer_sample$post_death_opd,
                     post_death_total                = icer_sample$post_death_total,
                     death_averted                   = icer_sample$death_averted,
                     death_averted_ipd               = icer_sample$death_averted_ipd,
                     death_averted_opd               = icer_sample$death_averted_opd,
                     yll_pre_ipd                     = icer_sample$yll_pre_ipd,
                     yll_pre_opd                     = icer_sample$yll_pre_opd,
                     yll_pre_total                   = icer_sample$yll_pre_total,
                     yll_post_ipd                    = icer_sample$yll_post_ipd,
                     yll_post_opd                    = icer_sample$yll_post_opd,
                     yll_post_total                  = icer_sample$yll_post_total,
                     yld_pre_ipd                     = icer_sample$yld_pre_ipd,
                     yld_pre_opd                     = icer_sample$yld_pre_opd,
                     yld_pre_total                   = icer_sample$yld_pre_total,
                     dalys_pre_total                 = icer_sample$dalys_pre_total,
                     yld_post_ipd                    = icer_sample$yld_post_ipd,
                     yld_post_opd                    = icer_sample$yld_post_opd,
                     yld_post_total                  = icer_sample$yld_post_total,
                     dalys_post_total                = icer_sample$dalys_post_total,
                     yld_averted_ipd                 = icer_sample$yld_averted_ipd,
                     yld_averted_opd                 = icer_sample$yld_averted_opd,
                     yld_averted_total               = icer_sample$yld_averted_total,
                     yll_averted_ipd                 = icer_sample$yll_averted_ipd,
                     yll_averted_opd                 = icer_sample$yll_averted_opd,
                     yll_averted_total               = icer_sample$yll_averted_total,
                     incremental_daly_ipd            = icer_sample$incremental_daly_ipd,
                     incremental_daly_opd            = icer_sample$incremental_daly_opd,
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

# Pre-vaccination burden of disease
prevacc <- quantile(icer_dt$prevacc, c(0.025, 0.5, 0.975))
prevacc <- as.data.table(prevacc)
prevacc$ui_interval <- c(2.5, 50, 97.5)

# Pre-vaccination burden of disease (hospitalization)
prevacc_ipd <- quantile(icer_dt$prevacc_ipd, c(0.025, 0.5, 0.975))
prevacc_ipd <- as.data.table(prevacc_ipd)
prevacc_ipd$ui_interval <- c(2.5, 50, 97.5)

# Pre-vaccination burden of disease (death by age-group)
pre_death    <- quantile(icer_dt$pre_death, c(0.025, 0.5, 0.975))
pre_death    <- as.data.table(pre_death)
pre_death$ui_interval <- c(2.5, 50, 97.5)

# vaccine impact
# case avert 
case_avert <- quantile(icer_dt$case_avert, c(0.025, 0.5, 0.975))
case_avert <- as.data.table(case_avert)
case_avert$ui_interval <- c(2.5, 50, 97.5)

#hospitalization avert 
ipd_avert <- quantile(icer_dt$case_avert_ipd, c(0.025, 0.5, 0.975))
ipd_avert <- as.data.table(ipd_avert)
ipd_avert$ui_interval <- c(2.5, 50, 97.5)

# Post-vaccination burden of disease
postvacc <- quantile(icer_dt$postvacc, c(0.025, 0.5, 0.975))
postvacc <- as.data.table(postvacc)
postvacc$ui_interval <- c(2.5, 50, 97.5)

# Post-vaccination burden of disease (hospitalization)
postvacc_ipd <- quantile(icer_dt$postvacc_ipd, c(0.025, 0.5, 0.975))
postvacc_ipd <- as.data.table(postvacc_ipd)
postvacc_ipd$ui_interval <- c(2.5, 50, 97.5)

# bar graph (Case averted)
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

case <- c(prevacc  = quantile(icer_dt$prevacc, c(0.5)),
          postvacc = quantile(icer_dt$postvacc, c(0.5)))

case_ui_low    <- c(ui_pre_low      = quantile(icer_dt$prevacc, c(0.025)),
               ui_post_low     = quantile(icer_dt$postvacc, c(0.025)))

case_ui_high    <- c(ui_pre_high      = quantile(icer_dt$prevacc, c(0.975)),
                ui_post_high     = quantile(icer_dt$postvacc, c(0.97)))

data_1 <- data.frame(case, prepost, case_ui_low, case_ui_high)

data_1 <- data.frame(prepost           = factor(c("pre-vaccination","post-vaccination"),
                                      levels = c("pre-vaccination", "post-vaccination")))
case_averted <- ggplot(data = data_1, aes(y = case, x = prepost))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = prepost, ymin = case_ui_low, ymax = case_ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Vaccination status", y = "case per 100,000 persons") +
  theme_bw()
case_averted

# bar graph (deaths averted)
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

death      <- c(predeath        = quantile(icer_dt$pre_death_total, c(0.5)),
               postdeath         = quantile(icer_dt$post_death_total, c(0.5)))

death_ui_low    <- c(ui_pre_low        = quantile(icer_dt$pre_death_total, c(0.025)),
               ui_post_low       = quantile(icer_dt$post_death_total, c(0.025)))

death_ui_high   <- c(ui_pre_high       = quantile(icer_dt$pre_death_total, c(0.975)),
               ui_post_high      = quantile(icer_dt$post_death_total, c(0.97)))

data_2 <- data.frame(death, prepost, death_ui_low, death_ui_high)

data_2 <- data.frame(prepost       = factor(c("pre-vaccination","post-vaccination"),
                                            levels = c("pre-vaccination", "post-vaccination")))
death_averted <- ggplot(data = data_2, aes(y = death, x = prepost))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = prepost, ymin = death_ui_low, ymax = death_ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Vaccination status", y = "Deaths per 100,000 persons") +
  theme_bw()

death_averted

# bar graph (hospitalisation averted)
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

ipd      <- c(pre_ipd           = quantile(icer_dt$prevacc_ipd, c(0.5)),
               post_ipd         = quantile(icer_dt$postvacc_ipd, c(0.5)))

ipd_ui_low    <- c(ui_pre_low        = quantile(icer_dt$prevacc_ipd, c(0.025)),
               ui_post_low       = quantile(icer_dt$postvacc_ipd, c(0.025)))

ipd_ui_high   <- c(ui_pre_high       = quantile(icer_dt$prevacc_ipd, c(0.975)),
               ui_post_high      = quantile(icer_dt$postvacc_ipd, c(0.97)))

data_3 <- data.frame(ipd, prepost, ipd_ui_low, ipd_ui_high)

data_3 <- data.frame(prepost       = factor(c("pre-vaccination","post-vaccination"),
                                          levels = c("pre-vaccination", "post-vaccination")))
ipd_averted <- ggplot(data = data_3, aes(y = ipd, x = prepost))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = prepost, ymin = ipd_ui_low, ymax = ipd_ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Vaccination status", y = "Hospitalisation per 100,000 persons") +
  theme_bw()
ipd_averted

# bar graph (DALYs averted)
prepost   <- rep(c("pre-vaccination", "post-vaccination"), 3)

daly      <- c(pre_daly          = quantile(icer_dt$dalys_pre_total, c(0.5)),
               post_daly         = quantile(icer_dt$dalys_post_total, c(0.5)))

daly_ui_low    <- c(ui_pre_low        = quantile(icer_dt$dalys_pre_total, c(0.025)),
               ui_post_low       = quantile(icer_dt$dalys_post_total, c(0.025)))

daly_ui_high   <- c(ui_pre_high       = quantile(icer_dt$dalys_pre_total, c(0.975)),
               ui_post_high      = quantile(icer_dt$dalys_post_total, c(0.97)))

data_4 <- data.frame(case, prepost, daly_ui_low, daly_ui_high)

data_4 <- data.frame(prepost       = factor(c("pre-vaccination","post-vaccination"),
                                          levels = c("pre-vaccination", "post-vaccination")))
daly_averted <- ggplot(data = data_4, aes(y = daly, x = prepost))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(x = prepost, ymin = daly_ui_low, ymax = daly_ui_high),
                width = 0.2, position=position_dodge(.9)) +
  labs(x = "Vaccination status", y = "DALYs per 100,000 persons") +
  theme_bw()
daly_averted

# combine multiple graphs on one page
# http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

figure <- ggarrange(case_averted, death_averted, ipd_averted, daly_averted,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
figure

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
      prevacc              <- param_dt [parameter == "prevacc",              mid ]
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
      # prevacc (IPD)
      
      # prevacc (IPD)
      
      prevacc_ipd <- prevacc * 0.011
      
      # prevacc (OPD)
      
      prevacc_opd <- prevacc - prevacc_ipd
      
      # post-vacc
      
      postvacc <-  (prevacc)*(1 - ve*(vacc_pop/tot_pop))
      
      # post-vacc (IPD)
      
      postvacc_ipd <- (prevacc_ipd)*(1 - ve*(vacc_pop/tot_pop))
      
      # post-vacc (OPD)
      
      postvacc_opd <- postvacc - postvacc_ipd
      
      # case-averted (per 100 000)
      
      case_avert <- (prevacc - postvacc)
      
      # case-averted IPD (per 100 000)
      
      case_avert_ipd <- (case_avert) * 0.011
      
      # case-averted OPD (per 100 000)
      
      case_avert_opd <- (case_avert) - (case_avert_ipd)
      
      # incremental cost of inpatients 
      
      incremental_cost_ipd   <- (v_cost + delivery_cost) * (vacc_pop) - (case_avert_ipd)*(facility_cost + dmc_ipd + dnmc_ipd)
      
      # incremental cost of outpatients
      
      incremental_cost_opd   <- (v_cost + delivery_cost) * (vacc_pop) - (case_avert_opd)*(facility_cost + dmc_opd + dnmc_opd)
      
      # incremental cost total (used for the final value as a delta C)
      
      incremental_cost_total <-  (v_cost + delivery_cost)* (vacc_pop) - (case_avert)*(facility_cost + dmc +dnmc + indirect) 
      
      # YLL calculation component: amount of years of life loss between life-expectancy and average age at death
      
      distance               <- (life_exp) - (age_death)
      
      # number of deaths in no-vaccination situation (inpatients)
      
      pre_death_ipd          <- (cfr_ipd)*(prevacc_ipd)
      
      # number of deaths in no-vaccination situation (outpatients)
      
      pre_death_opd          <- (cfr_opd)*(prevacc_opd)
      
      # number of deaths in no-vaccination situation (total)
      
      pre_death_total        <- (pre_death_ipd) + (pre_death_opd)
      
      # number of deaths post-vaccination situation (inpatients)
      
      post_death_ipd         <- (cfr_ipd)*(postvacc_ipd)
      
      # number of deaths post-vaccination situation (outpatients)
      
      post_death_opd         <- (cfr_opd)*(postvacc_opd)
      
      # number of deaths post-vaccination situation (total)
      
      post_death_total       <-  (post_death_ipd) + (post_death_opd) 
      
      # number of deaths averted (total)
      
      death_averted          <-  (pre_death_total) - (post_death_total)
      
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
      
      yld_pre_ipd            <- (prevacc_ipd)*(dw_ipd)*(illness_duration_ipd)
      
      # years of life with disease for outpatients in pre-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
      
      yld_pre_opd            <- (prevacc_opd)*(dw_opd)*(illness_duration_opd)
      
      # total YLD of pre-vaccination typhoid cases
      
      yld_pre_total          <- (yld_pre_ipd) + (yld_pre_opd)
      
      # total DALYs for pre-vaccination typhoid cases (YLL total + YLD total)
      
      dalys_pre_total        <- (yll_pre_total) + (yld_pre_total)
      
      # years of life with disease for inpatients in post-vacc situation (number of inpatient cases * disability weight * total illness duration presented as years)
      
      yld_post_ipd           <- (postvacc_ipd)*(dw_ipd)*(illness_duration_ipd)
      
      # years of life with disease for outpatients in post-vacc situation (number of outpatient cases * disability weight * total illness duration presented as years)
      
      yld_post_opd           <- (postvacc_opd)*(dw_opd)*(illness_duration_opd)
      
      # total YLD for post-vaccination typhoid cases (YLL total + YLD total)
      
      yld_post_total         <- (yld_post_ipd) + (yld_post_opd)
      
      # total DALYs for post-vaccination typhoid cases (YLL total + YLD total)
      
      dalys_post_total       <- (yll_post_total) + (yld_post_total)
      
      # YLD averted for inpatients
      
      yld_averted_ipd        <- (case_avert_ipd)* (dw_ipd) *(illness_duration_ipd)
      
      # YLD averted for outpatients
      
      yld_averted_opd        <- (case_avert_opd)* (dw_opd) *(illness_duration_opd)
      
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
param <- data.frame (v_cost               = c (mid = 3.274093,   low =  3.274093,   high =  3.274093), 
                     delivery_cost        = c (mid = 1.6296537,   low = 0.7567397,  high = 3.0051688),
                     tot_pop              = c (mid = 159831, low = 159831, high = 159831),
                     vacc_pop             = c (mid = 113420, low = 113420, high = 113420),
                     prevacc              = c (mid = 514.9838,     low = 377.4633,   high = 702.5708),
                     facility_cost        = c (mid = 96.80919,  low = 74.69514, high = 122.89700),
                     dmc_ipd              = c (mid = 145.122302, low = 2.415348, high = 964.741518),
                     dmc_opd              = c (mid = 60.1440990, low = 0.3186177, high = 536.6745889),
                     dmc                  = c (mid = 100.4637987, low = 0.7431742, high = 825.6930637),
                     dnmc_ipd             = c (mid = 39.39388,  low = 4.41536,  high = 144.22773), 
                     dnmc_opd             = c (mid = 12.3804299,  low = 0.1291842,  high = 93.3987597), 
                     dnmc                 = c (mid = 26.415861,  low = 1.207362,  high = 131.293796), 
                     indirect_ipd         = c (mid = 73.82400,  low = 30.99298, high = 145.19221),
                     indirect_opd         = c (mid = 71.06636,  low = 31.33397, high = 135.47664),
                     indirect             = c (mid = 72.67311,  low = 31.24815, high = 140.54778),
                     ve                   = c (mid = 0.8485963,   low = 0.4704072,   high = 0.9903380), 
                     dw_ipd               = c (mid = 0.2098825,  low = 0.1907369,  high = 0.2299149),
                     dw_opd               = c (mid = 0.05139286,  low = 0.03420124,  high = 0.07323530),
                     illness_duration_ipd = c (mid = 0.03848175,  low = 0.01618586, high = 0.07548403),
                     illness_duration_opd = c (mid = 0.04111426,  low = 0.01049451, high = 0.10599620),
                     cfr_ipd              = c (mid = 0.02708096,  low = 0.01326271,  high = 0.04792761),
                     cfr_opd              = c (mid = 0.01871302,  low = 0.01185146,  high = 0.02777084),
                     age_death            = c (mid = 7.547812,    low = 6.563674,  high = 8.625822),
                     life_exp             = c (mid = 72.68007,  low = 69.27367, high = 76.08483) )

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

# delete variables with no uncertainty ranges
param_table <- param_table[-c(1,3,4)]

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

param_table_1 <- param_table
Parameter <- c("Vaccine effectiveness", "CFR(OPD)", "Pre-vaccination disease burden", 
           "Direct medical cost (USD)", "Delivery cost (USD)", "Life-expectancy", 
           "Direct non-medical cost(USD)", "Indirect cost (USD)", "Facility cost (USD)",
           "Average age at death", "CFR (IPD)", "Illness duration (OPD)","Disease weight (OPD)","Illness duration (IPD)",
           "Disease weight (IPD)", "Direct medical cost (IPD)", "Direct medical cost (OPD)",
           "Direct non-medical cost (IPD)", "Direct non-medical cost (OPD)", "Indirect cost (IPD)","Indirect cost (OPD)")

param_table_2 <- cbind (Parameter, param_table_1) 


tornado_plot <- ggplot (param_table_2,
                        aes (reorder(Parameter, length_abs), 
                             ymin = icer_low, 
                             ymax = icer_high, 
                             color = Parameter)) +
  geom_linerange (size = 5) +
  geom_hline (yintercept = param_table_2$icer_mid [1], 
              linetype   = "dotted") +
  xlab ("Parameter") +
  ylab ("ICER (USD, cost per DALY averted)") + 
  coord_flip() +
  theme_bw () +
  theme(legend.text=element_text(size=5)) 

print (tornado_plot)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
toc ()
# main program (end)
# ------------------------------------------------------------------------------