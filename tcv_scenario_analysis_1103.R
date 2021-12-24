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

# clear workspace
rm (list = ls())

  
# scenario analysis (phase2 area)------------------------------------------------------------------------------
  
  # parameters
create_psa_sample <- function (sample_n) {
  
 my_params <-c ("v_cost",          # vaccine cost (including syringes): 2018 USD
                 "delivery_cost",   # delivery cost: 2018 USD
                 "tot_pop",         # total population in Phase 1 region
                 "vacc_pop",        # total number of vaccinated people
                 "postvacc_p2",     # post vaccination cases in Phase1 region
                 "facility_cost",   # health facility cost : 2021 USD
                 "dmc",             # direct medical cost : 2021 USD
                 "dnmc",            # direct non-medical cost : 2021 USD
                 "indirect",        # indirect cost : 2021 USD
                 "ve",              # vaccine-effectiveness
                 "dw",              # disability weight
                 "illness_duration",# illness duration: days
                 "cfr",             # case fatality rate
                 "age_death",       # average age at death
                 "life_exp"   )      # reference life expectancy
  
  
  # constant value
  # "postvacc"        # post-vaccination case
  # "vacc_pop"
  # "postvacc_p2"
  
  # trickle down uncertainty
  # "prevacc",         # pre-vaccination case
  # "totcost_vacc",    # total cost of vaccinated group
  # "totcost_unvacc")  # total cost of unvaccinated group
  
  
  # distributions for parameters
  my_dists <- c("constant",
                "log-normal",
                "constant",
                "constant",
                "constant",
                "log-normal",
                "log-normal",
                "log-normal",
                "log-normal",
                "truncated-normal",
                "truncated-normal",
                "log-normal",
                "truncated-normal",
                "truncated-normal",
                "truncated-normal")
  
  
  # parameter types
  my_parameterization_types <-c ("val",
                                 "mean, sd",
                                 "val",
                                 "val",
                                 "val",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul",
                                 "mean, sd",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul",
                                 "mean, sd, ll, ul")
  
  
  # values for parameters
  my_dists_params <- list( c (2.96),
                           c (1.49, 1.38),
                           c (331418),
                           c (113420), 
                           c (54),
                           c (97.33,  85.33),
                           c (183.07, 229.31),
                           c (30.13,  31.06),
                           c (73.09, 28.01),
                           c (.81,     2.05, .54, .926),
                           c (.052,  1.37, .031, .079),
                           c (.043, .002),
                           c (.019, .72, .01, .045),
                           c (7.5, 4.25,  0.75 , 15),
                           c (72.68, 3.07, 69.27, 81.54))
  
  # input samples for PSA
  my_psa <- gen_psa_samp (params                 = my_params,
                          dists                  = my_dists,
                          parameterization_types = my_parameterization_types,
                          dists_params           = my_dists_params,
                          n                      = sample_n)
  
  return (my_psa)
  
}

# function - computer icer for a given input parameter sample

# ------------------------------------------------------------------------------
# Latin Hypercube Sampling
# ------------------------------------------------------------------------------

set.seed (3)
runs <- 2000

# a design with n samples from k parameters
A <- randomLHS (n = runs, 
                k = 26) 
# It is common to transform the margins of the design (the columns) 
# into other distributions
# for gamma dist: a= shape (m^2/sigma^2), b= rate (m/sigma^2)
# for truncated norm dist: 

lhs_sample_s <- matrix (nrow = nrow(A), ncol = ncol(A))
lhs_sample_s [,1]  <- 2.96
lhs_sample_s [,2]  <- qgamma (p = A[,2], shape = ((1.49/1.38)^2), rate = (1.49/(1.38)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,3]  <- 331418
lhs_sample_s [,4]  <- 113420
lhs_sample_s [,5]  <- 458
lhs_sample_s [,6]  <- 1140
lhs_sample_s [,7]  <- 1598
lhs_sample_s [,8]  <- qgamma (p = A[,8], shape = ((97.33/85.33)^2), rate = (97.33/(85.33)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,9]  <- qgamma (p = A[,9], shape = ((234.77/265.99)^2), rate = (234.77/(265.99)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,10] <- qgamma (p = A[,10], shape = ((115.38/149.92)^2), rate = (115.38/(149.92)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,11] <- qgamma (p = A[,11], shape = ((183.34/229.85)^2), rate = (183.34/(229.85)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,12] <- qgamma (p = A[,12], shape = ((48.59/37.40)^2), rate = (48.59/(37.40)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,13] <- qgamma (p = A[,13],  shape = ((21.49/25.90)^2), rate = (21.49/(25.90)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,14] <- qgamma (p = A[,14],  shape = ((36.91/35.39)^2), rate = (36.91/(35.39)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,15] <- qgamma (p = A[,15],  shape = ((77.54/29.51)^2), rate = (77.54/(29.51)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,16] <- qgamma (p = A[,16],   shape = ((74.29/26.90)^2), rate = (74.29/(26.90)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,17] <- qgamma (p = A[,17],  shape = ((76.13/28.25)^2), rate = (76.13/(28.25)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,18] <- qtruncnorm (p = A[,18], a = 0.54 , b = 0.926 , mean = 0.81, sd = 2.05)
lhs_sample_s [,19] <- qtruncnorm (p = A[,19], a = 0.19,  b = 0.23  , mean = 0.21, sd = 0.052) 
lhs_sample_s [,20] <- qtruncnorm (p = A[,20], a = 0.031, b = 0.079 , mean = 0.052, sd = 1.37)
lhs_sample_s [,21] <- qgamma (p = A[,21], shape = ((0.04/0.153)^2), rate = (0.04/(0.153)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,22] <- qgamma (p = A[,22], shape = ((0.046/0.025)^2), rate = (0.046/(0.025)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,23] <- qtruncnorm (p = A[,23], a = 0.02, b = 0.036, mean = 0.028, sd = 0.211)
lhs_sample_s [,24] <- qtruncnorm (p = A[,24], a = 0.01, b = 0.045, mean = 0.019, sd = 0.72)
lhs_sample_s [,25] <- qgamma (p = A[,25], shape = ((7.56/6.53)^2), rate = (7.56/(6.53)^2), 
                            lower.tail = TRUE, log.p = FALSE)
lhs_sample_s [,26] <- qnorm (p = A[,26], mean = 72.68, sd = 3.07, lower.tail = TRUE, log.p = FALSE)

# change the name of the lhs_sample
cols <- c(  "v_cost", 
            "delivery_cost",    
            "tot_pop",            
            "vacc_pop",             
            "postvacc_p1_ipd",      
            "postvacc_p1_opd",     
            "postvacc_p1",          
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
colnames (lhs_sample_s) <- cols
# change the matrix into the data table. (R doesn't recognize atomic values inside matrix)
lhs_sample_s <- as.data.table (lhs_sample_s)

# function - computer icer for a given input parameter sample

icer_scenario <- function (v_cost           = 2.96,
                          delivery_cost,
                          tot_pop           = 331418,
                          vacc_pop          = 113420,
                          postvacc_p1_ipd   = 458,
                          postvacc_p1_opd   = 1140,
                          postvacc_p1       = 1598,
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
  
  
  totcost_vacc_ipd   <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_ipd) * (facility_cost + dmc_ipd + dnmc_ipd + indirect_ipd)
  
  totcost_vacc_opd   <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_opd) * (facility_cost + dmc_opd + dnmc_opd + indirect_opd)
  
  totcost_vacc_tot   <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1_ipd) * (facility_cost + dmc_ipd + dnmc_ipd + indirect_ipd) + (postvacc_p1_opd) * (facility_cost + dmc_opd + dnmc_opd + indirect_opd)
  
  prevacc_ipd        <- (postvacc_p1_ipd)/(1-ve*(vacc_pop/tot_pop))
  
  prevacc_opd        <- (postvacc_p1_opd)/(1-ve*(vacc_pop/tot_pop))
  
  prevacc            <- (prevacc_ipd)  + (prevacc_opd)
  
  totcost_unvacc     <- (prevacc)*(facility_cost + dmc + dnmc + indirect)
  
  incremental_cost_ipd   <- (totcost_vacc_ipd) - (totcost_unvacc)
  
  incremental_cost_opd   <- (totcost_vacc_opd) - (totcost_unvacc)
  
  incremental_cost_total <- (totcost_vacc_tot) - (totcost_unvacc)
  
  case_averted_ipd      <- (prevacc_ipd)      - (postvacc_p1_ipd)
  
  case_averted_opd      <- (prevacc_opd)      - (postvacc_p1_opd)
  
  case_averted_total    <- (case_averted_ipd) + (case_averted_opd)
  
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
  
  yll_pre_ipd            <- (pre_death_ipd)*(distance)
  
  yll_pre_opd            <- (pre_death_opd)*(distance)
  
  yll_pre_total          <- (yll_pre_ipd)  + (yll_pre_opd)
  
  yll_post_ipd           <- (post_death_ipd)*(distance)
  
  yll_post_opd           <- (post_death_opd)*(distance)
  
  yll_post_total         <- (yll_post_ipd) + (yll_post_opd)
  
  yld_pre_ipd            <- (prevacc_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  yld_pre_opd            <- (prevacc_opd)*(dw_opd)*(illness_duration_opd)
  
  yld_pre_total          <- (yld_pre_ipd) + (yld_pre_opd)
  
  yld_post_ipd           <- (postvacc_p1_ipd)*(dw_ipd)*(illness_duration_ipd)
  
  yld_post_opd           <- (postvacc_p1_opd)*(dw_opd)*(illness_duration_opd)
  
  yld_post_total         <- (yld_post_ipd) + (yld_post_opd)
  
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
  
  icer_daly_total          <- incremental_cost_total / incremental_daly_total
  
  
  
  return (list (totcost_vacc_ipd              = totcost_vacc_ipd,
                totcost_vacc_opd              = totcost_vacc_opd,
                totcost_vacc_tot              = totcost_vacc_tot,
                totcost_unvacc                = totcost_unvacc, 
                incremental_cost_ipd          = incremental_cost_ipd,
                incremental_cost_opd          = incremental_cost_opd,
                incremental_cost_total        = incremental_cost_total,
                case_averted_ipd              = case_averted_ipd,
                case_averted_opd              = case_averted_opd,
                case_averted_total            = case_averted_total,
                icer_ipd                      = icer_ipd,
                icer_opd                      = icer_opd,
                icer_tot                      = icer_tot,
                distance                      = distance,
                pre_death_ipd                 = pre_death_ipd,
                pre_death_opd                 = pre_death_opd,
                pre_death_total               = pre_death_total,
                post_death_ipd                = post_death_ipd,
                post_death_opd                = post_death_opd,
                post_death_total              = post_death_total,
                prevacc_ipd                   = prevacc_ipd,
                prevacc_opd                   = prevacc_opd,
                prevacc                       = prevacc,
                distance                      = distance,
                yll_pre_ipd                   = yll_pre_ipd,
                yll_pre_opd                   = yll_pre_opd,
                yll_pre_total                 = yll_pre_total,
                yll_pre_total                 = yll_pre_total,
                yll_post_ipd                  = yll_post_ipd,
                yll_post_opd                  = yll_post_total,
                yll_post_total                = yll_post_total,
                yld_pre_ipd                   = yll_pre_ipd,
                yld_pre_opd                   = yld_pre_opd,
                yld_pre_total                 = yld_pre_total,
                yld_post_ipd                  = yld_post_ipd,
                yld_post_opd                  = yld_post_opd,
                yld_post_total                = yld_post_total,
                yld_averted_ipd               = yld_averted_ipd,
                yld_averted_opd               = yld_averted_opd,
                yld_averted_total             = yld_averted_total,
                yll_averted_ipd               = yll_averted_ipd,
                yll_averted_opd               = yll_averted_opd,
                yll_averted_total             = yll_averted_total,
                incremental_daly_ipd          = incremental_daly_ipd,
                incremental_daly_opd          = incremental_daly_opd,
                incremental_daly_total        = incremental_daly_total,
                icer_daly_ipd                 = icer_daly_ipd,
                icer_daly_opd                 = icer_daly_opd,
                icer_daly_total               = icer_daly_total)) 
}

# start time
tic ()
print (Sys.time ())

# create empty icer data table 
icer_dt_s <- data.table (run_id = 1:runs)

# combine data tables -- run_id & lhs sample
icer_dt_s <- cbind (icer_dt_s, lhs_sample_s)

# loop through psa sample to generate icers
for (i in 1:runs) {
  
  icer_sample_s <- icer_scenario (lhs_sample_s$v_cost              [i],
                               lhs_sample_s$delivery_cost        [i],
                               lhs_sample_s$tot_pop              [i],
                               lhs_sample_s$vacc_pop             [i],
                               lhs_sample_s$postvacc_p1_ipd      [i],
                               lhs_sample_s$postvacc_p1_opd      [i],
                               lhs_sample_s$postvacc_p1          [i],
                               lhs_sample_s$facility_cost        [i],   
                               lhs_sample_s$dmc_ipd              [i],
                               lhs_sample_s$dmc_opd              [i],
                               lhs_sample_s$dmc                  [i],
                               lhs_sample_s$dnmc_ipd             [i],
                               lhs_sample_s$dnmc_opd             [i],
                               lhs_sample_s$dnmc                 [i],
                               lhs_sample_s$indirect_ipd         [i],
                               lhs_sample_s$indirect_opd         [i],
                               lhs_sample_s$indirect             [i],
                               lhs_sample_s$ve                   [i],
                               lhs_sample_s$dw_ipd               [i],
                               lhs_sample_s$dw_opd               [i],
                               lhs_sample_s$illness_duration_ipd [i],
                               lhs_sample_s$illness_duration_opd [i],
                               lhs_sample_s$cfr_ipd              [i],
                               lhs_sample_s$cfr_opd              [i],
                               lhs_sample_s$age_death            [i],
                               lhs_sample_s$life_exp             [i])
  
  icer_dt_s [i, `:=` ( totcost_unvacc                = icer_sample_s$totcost_unvacc,
                     totcost_vacc_ipd                = icer_sample_s$totcost_vacc_ipd,
                     totcost_vacc_opd                = icer_sample_s$totcost_vacc_opd,
                     incremental_cost_ipd            = icer_sample_s$incremental_cost_ipd,
                     incremental_cost_opd            = icer_sample_s$incremental_cost_opd,
                     prevacc_ipd                     = icer_sample_s$prevacc_ipd,
                     prevacc_opd                     = icer_sample_s$prevacc_opd,
                     prevacc                         = icer_sample_s$prevacc,
                     incremental_cost_total          = icer_sample_s$incremental_cost_total,
                     case_averted_ipd                = icer_sample_s$case_averted_ipd,
                     case_averted_opd                = icer_sample_s$case_averted_opd,
                     case_averted_total              = icer_sample_s$case_averted_total,
                     icer_daly_ipd                   = icer_sample_s$icer_daly_ipd,
                     icer_daly_opd                   = icer_sample_s$icer_daly_opd,
                     icer_daly_total                 = icer_sample_s$icer_daly_total,
                     incremental_daly_total          = icer_sample_s$incremental_daly_total
                     
                     
  )]
}

# central value

mean_cost                   <- mean(icer_dt_s$incremental_cost_total)
mean_case_averted           <- mean(icer_dt_s$case_averted_total)
mean_prevacc                <- mean(icer_dt_s$prevacc)
mean_icer_daly_ipd          <- mean(icer_dt_s$icer_daly_ipd)
mean_icer_daly_opd          <- mean(icer_dt_s$icer_daly_opd)
mean_icer_daly_total        <- mean(icer_dt_s$icer_daly_total)
mean_incremental_daly_total <- mean(icer_dt_s$incremental_daly_total)

# icer plot 
plot ( x = icer_dt$incremental_daly_total, y =icer_dt$incremental_cost_total)


# ceac

print(summary(icer_dt_s$icer_daly_total))
print(quantile(icer_dt_S$icer_daly_total, c(.5 , .025, .975 )))

probability_cea <- seq(from=0, to= 1, by= .01)

wtp_scenario <- quantile(icer_dt$icer_daly_total, probability_cea)

wtp_dt_scenario <- as.data.table(wtp_scenario)

wtp_prob_scenario <- cbind (wtp_dt_scenario, probability_cea)

# keep only positive wtp values

wtp_post_scenario <- wtp_prob_scenario %>% filter(wtp_scenario > 0)

# ceac ggplot

ceac_scenario <- ggplot(data = wtp_post_scenario, aes(x = wtp_scenario, y = probability_cea)) +
  geom_line(color = 'red', size = 0.8) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))



