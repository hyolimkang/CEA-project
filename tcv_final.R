# TCV cost-effectiveness

# things to do 
# update to latin hypercube sampling 
# ggplot for CEA Acceptability curve
# ggplot for tornado plot
# scenario analysis (how many scenarios?)
# correction factor
# any adjustments?


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

# function - computer icer for a given input parameter sample

compute_icer <- function (v_cost            = 2.96,
                          delivery_cost,
                          tot_pop           = 159831,
                          vacc_pop          = 113420,
                          postvacc          = c (ipd = 200, opd = 459, tot = 659),
                          facility_cost,
                          dmc               = c (ipd = dmc_ipd, opd = dmc_opd, tot = dmc)
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
  
  
  compute_icer <- function (param_dt)
    
  param_dt [, icer_ipd  := 0]
  param_dt [, icer_opd  := 0]
  param_dt [, icer_all  := 0]
  
  
  
  for(i in nrow(param_dt)) {
    
    # loop through mid, low, and high values (which are in columns 2, 3, 4)
    for (j in 2:4) {
    
    # select the mid value for each parameter
    v_cost           <- param_dt [parameter == "v_cost",           mid ]
    delivery_cost    <- param_dt [parameter == "delivery_cost",    mid ]
    tot_pop          <- param_dt [parameter == "tot_pop",          mid ]
    vacc_pop         <- param_dt [parameter == "vacc_pop",         mid ]
    postvacc_p1      <- param_dt [parameter == "postvacc_p1",      mid ]
    facility_cost    <- param_dt [parameter == "facility_cost",    mid ]
    dmc              <- param_dt [parameter == "dmc",              mid ]
    dnmc             <- param_dt [parameter == "dnmc",             mid ]
    indirect         <- param_dt [parameter == "indirect",         mid ]
    ve               <- param_dt [parameter == "ve",               mid ]
    dw               <- param_dt [parameter == "dw",               mid ]
    illness_duration <- param_dt [parameter == "illness_duration", mid ]
    cfr              <- param_dt [parameter == "cfr",              mid ]
    age_death        <- param_dt [parameter == "age_death",        mid ]
    life_exp         <- param_dt [parameter == "life_exp",         mid ]
    
  }
  
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
  
  incremental_daly_total <- (yld_averted_total) + (yll_averted_total)
  
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

# total number of PSA runs (or samples)
sample_size <- 1000

# create psa sample
psa_sample <- create_psa_sample (sample_n = sample_size)
psa_sample <- setDT (psa_sample)
psa_sample [, nsamp := NULL]

# create empty icer data table
icer_dt <- data.table (run_id = 1:sample_size)

# combine data tables -- run_id & psa sample
icer_dt <- cbind (icer_dt, psa_sample)

# create new columns
icer_dt [, c("incremental_cost")] <- 0  

# loop through psa sample to generate icers
for (i in 1:sample_size) {
  
  icer_sample <- compute_icer (psa_sample$v_cost               [i],
                               psa_sample$delivery_cost        [i],
                               psa_sample$tot_pop              [i],
                               psa_sample$vacc_pop             [i],
                               psa_sample$postvacc_p1_ipd      [i],
                               psa_sample$postvacc_p1_opd      [i],
                               psa_sample$postvacc_p1          [i],
                               psa_sample$facility_cost        [i],   
                               psa_sample$dmc_ipd              [i],
                               psa_sample$dmc_opd              [i],
                               psa_sample$dmc                  [i],
                               psa_sample$dnmc_ipd             [i],
                               psa_sample$dnmc_opd             [i],
                               psa_sample$dnmc                 [i],
                               psa_sample$indirect_ipd         [i],
                               psa_sample$indirect_opd         [i],
                               psa_sample$indirect             [i],
                               psa_sample$ve                   [i],
                               psa_sample$dw_ipd               [i],
                               psa_sample$dw_opd               [i],
                               psa_sample$illness_duration_ipd [i],
                               psa_sample$illness_duration_opd [i],
                               psa_sample$cfr_ipd              [i],
                               psa_sample$cfr_opd              [i],
                               psa_sample$age_death            [i],
                               psa_sample$life_exp             [i])
  
  icer_dt [i, `:=` ( totcost_unvacc                  = icer_sample$totcost_unvacc,
                     totcost_vacc_ipd                = icer_sample$totcost_vacc_ipd,
                     totcost_vacc_opd                = icer_sample$totcost_vacc_opd,
                     incremental_cost_ipd            = icer_sample$incremental_cost_ipd,
                     incremental_cost_opd            = icer_sample$incremental_cost_opd,
                     prevacc_ipd                     = icer_sample$prevacc_ipd,
                     prevacc_opd                     = icer_sample$prevacc_opd,
                     prevacc                         = icer_sample$prevacc,
                     incremental_cost_total          = icer_sample$incremental_cost_total,
                     case_averted_ipd                = icer_sample$case_averted_ipd,
                     case_averted_opd                = icer_sample$case_averted_opd,
                     case_averted_total              = icer_sample$case_averted_total,
                     icer_daly_ipd                   = icer_sample$icer_daly_ipd,
                     icer_daly_opd                   = icer_sample$icer_daly_opd,
                     icer_daly_total                 = icer_sample$icer_daly_total,
                     incremental_daly_total          = icer_sample$incremental_daly_total
                     
                     
  )]
}

# central value

mean_cost                   <- mean(icer_dt$incremental_cost_total)
mean_case_averted           <- mean(icer_dt$case_averted_total)
mean_prevacc                <- mean(icer_dt$prevacc)
mean_icer_daly_ipd          <- mean(icer_dt$icer_daly_ipd)
mean_icer_daly_opd          <- mean(icer_dt$icer_daly_opd)
mean_icer_daly_total        <- mean(icer_dt$icer_daly_total)
mean_incremental_daly_total <- mean(icer_dt$incremental_daly_total)

# icer plot 
plot ( x = icer_dt$incremental_daly_total, y =icer_dt$incremental_cost_total)


# ceac

print(summary(icer_dt$icer_daly_total))
print(quantile(icer_dt$icer_daly_total, c(.5 , .025, .975 )))

probability_cea <- seq(from=0, to= 1, by= .01)
wtp<- quantile(icer_dt$icer_daly_total, probability_cea)


plot ( x= wtp, y = probability_cea)

# cea plane with GDP per capita (2021: USD 2081.28 ~ 3*GDP: 6243.83)
ggplot() +
  geom_point(mapping = aes(x = incremental_daly_total,
                           y = incremental_cost_total),
             data = icer_dt)



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
      v_cost           <- param_dt [parameter == "v_cost",           mid ]
      delivery_cost    <- param_dt [parameter == "delivery_cost",    mid ]
      tot_pop          <- param_dt [parameter == "tot_pop",          mid ]
      vacc_pop         <- param_dt [parameter == "vacc_pop",         mid ]
      postvacc_p1      <- param_dt [parameter == "postvacc_p1",      mid ]
      facility_cost    <- param_dt [parameter == "facility_cost",    mid ]
      dmc              <- param_dt [parameter == "dmc",              mid ]
      dnmc             <- param_dt [parameter == "dnmc",             mid ]
      indirect         <- param_dt [parameter == "indirect",         mid ]
      ve               <- param_dt [parameter == "ve",               mid ]
      dw               <- param_dt [parameter == "dw",               mid ]
      illness_duration <- param_dt [parameter == "illness_duration", mid ]
      cfr              <- param_dt [parameter == "cfr",              mid ]
      age_death        <- param_dt [parameter == "age_death",        mid ]
      life_exp         <- param_dt [parameter == "life_exp",         mid ]
      
      # ------------------------------------------------------------------------
      # assign value for a specific parameter to mid/low/high values
      assign (x     = param_dt [i, parameter], 
              value = param_dt [i, j, with = FALSE] )
      # ------------------------------------------------------------------------
      
      # ------------------------------------------------------------------------
      # formulas for icer
      totcost_vacc       <- (v_cost + delivery_cost) * (vacc_pop) + 
        (postvacc_p1) * (facility_cost + dmc + dnmc + indirect)
      
      prevacc            <- (postvacc_p1)/(1-ve*(vacc_pop/tot_pop))
      
      totcost_unvacc     <- (prevacc)*(facility_cost+ dmc + dnmc + indirect)
      
      incremental_cost   <- (totcost_vacc) - (totcost_unvacc)
      
      incremental_effect <- (prevacc)      - (postvacc_p1) 
      
      icer               <- incremental_cost / incremental_effect
      
      distance           <- (life_exp) - (age_death)
      
      pre_death          <- (cfr)*(prevacc)
      
      post_death         <- (cfr)*(postvacc_p1)
      
      yll_pre            <- (pre_death)*(distance)
      
      yll_post           <- (post_death)*(distance)
      
      yld_pre            <- (prevacc)*(dw)*(illness_duration)
      
      yld_post           <- (postvacc_p1)*(dw)*(illness_duration)
      
      yld_averted        <- (prevacc - postvacc_p1)* (dw) *(illness_duration)
      
      yll_averted        <- (yll_pre)  -(yll_post)
      
      daly_total         <- (yll_averted)  + (yld_averted)
      
      icer_daly          <- incremental_cost / daly_total
      
      # ----------------------------------------------------------------------
      
      # ------------------------------------------------------------------------
      # save icer_daly value to corresponding row and column
      # note: icer (mid, low, high) columns are 3 columns after parameter (mid, low, high) columns
      param_dt [i, j + 3] <- icer_daly
      
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

# initialise parameter data table
param <- data.frame (v_cost           = c (mid = 2.96,   low = 2.96,   high = 2.96    ), 
                     delivery_cost    = c (mid = 1.49,   low = 0.246,  high = 4.752   ),
                     tot_pop          = c (mid = 159831, low = 159831, high = 159831  ),
                     vacc_pop         = c (mid = 113420, low = 113420, high = 113420  ),
                     postvacc_p1      = c (mid = 23,     low = 23,     high = 23      ),
                     facility_cost    = c (mid = 97.33,  low = 16.676, high = 301.209 ),
                     dmc              = c (mid = 183.07, low = 18.507, high = 532.241 ),
                     dnmc             = c (mid = 30.13,  low = 4.361,  high = 118.718 ), 
                     indirect         = c (mid = 73.09,  low = 32.370, high = 135.262 ),
                     ve               = c (mid = 0.81,   low = 0.55,   high = 0.916   ), 
                     dw               = c (mid = 0.052,  low = 0.032,  high = 0.078   ),
                     illness_duration = c (mid = 0.043,  low = 0.0392, high = 0.0469  ),
                     cfr              = c (mid = 0.019,  low = 0.011,  high = 0.0441  ),
                     age_death        = c (mid = 7.5,    low = 1.840,  high = 13.95   ),
                     life_exp         = c (mid = 72.68,  low = 69.576, high = 78.70   ) )

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
tornado_plot <- ggplot (param_table,
                        aes (parameter, 
                             ymin = icer_low, 
                             ymax = icer_high, 
                             color = parameter)) +
  geom_linerange (size = 5) +
  geom_hline (yintercept = param_table$icer_mid [1], 
              linetype   = "dotted") +
  xlab ("paramter") +
  ylab ("ICER (cost per DALY averted)") + 
  coord_flip() +
  theme_bw () 

print (tornado_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# main program (end)
# ------------------------------------------------------------------------------





