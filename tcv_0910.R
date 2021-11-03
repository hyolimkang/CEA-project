# TCV cost-effectiveness

# things to do 
# update to latin hypercube sampling 

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
  my_params <-c ("v_cost",          # vaccine cost (including syringes): 2018 USD
                 "delivery_cost",   # delivery cost: 2018 USD
                 "tot_pop",         # total population in Phase 1 region
                 "vacc_pop",        # total number of vaccinated people
                 "postvacc_p1",     # post vaccination cases in Phase1 region
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
  # "postvacc_p1"
  
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
                           c (159831),
                           c (113420), 
                           c (23),
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

  compute_icer <- function (v_cost        = 2.96,
                            delivery_cost,
                            tot_pop       = 159831,
                            vacc_pop      = 113420,
                            postvacc_p1   = 23,
                            facility_cost,   
                            dmc,
                            dnmc,
                            indirect,
                            ve,
                            dw,
                            illness_duration,
                            cfr,
                            age_death,
                            life_exp) {
    
  totcost_vacc       <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1) * (facility_cost + dmc + dnmc + indirect)
  
  prevacc            <- (postvacc_p1)/(1-ve*(vacc_pop/tot_pop))
  
  totcost_unvacc     <- (prevacc)*(facility_cost + dmc + dnmc + indirect)
  
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
  
  
  return (list (totcost_vacc       = totcost_vacc, 
                totcost_unvacc     = totcost_unvacc, 
                incremental_cost   = incremental_cost,
                incremental_effect = incremental_effect,
                icer               = icer,
                prevacc            = prevacc,
                distance           = distance,
                pre_death          = pre_death,
                post_death         = post_death,
                yll_pre            = yll_pre,
                yll_post           = yll_post,
                yld_pre            = yld_pre,
                yld_post           = yld_post,
                yld_averted        = yld_averted,
                yll_averted        = yll_averted,
                daly_total         = daly_total,
                icer_daly          = icer_daly)) 
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
    
    icer_sample <- compute_icer (psa_sample$v_cost           [i],
                                 psa_sample$delivery_cost    [i],
                                 psa_sample$tot_pop          [i],
                                 psa_sample$vacc_pop         [i],
                                 psa_sample$postvacc_p1      [i],
                                 psa_sample$facility_cost    [i],   
                                 psa_sample$dmc              [i],
                                 psa_sample$dnmc             [i],
                                 psa_sample$indirect         [i],
                                 psa_sample$ve               [i],
                                 psa_sample$dw               [i],
                                 psa_sample$illness_duration [i],
                                 psa_sample$cfr              [i],
                                 psa_sample$age_death        [i],
                                 psa_sample$life_exp         [i])
    
    icer_dt [i, `:=` ( totcost_unvacc           = icer_sample$totcost_unvacc,
                       totcost_vacc             = icer_sample$totcost_vacc,
                      incremental_cost          = icer_sample$incremental_cost,
                      incremental_effect        = icer_sample$incremental_effect,
                      icer                      = icer_sample$icer,
                      prevacc                   = icer_sample$prevacc,
                      distance                  = icer_sample$distance,
                      pre_death                 = icer_sample$pre_death,
                      post_death                = icer_sample$post_death,
                      yll_pre                   = icer_sample$yll_pre,
                      yll_post                  = icer_sample$yll_post,
                      yld_pre                   = icer_sample$yld_pre,
                      yld_post                  = icer_sample$yld_post,
                      yld_averted               = icer_sample$yld_averted,
                      yll_averted               = icer_sample$yll_averted,
                      daly_total                = icer_sample$daly_total,
                      icer_daly                 = icer_sample$icer_daly
                     )]
  }
  
# central value

mean_icer         <- mean(icer_dt$icer)
mean_cost         <- mean(icer_dt$incremental_cost)
mean_case_averted <- mean(icer_dt$incremental_effect)
mean_prevacc      <- mean(icer_dt$prevacc)
mean_daly_total   <- mean(icer_dt$daly_total)
mean_yll_averted  <- mean(icer_dt$yll_averted)
mean_yld_averted  <- mean(icer_dt$yld_averted)
mean_icer_daly    <- mean(icer_dt$icer_daly)

# icer plot 
plot ( x = icer_dt$incremental_effect, y =icer_dt$incremental_cost)


# ceac

print(summary(icer_dt$icer_daly))
print(quantile(icer_dt$icer_daly, c(.5 , .025, .975 )))

probability_cea <- seq(from=0, to= 1, by= .01)
wtp<- quantile(icer_dt$icer_daly, probability_cea)

plot ( x= wtp, y = probability_cea)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
owsa <- function (v_cost           , 
                  delivery_cost    ,
                  tot_pop          ,
                  vacc_pop         ,
                  postvacc_p1      ,
                  facility_cost    ,
                  dmc              ,
                  dnmc             , 
                  indirect         ,
                  ve               , 
                  dw               ,
                  illness_duration ,
                  cfr              ,  
                  age_death        ,
                  life_exp         , 
                  parameters_change, 
                  change) {
  
  assign (x     = parameters_change, 
          value = get (parameter) * (1 - change))
  
  
  # formulas 
  
  totcost_vacc       <- (v_cost+ delivery_cost) * (vacc_pop) + (postvacc_p1) * (facility_cost + dmc + dnmc + indirect)
  
  prevacc            <- (postvacc_p1)/(1-ve*(vacc_pop/tot_pop))
  
  totcost_unvacc     <- (prevacc)*(facility_cost + dmc + dnmc + indirect)
  
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
  
  
  # return data table with values for tornado plot 
  
  
  return (list (icer_daly = icer_daly)) 
}
# end of function -- owsa

# one-way sensitivity analysis of parameters
# v_cost, delivery_cost, facility_cost, dmc, dnmc, indirect, ve, dw, illness_duration, cfr, age_death, life_exp: changing vars

parameters_change <- c ("v_cost",
                        "delivery_cost", 
                        "facility_cost", 
                        "dmc", 
                        "dnmc", 
                        "indirect", 
                        "ve", 
                        "dw", 
                        "illness_duration", 
                        "cfr", 
                        "age_death", 
                        "life_exp")

for (parameter in parameters_change) {
  
  change = 0.1  # +10% or -10%

  results_ll <- owsa (v_cost            = 2.96,
                      delivery_cost     = 1.49 ,
                      tot_pop           = 159831,
                      vacc_pop          = 113420,
                      postvacc_p1       = 23,
                      facility_cost     = 97.33 ,
                      dmc               = 183.07 ,
                      dnmc              = 30.13, 
                      indirect          = 73.09 ,
                      ve                = .81 ,
                      dw                = .052,
                      illness_duration  = .043  , 
                      cfr               = .019 ,
                      age_death         =  7.5  , 
                      life_exp          =  72.68 ,
                      parameters_change = parameter,
                      change            = change )
  
  results_ul <- owsa (v_cost            = 2.96,
                      delivery_cost     = 1.49 ,
                      tot_pop           = 159831,
                      vacc_pop          = 113420,
                      postvacc_p1       = 23,
                      facility_cost     = 97.33 ,
                      dmc               = 183.07 ,
                      dnmc              = 30.13, 
                      indirect          = 73.09 ,
                      ve                = .81 ,
                      dw                = .052,
                      illness_duration  = .043  , 
                      cfr               = .019 ,
                      age_death         =  7.5  , 
                      life_exp          =  72.68 ,
                      parameters_change = parameter, 
                      change            = -change )
  
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# main program
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# seed for random number generator
set.seed (1)  

# 2.5 % and 97.5% of each parameter (95% CI)

sapply (psa_sample, function(x) quantile(x, probs = seq(0, 1, 0.025)))


# data table for central ICER

parameters_change <- c ("v_cost",
                        "delivery_cost", 
                        "facility_cost", 
                        "dmc", 
                        "dnmc", 
                        "indirect", 
                        "ve", 
                        "dw", 
                        "illness_duration", 
                        "cfr", 
                        "age_death", 
                        "life_exp")

icer_dt_owsa_central <- data.table (run_id = 1 : length(1))

i <- 1

for (parameter in parameters_change) {
  
  
  owsa_sample <-        owsa(v_cost            = 2.96             ,
                             delivery_cost     = 1.49             ,
                             tot_pop           = 159831           ,
                             vacc_pop          = 113420           ,
                             postvacc_p1       = 23               ,
                             facility_cost     = 97.33            ,
                             dmc               = 183.07           ,
                             dnmc              = 30.13            , 
                             indirect          = 73.09            ,
                             ve                = .81              ,
                             dw                = .052             ,
                             illness_duration  = .043             , 
                             cfr               = .019             ,
                             age_death         =  7.5             , 
                             life_exp          =  72.68           ,
                             parameters_change = parameter       , 
                             change            = 0          ) 
  
  
  
  icer_dt_owsa_central [i, `:=` (icer_daly                 = owsa_sample$icer_daly
  )]
  
}
# data table for lower limit values
icer_dt_owsa_ll <- data.table (run_id = 1 : length(parameters_change))

i <- 0

for (parameter in parameters_change) {
  
i <- i + 1
    
  owsa_sample_ll <- owsa(v_cost        = 2.96             ,
                      delivery_cost    = 1.49             ,
                      tot_pop          = 159831           ,
                      vacc_pop         = 113420           ,
                      postvacc_p1      = 23               ,
                      facility_cost    = 97.33            ,
                      dmc              = 183.07           ,
                      dnmc             = 30.13            , 
                      indirect         = 73.09            ,
                      ve               = .81              ,
                      dw               = .052             ,
                      illness_duration = .043             , 
                      cfr              = .019             ,
                      age_death        =  7.5             , 
                      life_exp         =  72.68           ,
                      parameters_change = parameter       , 
                      change            = change          ) 
  

  
  
  icer_dt_owsa_ll [i, `:=` (icer_daly                 = owsa_sample_ll$icer_daly
  )]
  
}
# data table for upper limit values

icer_dt_owsa_ul <- data.table (run_id = 1 : length(parameters_change))

  i <- 0

for (parameter in parameters_change) {
  
  i <- i + 1
  
  owsa_sample_ul <- owsa(v_cost           = 2.96             ,
                         delivery_cost    = 1.49             ,
                         tot_pop          = 159831           ,
                         vacc_pop         = 113420           ,
                         postvacc_p1      = 23               ,
                         facility_cost    = 97.33            ,
                         dmc              = 183.07           ,
                         dnmc             = 30.13            , 
                         indirect         = 73.09            ,
                         ve               = .81              ,
                         dw               = .052             ,
                         illness_duration = .043             , 
                         cfr              = .019             ,
                         age_death        =  7.5             , 
                         life_exp         =  72.68           ,
                         parameters_change = parameter        , 
                         change            = -change          ) 
  
  
  icer_dt_owsa_ul [i, `:=` (icer_daly                 = owsa_sample_ul$icer_daly
  )]
  
}
  

# rename the icers in Lower bound and Upper bound   
  setnames(icer_dt_owsa_ll, "icer_daly", "icer_daly_ll")
  setnames(icer_dt_owsa_ul, "icer_daly", "icer_daly_ul")
  icer_dt_owsa_ul[, run_id := NULL]
  
# merge two data tables and remove the run_id to make a single dataset
  merge_icers  <- cbind(icer_dt_owsa_ll, icer_dt_owsa_central, icer_dt_owsa_ul)


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
# ------------------------------------------------------------------------
        
# ------------------------------------------------------------------------
# save icer_daly value to corresponding row and column
# note: icer (mid, low, high) columns are 3 columns after parameter (mid, low, high) columns
param_dt [i, j + 3] <- icer_daly
        
# ------------------------------------------------------------------------
        
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
  
  


 
 