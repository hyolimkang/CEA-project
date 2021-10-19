# TCV cost-effectiveness

# load libraries
library (dampack)
library (data.table)
library (ggplot2)
library (tictoc)
library (caret)
library (rriskDistributions)
library (BCEA)

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
                           c (.052,  .031, .079, 1.37),
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
  
  temp_distance      <- (life_exp) - (age_death)
  
  pre_death          <- (cfr)*(prevacc)
  
  post_death         <- (cfr)*(postvacc_p1)
  
  yll_pre            <- (pre_death)*(temp_distance)
  
  yll_post           <- (post_death)*(temp_distance)
  
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
                temp_distance      = temp_distance,
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
                     temp_distance             = icer_sample$temp_distance,
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
owsa <- function (v_cost           = 2.96,
                  delivery_cost    = 1.49 ,
                  tot_pop          = 159831,
                  vacc_pop         = 113420,
                  postvacc_p1      = 23,
                  facility_cost    = 97.33 ,
                  dmc              = 183.07 ,
                  dnmc             = 30.13, 
                  indirect         = 73.09 ,
                  ve               = .81 ,
                  dw               = .052,
                  illness_duration = .043  , 
                  cfr              = .019 ,
                  age_death        =  7.5  , 
                  life_exp         =  72.68 ,
                  parameters_change, 
                  change) {
  
  assign (x     = parameters_change, 
          value = get (parameters_change) * (1 - change))
  
  
  # lower bound
  totcost_vacc       <- (v_cost + delivery_cost) * (vacc_pop) + (postvacc_p1) * (facility_cost + dmc + dnmc + indirect)
  
  prevacc            <- (postvacc_p1)/(1-ve*(vacc_pop/tot_pop))
  
  totcost_unvacc     <- (prevacc)*(facility_cost + dmc + dnmc + indirect)
  
  incremental_cost   <- (totcost_vacc) - (totcost_unvacc)
  
  incremental_effect <- (prevacc)      - (postvacc_p1) 
  
  icer               <- incremental_cost / incremental_effect
  
  temp_distance      <- (life_exp) - (age_death)
  
  pre_death          <- (cfr)*(prevacc)
  
  post_death         <- (cfr)*(postvacc_p1)
  
  yll_pre            <- (pre_death)*(temp_distance)
  
  yll_post           <- (post_death)*(temp_distance)
  
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
                temp_distance      = temp_distance,
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
# end of function -- owsa

# one-way sensitivity analysis of parameters
# v_cost, delivery_cost, facility_cost, dmc, dnmc, indirect, ve, dw, illness_duration, cfr, age_death, life_exp: changing vars

# changing vaccine cost (+/- 10%) with all other variables constant

parameters_change1 <- c ("v_cost")
                        
for (parameter in parameters_change1) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_vcost <- owsa (v_cost            = 2.96,
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
  
  results_ul_vcost <- owsa (v_cost            = 2.96,
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

# changing delivery cost (+/- 10%) with all other variables constant

parameters_change2 <-   c("delivery_cost") 

for (parameter in parameters_change2) {
  
change = 0.1  # +10% or -10%

results_ll_deliv <- owsa (v_cost      = 2.96,
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

results_ul_deliv <- owsa (v_cost      = 2.96,
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

# changing facility cost (+/- 10%) with all other variables constant

parameters_change3 <-   c("facility_cost") 

for (parameter in parameters_change3) {
  
change = 0.1  # +10% or -10%

results_ll_facil <- owsa (v_cost      = 2.96,
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

results_ul_facil <- owsa (v_cost            = 2.96,
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

# changing direct medical cost (+/- 10%) with all other variables constant

parameters_change4 <-     c("dmc")

for (parameter in parameters_change4) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_dmc <- owsa (  v_cost            = 2.96,
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
  
  results_ul_dmc  <- owsa (v_cost             = 2.96,
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

# changing direct non-medical cost (+/- 10%) with all other variables constant

parameters_change5<-      c("dnmc") 

for (parameter in parameters_change5) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_dnmc <- owsa (v_cost           = 2.96,
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
  
  results_ul_dnmc  <- owsa (v_cost             = 2.96,
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

# changing indirect cost (+/- 10%) with all other variables constant

parameters_change6 <-     c("indirect") 

for (parameter in parameters_change6) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_indirect <- owsa (v_cost            = 2.96,
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
  
  results_ul_indirect  <- owsa (v_cost             = 2.96,
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

# changing indirect cost (+/- 10%) with all other variables constant
parameters_change7 <-      c("ve")

for (parameter in parameters_change7) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_ve<- owsa (       v_cost            = 2.96,
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
  
  results_ul_ve  <- owsa (v_cost            = 2.96,
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

# changing indirect cost (+/- 10%) with all other variables constant
parameters_change8 <-      c("dw") 

for (parameter in parameters_change8) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_dw<- owsa (       v_cost            = 2.96,
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
  
  results_ul_dw  <- owsa (v_cost            = 2.96,
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

# changing indirect cost (+/- 10%) with all other variables constant
parameters_change9 <-     c("illness_duration") 

for (parameter in parameters_change9) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_illness<- owsa (       v_cost            = 2.96,
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
  
  results_ul_illness <- owsa (v_cost            = 2.96,
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

# changing CFR (+/- 10%) with all other variables constant
parameters_change10 <-     c("cfr")

for (parameter in parameters_change10) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_cfr <- owsa (       v_cost            = 2.96,
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
  
  results_ul_cfr <- owsa (v_cost            = 2.96,
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

# changing age at death (+/- 10%) with all other variables constant

parameters_change11 <-     c("age_death")

for (parameter in parameters_change11) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_age <- owsa (       v_cost            = 2.96,
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
  
  results_ul_age <- owsa (v_cost            = 2.96,
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

# changing life expectancy (+/- 10%) with all other variables constant

parameters_change12 <-     c("life_exp")

for (parameter in parameters_change12) {
  
  change = 0.1  # +10% or -10%
  
  results_ll_le <- owsa (        v_cost            = 2.96,
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
  
  results_ul_le <- owsa (v_cost            = 2.96,
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
  
  
  
  
  icer_dt_owsa_central [i, `:=` ( totcost_unvacc          = owsa_sample$totcost_unvacc,
                                totcost_vacc              = owsa_sample$totcost_vacc,
                                incremental_cost          = owsa_sample$incremental_cost,
                                incremental_effect        = owsa_sample$incremental_effect,
                                icer                      = owsa_sample$icer,
                                prevacc                   = owsa_sample$prevacc,
                                temp_distance             = owsa_sample$temp_distance,
                                pre_death                 = owsa_sample$pre_death,
                                post_death                = owsa_sample$post_death,
                                yll_pre                   = owsa_sample$yll_pre,
                                yll_post                  = owsa_sample$yll_post,
                                yld_pre                   = owsa_sample$yld_pre,
                                yld_post                  = owsa_sample$yld_post,
                                yld_averted               = owsa_sample$yld_averted,
                                yll_averted               = owsa_sample$yll_averted,
                                daly_total                = owsa_sample$daly_total,
                                icer_daly                 = owsa_sample$icer_daly
  )]
  
}


# data table for p1 change (+/- 10 %)
icer_dt_owsa_p1_ll <- data.table (run_id = 1 : length(parameters_change1))

i <- 0

for (parameter in parameters_change1) {
  
  i <- i + 1
  
  # lower limit values
  
  owsa_sample_p1_ll <- owsa(v_cost        = 2.96             ,
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
  
  
  
  
  icer_dt_owsa_p1_ll [i, `:=` ( totcost_unvacc      = owsa_sample_p1_ll$totcost_unvacc,
                          totcost_vacc              = owsa_sample_p1_ll$totcost_vacc,
                          incremental_cost          = owsa_sample_p1_ll$incremental_cost,
                          incremental_effect        = owsa_sample_p1_ll$incremental_effect,
                          icer                      = owsa_sample_p1_ll$icer,
                          prevacc                   = owsa_sample_p1_ll$prevacc,
                          temp_distance             = owsa_sample_p1_ll$temp_distance,
                          pre_death                 = owsa_sample_p1_ll$pre_death,
                          post_death                = owsa_sample_p1_ll$post_death,
                          yll_pre                   = owsa_sample_p1_ll$yll_pre,
                          yll_post                  = owsa_sample_p1_ll$yll_post,
                          yld_pre                   = owsa_sample_p1_ll$yld_pre,
                          yld_post                  = owsa_sample_p1_ll$yld_post,
                          yld_averted               = owsa_sample_p1_ll$yld_averted,
                          yll_averted               = owsa_sample_p1_ll$yll_averted,
                          daly_total                = owsa_sample_p1_ll$daly_total,
                          icer_daly                 = owsa_sample_p1_ll$icer_daly
  )]
  
}
  # upper limit values
icer_dt_owsa_p1_ul <- data.table (run_id = 1 : length(parameters_change1))

i <- 0

for (parameter in parameters_change1) {
  
i <- i + 1

  owsa_sample_p1_ul <- owsa(v_cost        = 2.96             ,
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
  
  
  icer_dt_owsa_p1_ul [i, `:=` ( totcost_unvacc         = owsa_sample_p1_ul$totcost_unvacc,
                          totcost_vacc              = owsa_sample_p1_ul$totcost_vacc,
                          incremental_cost          = owsa_sample_p1_ul$incremental_cost,
                          incremental_effect        = owsa_sample_p1_ul$incremental_effect,
                          icer                      = owsa_sample_p1_ul$icer,
                          prevacc                   = owsa_sample_p1_ul$prevacc,
                          temp_distance             = owsa_sample_p1_ul$temp_distance,
                          pre_death                 = owsa_sample_p1_ul$pre_death,
                          post_death                = owsa_sample_p1_ul$post_death,
                          yll_pre                   = owsa_sample_p1_ul$yll_pre,
                          yll_post                  = owsa_sample_p1_ul$yll_post,
                          yld_pre                   = owsa_sample_p1_ul$yld_pre,
                          yld_post                  = owsa_sample_p1_ul$yld_post,
                          yld_averted               = owsa_sample_p1_ul$yld_averted,
                          yll_averted               = owsa_sample_p1_ul$yll_averted,
                          daly_total                = owsa_sample_p1_ul$daly_total,
                          icer_daly                 = owsa_sample_p1_ul$icer_daly
  )]
}

# append ul and ll of p1
append_icer_p1 <- data.table :: rbindlist(list(icer_dt_owsa_p1_ll, icer_dt_owsa_p1_ul))

# data table for p2 change (+/- 10 %)
icer_dt_owsa_p2_ll <- data.table (run_id = 1 : length(parameters_change2))

i <- 0

for (parameter in parameters_change2) {
  
  i <- i + 1
  
  # lower limit values
  
  owsa_sample_p2_ll <- owsa(v_cost           = 2.96             ,
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
  
  
  
  
  icer_dt_owsa_p2_ll [i, `:=` ( totcost_unvacc            = owsa_sample_p2_ll$totcost_unvacc,
                                totcost_vacc              = owsa_sample_p2_ll$totcost_vacc,
                                incremental_cost          = owsa_sample_p2_ll$incremental_cost,
                                incremental_effect        = owsa_sample_p2_ll$incremental_effect,
                                icer                      = owsa_sample_p2_ll$icer,
                                prevacc                   = owsa_sample_p2_ll$prevacc,
                                temp_distance             = owsa_sample_p2_ll$temp_distance,
                                pre_death                 = owsa_sample_p2_ll$pre_death,
                                post_death                = owsa_sample_p2_ll$post_death,
                                yll_pre                   = owsa_sample_p2_ll$yll_pre,
                                yll_post                  = owsa_sample_p2_ll$yll_post,
                                yld_pre                   = owsa_sample_p2_ll$yld_pre,
                                yld_post                  = owsa_sample_p2_ll$yld_post,
                                yld_averted               = owsa_sample_p2_ll$yld_averted,
                                yll_averted               = owsa_sample_p2_ll$yll_averted,
                                daly_total                = owsa_sample_p2_ll$daly_total,
                                icer_daly                 = owsa_sample_p2_ll$icer_daly
  )]
  
}
# upper limit values
icer_dt_owsa_p2_ul <- data.table (run_id = 1 : length(parameters_change2))

i <- 0

for (parameter in parameters_change2) {
  
  i <- i + 1
  
  owsa_sample_p2_ul <- owsa(v_cost           = 2.96             ,
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
  
  
  icer_dt_owsa_p2_ul [i, `:=` ( totcost_unvacc            = owsa_sample_p2_ul$totcost_unvacc,
                                totcost_vacc              = owsa_sample_p2_ul$totcost_vacc,
                                incremental_cost          = owsa_sample_p2_ul$incremental_cost,
                                incremental_effect        = owsa_sample_p2_ul$incremental_effect,
                                icer                      = owsa_sample_p2_ul$icer,
                                prevacc                   = owsa_sample_p2_ul$prevacc,
                                temp_distance             = owsa_sample_p2_ul$temp_distance,
                                pre_death                 = owsa_sample_p2_ul$pre_death,
                                post_death                = owsa_sample_p2_ul$post_death,
                                yll_pre                   = owsa_sample_p2_ul$yll_pre,
                                yll_post                  = owsa_sample_p2_ul$yll_post,
                                yld_pre                   = owsa_sample_p2_ul$yld_pre,
                                yld_post                  = owsa_sample_p2_ul$yld_post,
                                yld_averted               = owsa_sample_p2_ul$yld_averted,
                                yll_averted               = owsa_sample_p2_ul$yll_averted,
                                daly_total                = owsa_sample_p2_ul$daly_total,
                                icer_daly                 = owsa_sample_p2_ul$icer_daly
  )]
}

# append ul and ll of p1
append_icer_p2 <- data.table :: rbindlist(list(icer_dt_owsa_central, icer_dt_owsa_p1_ll, 
                                               icer_dt_owsa_p1_ul,icer_dt_owsa_p2_ll, 
                                               icer_dt_owsa_p2_ul))

# 
data.central <- barplot(append_icer_p2$icer_daly)

