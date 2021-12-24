owsa_new   <-function(v_cost        = c(mid = 2.96,  low = 2.96,   high = 2.96), 
                   delivery_cost    = c(mid = 1.49,  low = 0.246,  high = 4.752),
                   tot_pop          = c(mid = 159831,low = 159831, high = 159831),
                   vacc_pop         = c(mid =113420, low = 113420, high = 113420),
                   postvacc_p1      = c(mid = 23,    low = 23,     high = 23),
                   facility_cost    = c(mid = 97.33, low = 16.676, high = 301.209),
                   dmc              = c(mid = 183.07,low = 18.507, high = 532.241),
                   dnmc             = c(mid =30.13,  low =  4.361, high = 118.718), 
                   indirect         = c(mid = 73.09, low = 32.370, high = 135.262),
                   ve               = c(mid = 0.81,  low =  0.55,  high =  0.916), 
                   dw               = c(mid = 0.052, low =  0.080, high = 0.128),
                   illness_duration = c(mid =0.043,  low = 0.0392, high = 0.0469),
                   cfr              = c(mid = 0.019, low = 0.011,  high =  0.0441),
                   age_death        = c(mid = 7.5,   low = 1.840,  high = 13.95),
                   life_exp         = c(mid = 72.68, low =  69.576,high = 78.70)) 
  {
  
  # values
  
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

parameters <- c   (v_cost           = c(mid = 2.96,  low = 2.96,   high = 2.96), 
                   delivery_cost    = c(mid = 1.49,  low = 0.246,  high = 4.752),
                   tot_pop          = c(mid = 159831,low = 159831, high = 159831),
                   vacc_pop         = c(mid =113420, low = 113420, high = 113420),
                   postvacc_p1      = c(mid = 23,    low = 23,     high = 23),
                   facility_cost    = c(mid = 97.33, low = 16.676, high = 301.209),
                   dmc              = c(mid = 183.07,low = 18.507, high = 532.241),
                   dnmc             = c(mid =30.13,  low =  4.361, high = 118.718), 
                   indirect         = c(mid = 73.09, low = 32.370, high = 135.262),
                   ve               = c(mid = 0.81,  low =  0.55,  high =  0.916), 
                   dw               = c(mid = 0.052, low =  0.080, high = 0.128),
                   illness_duration = c(mid =0.043,  low = 0.0392, high = 0.0469),
                   cfr              = c(mid = 0.019, low = 0.011,  high =  0.0441),
                   age_death        = c(mid = 7.5,   low = 1.840,  high = 13.95),
                   life_exp         = c(mid = 72.68, low =  69.576,high = 78.70))  

param_low   <-  parameters[c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,44)]

icer_dt_owsa_ul <- data.table (run_id = 1 : length(parameters))

param_low_samp   <- owsa_new(parameters)


icer_dt_owsa_ul [i, `:=`  ( totcost_unvacc            = param_low_samp$totcost_unvacc,
                           totcost_vacc              = param_low_samp$totcost_vacc,
                           incremental_cost          = param_low_samp$incremental_cost,
                           incremental_effect        = param_low_samp$incremental_effect,
                           icer                      = param_low_samp$icer,
                           prevacc                   = param_low_samp$prevacc,
                           temp_distance             = param_low_samp$temp_distance,
                           pre_death                 = param_low_samp$pre_death,
                           post_death                = param_low_samp$post_death,
                           yll_pre                   = param_low_samp$yll_pre,
                           yll_post                  = param_low_samp$yll_post,
                           yld_pre                   = param_low_samp$yld_pre,
                           yld_post                  = param_low_samp$yld_post,
                           yld_averted               = param_low_samp$yld_averted,
                           yll_averted               = param_low_samp$yll_averted,
                           daly_total                = param_low_samp$daly_total,
                           icer_daly                 = param_low_samp$icer_daly
)]


