# i need to put v_cost_mid, v_cost_low, v_cost_high... FOR 15 parameters ->45 inputs needed
# or if i put only put the base value here, i will need to tell the function to take the 
# low and high value when all other parameters take the base values
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
                  life_exp,
                  parameters_list) {
 
  for(parameter in parameters_list){
    # if i want to get the low values and control all the other variables into mid
    if(parameter ) {
  
# mid values
totcost_vacc   <-     (v_cost + delivery_cost) * (vacc_pop) + 
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

# low values


return(list(icer_daly_low  = icer_daly_low,
            icer_daly_mid  = icer_daly_mid,
            icer_daly_high = icer_daly_high))
       

}

# results: i want the lower value of icer when controlling all the other variables
# for instance, if i wan the lower value icer of v_cost, i want to put lower value for v_cost 
# and controlling all other variables to the mid values

    parameters_list <-   c (v_cost         = c(mid = 2.96,  low = 2.96,   high = 2.96), 
                          delivery_cost    = c(mid = 1.49,  low = 0.246,  high = 4.752),
                          tot_pop          = c(mid = 159831,low = 159831, high = 159831),
                          vacc_pop         = c(mid =113420, low = 113420, high = 113420),
                          postvacc_p1      = c(mid = 23,    low = 23,      high = 23),
                          facility_cost    = c(mid = 97.33, low = 16.676, high = 301.209),
                          dmc              = c(mid = 183.07,low = 18.507, high =532.241),
                          dnmc             = c(mid =30.13,  low =  4.361, high = 118.718), 
                          indirect         = c(mid = 73.09, low = 32.370, high = 135.262),
                          ve               = c(mid = 0.81,  low =  0.55,  high =  0.916), 
                          dw               = c(mid = 0.052, low =  0.080, high = 0.128),
                          illness_duration = c(mid =0.043,  low = 0.0392, high = 0.0469),
                          cfr              = c(mid = 0.019, low = 0.011,  high =  0.0441),
                          age_death        = c(mid = 7.5,   low = 1.840,  high = 13.95),
                          life_exp         = c(mid = 72.68, low =  69.576,high = 78.70),
                          parameters_change, 
                          change)
    
icer_low   <- owsa (v_cost            = 2.96,
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
                    life_exp          =  72.68)


