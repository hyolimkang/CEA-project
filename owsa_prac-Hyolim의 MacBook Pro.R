# input parameters datatable
param      <- data.frame   (v_cost          = c(mid = 2.96,  low = 2.96,   high = 2.96), 
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
param_table <- data.frame(t(param))
names(param_table) <- c("mid", "low", "high")
# input parameters datatable
param_table$name <-   c("v_cost", "delivery_cost", "tot_pop", "vacc_pop", "postvacc_p1",
                      "facility_cost", "dmc", "dnmc", "indirect", "ve", "dw", "illness_duration",
                      "cfr", "age_death", "life_exp")

# icer datatable
icer <- data.frame(matrix(NA, nrow=nrow(param_table), ncol = 3))
names(icer) <- c("mid", "low", "high")

# output calculating function 
  owsa_icer <- function(param) {
  
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
  
  return (list (icer_daly = icer_daly)) 
}


# for (i in 1:nrow(param_table)){
for (nm in param_table$name){ 
  # mid <- param_table[i, "mid"]
  # low <- param_table[i, "low"]
  # high <- param_table[i, "high"]
  mid <- param_table[param_table$name == nm, c("name","mid")]
  names(mid) <- c("name", "val")
  low <- param_table[param_table$name == nm, c("name","low")]
  names(low) <- c("name", "val")
  high <- param_table[param_table$name == nm, c("name","high")]
  names(high) <- c("name", "val")
  
  rest <- param_table[param_table$name != nm, c("name","mid")]
  names(rest) <- c("name", "val")
  
  newparam <- rbind(mid, rest)
  
  icer[, "mid"] <- calc_icer(newparam)
  newparam <- rbind(low, rest)
  icer[, "low"] <- calc_icer(newparam)
  newparam <- rbind(high, rest)
  icer[, "high"] <- calc_icer(newparam)
}


