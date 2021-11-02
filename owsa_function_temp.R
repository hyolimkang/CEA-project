rm (list = ls())

# icer calculating function

  icer_calc <- function (params){
  
  # select a value in a second column in param data table if param$name is "name"
  
  v_cost <- params[params$name == "v_cost", 2]
  delivery_cost <- params[params$name == "delivery_cost", 2]
  tot_pop <- params[params$name == "tot_pop", 2]
  vacc_pop <- params[params$name == "vacc_pop", 2]
  postvacc_p1  <- params[params$name == "postvacc_p1", 2]
  facility_cost <- params[params$name == "facility_cost", 2]
  dmc <- params[params$name == "dmc", 2]
  dnmc <- params[params$name == "dnmc", 2]
  indirect <- params[params$name == "indirect", 2]
  ve <- params[params$name == "ve", 2]
  dw <- params[params$name == "dw", 2]
  illness_duration <- params[params$name == "illness_duration", 2]
  cfr <- params[params$name == "cfr", 2]
  age_death <- params[params$name == "age_death", 2]
  life_exp <- params[params$name == "life_exp", 2]
  
  
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


return(list(icer_daly = icer_daly))


}


# parameter data table
param       <- data.frame (v_cost           = c(mid = 2.96,  low = 2.96,   high = 2.96), 
                           delivery_cost    = c(mid = 1.49,  low = 0.246,  high = 4.752),
                           tot_pop          = c(mid = 159831,low = 159831, high = 159831),
                           vacc_pop         = c(mid = 113420, low = 113420, high = 113420),
                           postvacc_p1      = c(mid = 23,    low =23,      high =23),
                           facility_cost    = c(mid = 97.33, low = 16.676, high = 301.209),
                           dmc              = c(mid = 183.07,low = 18.507, high =532.241),
                           dnmc             = c(mid =30.13,  low =  4.361, high = 118.718), 
                           indirect         = c(mid = 73.09, low = 32.370, high = 135.262),
                           ve               = c(mid = 0.81,  low =  0.55,  high =  0.916), 
                           dw               = c(mid = 0.052, low =  0.032, high = 0.078),
                           illness_duration = c(mid =0.043,  low = 0.0392, high = 0.0469),
                           cfr              = c(mid = 0.019, low = 0.011,  high =  0.0441),
                           age_death        = c(mid = 7.5,   low = 1.840,  high = 13.95),
                           life_exp         = c(mid = 72.68, low =  69.576,high = 78.70))
# transpose the table
param_table <- data.frame(t(param))
# give name to each columns 
names(param_table) <- c("mid", "low", "high")
param_table$name <-   c("v_cost", "delivery_cost", "tot_pop", "vacc_pop", "postvacc_p1", 
                         "facility_cost", "dmc", "dnmc", "indirect", "ve", "dw", "illness_duration",
                         "cfr", "age_death", "life_exp")


# i need to put v_cost_mid, v_cost_low, v_cost_high... FOR 15 parameters ->45 inputs needed
# or if i put only put the base value here, i will need to tell the function to take the 
# low and high value when all other parameters take the base values

#empty icer table
icer <- data.frame(matrix(NA, nrow = nrow(param_table), ncol = 4))
names(icer) <- c("name","mid", "low","high")
#bring names in the icer_name column from the param_table name column
icer$name<- param_table$name


# ref: https://www.infoworld.com/article/3530348/r-datatable-symbols-and-operators-you-should-know.html (datatable with quotes)
# datatble[i,j] => [row, column]

#central icer value
icer_mid <- icer_calc(param_table[,c("name","mid")])

base_input <- param_table[,c("name","mid")]
names(base_input) <- c("name","val")






