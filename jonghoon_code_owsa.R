params <- data.

    v_cost           = c(mid = 2.96,  low = 2.96,   high = 2.96), 
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

    for (nm in icer$name) {
      # nm = "v_cost"
      temp_input <- base_input
      low <- param_table[param_table$name == nm, "low"]
      temp_input[temp_input$name == nm, "val"] <- low
      low_icer <- icer_calc(temp_input)
      icer[icer$name == nm, "low"] <- low_icer$icer_daly
      
      high <- param_table[param_table$name == nm, "high"]
      temp_input[temp_input$name == nm, "val"] <- high
      high_icer <- icer_calc(temp_input)
      icer[icer$name == nm, "high"] <- high_icer$icer_daly
    }
    
    
    
    low_icer <- icer_calc(temp_input)
    
    for (i in 1:3) {
      mid <- param_table[param_table$name == nm, c("name","mid")]
      names(mid) <- c("name", "val")
      
    }