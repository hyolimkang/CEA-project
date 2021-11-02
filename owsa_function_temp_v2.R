# one way sensitivity analysis

# load libraries
library (data.table)
library (ggplot2)

rm (list = ls())

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







