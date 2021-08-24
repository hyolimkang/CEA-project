# cost-effectiveness

# load library 
library ("dampack")

# clear workspace
rm (list = ls())

# parameters
my_params <-c ("v_cost",          # vaccine cost 
               "facility_cost",   # health facility cost
               "dmc", 
               "dnmc", 
               "indirect", 
               "ve", 
               "prevacc", 
               "postvacc", 
               "totcost_vacc", 
               "totcost_unvacc", 
               "icer")

# distributions for parameters
my_dists <- c("log-normal", 
              "log-normal", 
              "log-normal", 
              "log-normal",
              "log-normal", 
              "beta", 
              "beta", 
              "beta", 
              "log-normal", 
              "log-normal", 
              "log-normal")

# paramter types
my_parameterization_types <-c ("mean, sd", 
                               "mean, sd", 
                               "mean, sd",
                               "mean, sd",
                               "mean, sd", 
                               "a, b", 
                               "a, b", 
                               "a, b", 
                               "mean, sd", 
                               "mean, sd", 
                               "mean, sd")

# values for parameters
my_dists_params <- list( c (4.45,   0.22673), 
                         c (58.64,  0.1136), 
                         c (183.07, 229.31), 
                         c (30.13,  31.06), 
                         c (110.95, 26.58), 
                         c (87,     0.042), 
                        c(384, 0.487), 
                        c(50, 0.5), 
                        c(1044187.37, 0.23), 
                        c(147226.92, 0.59), 
                        c(2680.57, 0.285))

# input samples for PSA
my_psa <- gen_psa_samp (params                 = my_params, 
                        dists                  = my_dists, 
                        parameterization_types = my_parameterization_types, 
                        dists_params           = my_dists_params, 
                        n                      = 100)

# test sample visualization
hist (my_psa$v_cost)

