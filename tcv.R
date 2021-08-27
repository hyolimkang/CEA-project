# TCV cost-effectiveness

# load libraries
library (dampack)
library (data.table)
library (ggplot2)
library (tictoc)

# clear workspace
rm (list = ls())

# ------------------------------------------------------------------------------
# function
# create psa samples for input parameters
# ------------------------------------------------------------------------------
create_psa_sample <- function (sample_n) {

  # parameters
  my_params <-c ("v_cost",          # vaccine cost
                 "facility_cost",   # health facility cost
                 "dmc",             # direct medical cost
                 "dnmc",            # direct non-medical cost
                 "indirect",        # indirect cost
                 "ve")              # vaccine-effectiveness

  # constant value
  # "postvacc"        # post-vaccination case

  # trickle down uncertainty
  # "prevacc",         # pre-vaccination case
  # "totcost_vacc",    # total cost of vaccinated group
  # "totcost_unvacc")  # total cost of unvaccinated group


  # distributions for parameters
  my_dists <- c("log-normal",
                "log-normal",
                "log-normal",
                "log-normal",
                "log-normal",
                "truncated-normal")


  # parameter types
  my_parameterization_types <-c ("mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd",
                                 "mean, sd, ll, ul")


  # values for parameters
  my_dists_params <- list( c (4.45,   0.22673),
                           c (58.64,  0.1136),
                           c (183.07, 229.31),
                           c (30.13,  31.06),
                           c (110.95, 26.58),
                           c (.87,     0.02666667, .79, .95))

  # c(384, 0.487),
  # c(50, 0.5),
  # c(1044187.37, 0.23),
  # c(147226.92, 0.59))


  # input samples for PSA
  my_psa <- gen_psa_samp (params                 = my_params,
                          dists                  = my_dists,
                          parameterization_types = my_parameterization_types,
                          dists_params           = my_dists_params,
                          n                      = sample_n)

  return (my_psa)

}
# end of function -- create_psa_sample
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# function - computer icer for a given input parameter sample
compute_icer <- function (v_cost,          # vaccine cost
                          facility_cost,   # health facility cost
                          dmc,
                          dnmc,
                          indirect,
                          ve) {


  # compute incremental_effectiveness
  incremental_effectiveness <- runif (1, 0, 10)

  # incremental_cost
  incremental_cost <- runif (1, 0, 50)


  # estimate icer
  icer <- incremental_cost / incremental_effectiveness


  return (list (incremental_cost          = incremental_cost,
                incremental_effectiveness = incremental_effectiveness,
                icer                      = icer))
}
# end of function -- compute_icer
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# create ICER plot
create_icer_plot <- function (dt,
                              plot_file) {

  # open plot file
  pdf (plot_file)

  ## @knitr custom-ceplane-plot
  xlim <- max (icer_dt$incremental_effectiveness)
  ylim <- max (icer_dt$incremental_cost)

  # create plot
  p <- ggplot (dt,
         aes(x = incremental_effectiveness, y = incremental_cost)) +
    geom_jitter(size = .5)  +
    xlab ("DALYs averted") +
    ylab ("Incremental cost") +
    scale_y_continuous (limit = c(-xlim, xlim)) +
    scale_x_continuous (limits = c(-xlim, xlim)) +
    theme (legend.position = "bottom") +
    scale_colour_discrete (name = "TCV vaccination strategy") +
    geom_abline (slope = 5, linetype = "dashed") +
    geom_hline (yintercept = 0) +
    geom_vline (xintercept = 0)

  # print plot
  print (p)

  # close plot file
  dev.off ()
}
# end of function -- create_icer_plot
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------
# main program
# ------------------------------------------------------------------------------

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
icer_dt [, c("incremental_cost",
             "incremental_effectiveness",
             "icer")] <- 0

# loop through psa sample to generate icers
for (i in 1:sample_size) {

  icer_sample <- compute_icer (psa_sample$v_cost        [i],   # vaccine cost
                               psa_sample$facility_cost [i],   # health facility cost
                               psa_sample$dmc           [i],
                               psa_sample$dnmc          [i],
                               psa_sample$indirect      [i],
                               psa_sample$ve            [i])

  icer_dt [i, `:=` (incremental_cost          = icer_sample$incremental_cost,
                    incremental_effectiveness = icer_sample$incremental_effectiveness,
                    icer                      = icer_sample$icer)]
}

# estimate icer -- median and 95% credible intervals
icer_estimate <- quantile (x = icer_dt$icer,
                           probs = c (0.5, 0.025, 0.975))

# print icer -- median and 95% credible intervals
print ("icer (median and 95% credible interval)")
print (icer_estimate)

# create ICER plot
create_icer_plot (dt        = icer_dt,
                  plot_file = "plot.pdf")


# create plot
# http://r-statistics.co/ggplot2-Tutorial-With-R.html

# test sample visualization
# hist (psa_sample$ve)



# stop time
toc ()
print (Sys.time ())

# end of -- main program
# ------------------------------------------------------------------------------







# ------------------------------------------------------------------------------
# things to do:
# - use same year as reference year for all cost calculations
# - CHEERS checklist: https://www.equator-network.org/wp-content/uploads/2013/04/Revised-CHEERS-Checklist-Oct13.pdf

# learning resource
# https://github.com/hesim-dev/rcea
# http://r-statistics.co/ggplot2-Tutorial-With-R.html
# ------------------------------------------------------------------------------

