

#Script 4 - 
#Here we work directly with the full posterior trajectory samples (from a fitted model object called ‘out’) 
#to calculate three types of trend metrics, summarize their uncertainty correctly, 
#and run a breakpoint analysis on every posterior trajectory.

#For the publication in review at Proc B 

###############################################################################

# Example workflow for colony-level abundance, regional totals, trends,

# and breakpoint analysis.

###############################################################################




library(tidyverse)
library(segmented)


# `out` — This is the fitted JAGS model object.
# It includes MCMC posterior samples stored in out$sims.list.
# Here, we assume `out` contains posterior samples of colony-level abundance:
#   N[posterior_draw, colony, year]

#save(out, file = "fitted_model_00-24_Nov2025.RData")
setwd("C:/Users/rfo37/OneDrive - University of Canterbury/EP Manuscript/revised versions/Supplementary material/Model 2005-2024") #data")

out <- load(file = "fitted_model_05-24_Nov2025.RData")

out


###############################################################################

# SECTION 1 — EXTRACT COLONY-LEVEL POSTERIOR SAMPLES

###############################################################################

# The JAGS model has returned posterior samples for abundance at each colony
# for each year. These are stored in a 3-D array:

#
#   N_samps[ sample_index , colony_index , year_index ]
#

# Example dimensions typically look like:
#  

N_samps <- out$sims.list$N
dim(N_samps)  #30000 posterior samples × 7 colonies × 20 years



###############################################################################

# SECTION 2 — SUM POSTERIOR SAMPLES TO GET THE REGIONAL TOTAL

###############################################################################

# Each posterior draw gives abundance at all colonies.
# sum across colonies to get "regional abundance" for each draw × year:



# The result is a matrix:
#   N_regional_samps[sample_index , year_index]

N_regional_samps <- apply(N_samps, c(1, 3), sum)
dim(N_regional_samps) # 30000 posterior samples × 20 years




###############################################################################

# SECTION 3 — PLOT EXAMPLE POSTERIOR TRAJECTORIES

###############################################################################

# It is helpful to visualize several posterior trajectories to illustrate
# uncertainty in population trends over time. We select 25 posterior draws
# at random (or in order) and plot them.


N_samples_to_plot <- 25


# Convert the selected draws into a long-format tibble for ggplot

traj_df <- as_tibble(N_regional_samps[1:N_samples_to_plot, ]) %>%
  mutate(sample = 1:N_samples_to_plot) %>%
  pivot_longer(
    cols = -sample,
    names_to = "year_index",
    values_to = "abundance"
  ) %>%
  
  mutate(
    year_index = as.integer(str_remove(year_index, "V")),
    year       = year_index
  )



# Plot the trajectories

p1 <- ggplot(traj_df, aes(x = year, y = abundance, group = sample)) +
  geom_line(color = "black", alpha = 0.5, size = 1) +
  labs(
    x = "Year",
    y = "Regional abundance (posterior samples)",
    title = "Posterior trajectories of regional abundance"
  ) + #scale_y_continuous(label = c("40000", "60000", "80000", "100000"))  +
  theme_bw()

print(p1)
ggsave("posterior_traj_correct.jpg",plot = p1,
       scale = 1, width = 7, height = 5, units = c("in"), bg = "white", dpi = 600)




###############################################################################

# SECTION 4 — CALCULATE POSTERIOR TRENDS

###############################################################################

# The power of Bayesian analysis is that every derived quantity can be computed
# *for every posterior sample*, producing a posterior distribution for:
#   - total change
#   - percent change
#   - average annual change (geometric mean)
#
# Below we compute three different trend measures. These represent uncertainty
# naturally because each is computed for all 30,000 posterior draws.



### 4a. Total change from first to last year
N_change <- N_regional_samps[, 20] - N_regional_samps[, 1]
length(N_change)           # number of posterior samples
mean(N_change)# > 0)         # posterior probability of an increase
#0.1015333
#-17125.99 - using new 'out' file... hmm
quantile(N_change, c(0.025, 0.5, 0.95, 0.975))  # median + 95% credible interval
#      2.5%        50%        95%      97.5% 
#  -41799.265 -18155.231   6668.248  13174.258 



### 4b. Percent change between first and last year
perc_change <- 100 * (N_regional_samps[, 20] - N_regional_samps[, 1]) /
  N_regional_samps[, 1]

mean(perc_change)# > 0) 
#-22.58727

quantile(perc_change, c(0.025, 0.5, 0.95, 0.975))
#       2.5%        50%        95%      97.5% 
#  -49.256600 -25.189082   9.745233  19.591561 



### 4c. Geometric mean annual percent change (over the 20 year period)
trend_annual <- 100 * ((N_regional_samps[, 20] /
                        N_regional_samps[, 1])^(1 / (20 - 1)) - 1)
mean(trend_annual)# > 0) 
#-1.457762
#0.1015333
quantile(trend_annual, c(0.025, 0.5, 0.95, 0.975))
#     2.5%        50%        95%      97.5% 
#  -3.5074772 -1.5157962  0.4906283  0.9460900 


###############################################################################

# SECTION 5 — BREAKPOINT ANALYSIS ON POSTERIOR TRAJECTORIES

###############################################################################

# A breakpoint model tries to find a year at which the population trend
# (when plotted on a log scale) changes slope.
#

# Why log scale?
#   logN transforms exponential or multiplicative trends into linear slopes,
#   making breakpoints easier to interpret.
#

# We use the `segmented` package, which fits piecewise linear regressions.
# For each posterior sample, we:
#   1. Fit a regression log(N_t) ~ year
#   2. Ask segmented() to estimate a single breakpoint
#   3. Store the breakpoint year
#

# Some fits may fail, so we use tryCatch to avoid stopping the script.


years <- 1:ncol(N_regional_samps)
T     <- length(years)
S     <- nrow(N_regional_samps)


### 5a. Example: fit breakpoint model to ONE posterior trajectory
df <- data.frame(
  year = years,
  logN = log(N_regional_samps[1, ])
)


# psi = median(years) gives midpoint as starting value
seg_example <- segmented(
  lm(logN ~ year, data = df),
  seg.Z = ~year,
  psi   = median(years)
)

seg_example$psi   # estimated breakpoint location (example)



### 5b. Fit breakpoint model to *all* posterior trajectories
breakpoints <- rep(NA, S)

# Note: A for-loop is shown for clarity, but apply()/purrr or parallelization would be faster

for (sample_number in 1:S) {
  df <- data.frame(
    year = years,
    logN = log(N_regional_samps[sample_number, ])
  )

  # Attempt segmented fit; if it fails, return NA instead of error
  seg_fit <- tryCatch(
    segmented(lm(logN ~ year, data = df),
              seg.Z = ~year,
              psi   = median(years)),
    error = function(e) return(NA)
  )

  # Extract breakpoint if model fit succeeded
  if (is.list(seg_fit)) {
    breakpoints[sample_number] <- seg_fit$psi[2]
    }
  
  print(sample_number)  # progress indicator
  }


# Summaries of breakpoint estimates
dev.off()
mean(!is.na(breakpoints))               # % of successful fits
#1
hist(breakpoints, main = "Breakpoint distribution",
     xlab = "Estimated breakpoint year", xlim = c(0,20))
median(breakpoints, na.rm = TRUE)
#15.89599 - round up to 16
quantile(breakpoints, c(0.025, 0.975), na.rm = TRUE)
#     2.5%     97.5% 
#  4.349016 17.895890 
summary(breakpoints)
mean(breakpoints>14) #probability of break points occurred after 2018
mean(breakpoints>16) #probability of break points occured after 2020

#median breakpoint = year 16 corresponding to 2020...
Ross$year

#95% CI: 2008 to 2022...
# 0.7816 probability break point occurred after 2018 


### 6. trends on each sector before and after break point  -----
#Ok now can I run the trend analysis again on each sector pre and post 2020 


### 6a. Total change from first to break point (2020) #last year
N_change05_20 <- N_regional_samps[, 16] - N_regional_samps[, 1]
length(N_change05_20)           # number of posterior samples
mean(N_change05_20)# > 0)         # posterior probability of an increase
#6183.701
quantile(N_change05_20, c(0.025, 0.5, 0.95, 0.975))  # median + 95% credible interval
#      2.5%        50%        95%      97.5% 
#  -17328.658   5993.852  26655.416  31060.349 



### 6b. Percent change between first and break point (2020) #last year
perc_change05_20 <- 100 * (N_regional_samps[, 16] - N_regional_samps[, 1]) /
  N_regional_samps[, 1]

mean(perc_change05_20)# > 0) 
##9.491056
quantile(perc_change05_20, c(0.025, 0.5, 0.95, 0.975))
#     2.5%        50%        95%      97.5% 
#  -20.572029   8.263758  39.902920  47.106373 



### 6c. Geometric mean annual percent change (before breakpoint in 2020)
trend_annual05_20 <- 100 * ((N_regional_samps[, 16] /
                          N_regional_samps[, 1])^(1 / (16 - 1)) - 1)
mean(trend_annual05_20)# > 0) 
#0.5296989
quantile(trend_annual05_20, c(0.025, 0.5, 0.95, 0.975))
#      2.5%        50%        95%      97.5% 
#  -1.5237358  0.5307386  2.2637667  2.6066320 


#And now we run the same for the last years after teh break point - change the other value 

### 6d. Total change from first to break point (2020) #last year
N_change20_24 <- N_regional_samps[, 20] - N_regional_samps[, 16]
length(N_change20_24)           # number of posterior samples
mean(N_change20_24)# > 0)         # posterior probability of an increase
#-23399.69
quantile(N_change20_24, c(0.025, 0.5, 0.95, 0.975))  # median + 95% credible interval
#      2.5%        50%        95%      97.5% 
#  -46973.153 -23744.921  -2772.712   2150.727 




### 6e. Percent change between first and break point (2020) #last year
perc_change20_24  <- 100 * (N_regional_samps[, 20] - N_regional_samps[, 16]) /
  N_regional_samps[, 1]

mean(perc_change20_24 )# > 0) 
##-32.07832
quantile(perc_change20_24 , c(0.025, 0.5, 0.95, 0.975))
#      2.5%        50%        95%      97.5% 
#  -65.533776 -32.210760  -3.803805   2.899079 



### 6f. Geometric mean annual percent change (over the 20 year period)
trend_annual20_24  <- 100 * ((N_regional_samps[, 20] /
                               N_regional_samps[, 16])^(1 / ( 5 - 1)) - 1)
mean(trend_annual20_24 )# > 0) 
#-8.476554

quantile(trend_annual20_24 , c(0.025, 0.5, 0.95, 0.975))
#2.5%         50%         95%       97.5% 
#  -16.2102358  -8.7541116  -0.9312702   0.7008894



###between 2021 and 2022 ------------
### 7d. Total change from first to break point (2020) #last year
N_change21_22 <- N_regional_samps[, 18] - N_regional_samps[, 17]
length(N_change21_22)           # number of posterior samples
mean(N_change21_22)# > 0)         # posterior probability of an increase
#-11517.9
quantile(N_change21_22, c(0.025, 0.5, 0.95, 0.975))  # median + 95% credible interval
#     2.5%         50%         95%       97.5% 
#  -27208.3390 -11242.9689    164.3906   2470.6104 


### 7e. Percent change between first and break point (2020) #last year
perc_change21_22  <- 100 * (N_regional_samps[, 18] - N_regional_samps[, 17]) /
  N_regional_samps[, 1]

mean(perc_change21_22 )# > 0) 
## -15.81292
quantile(perc_change21_22 , c(0.025, 0.5, 0.95, 0.975))
#     2.5%        50%        95%      97.5% 
#  -38.069969 -15.248281   0.230751   3.371046 


### 7f. Geometric mean annual percent change (over the 20 year period)
trend_annual21_22  <- 100 * ((N_regional_samps[, 18] /
                                N_regional_samps[, 17])^(1 / ( 2 - 1)) - 1)
mean(trend_annual21_22 )# > 0) 
#-8.476554

quantile(trend_annual21_22 , c(0.025, 0.5, 0.95, 0.975))
#2.5%         50%         95%       97.5% 
#  -16.2102358  -8.7541116  -0.9312702   0.7008894




