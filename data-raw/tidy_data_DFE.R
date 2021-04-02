library(readr)
library(dplyr)
library(usethis)
library(here)

# What : - Filters data to keep only common genotypes across treatments
#        - Transforms fitnesses in selection coefficients (s)
#        - Summarises : - mean of s over replicates per EnvironmentXPositionXAA
#                       - approximated standard deviation of the measurement error over the replicate per EnvironmentXPositionXAA
#        - Saves summarized data in a format for lazy loading.

# Reads the raw data
raw_data <- read_csv(here("data-raw", "raw_data_DFE.csv"))

# Finds the common combination of Position and AA across all the Treatment
common_Position_AA <- raw_data %>%
  group_by(Position, AA) %>%
  tally() %>%
  filter(n == n_distinct(raw_data$Treatment))

# Extracts only the rows with common Position and AA across all the Treatment
data_filtered <- semi_join(raw_data, common_Position_AA, by = c("Position", "AA"))

# Compute the mean of replicate and the variance of the error for each combination Position-AA in each Environment
average_sd <- function(x) {sqrt(sum(x^2)) / length(x)} # function to compute the standard deviation of the mean.
data_DFE <- data_filtered %>%
  mutate(sd = abs(Upper_CI - Low_CI) / 4) %>% #approximation of the standard deviation by the standard error.
  filter(AA != "WILDTYPE") %>% # remove the wildtype from the data
  # group_by(Treatment) %>%
  # arrange(AA != "WILDTYPE", AA) %>% #position wildtype as the first row of each group
  # mutate(s = s - first(s)) %>% # shift all fitnesses by the wt fitness -> becomes selection coefficients
  # ungroup() %>%
  group_by(Environment, Position, AA, Expression_level) %>%
  mutate(mean_s = mean(s), # mean of selection coefficients over the replicates per environmentXpositionXAA
         mean_sd = average_sd(sd)) # standard deviation of the measurement error over the replicates per environmentXpositionXAA

write_tsv(data_DFE, here("data-raw", "tidy_data.tsv"))
use_data(data_DFE, overwrite = TRUE)
