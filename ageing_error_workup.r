library(ggplot2)
library(tidyr)

# Define the number of age bins
n_bins <- 25
n_samples <- n_bins*10
n_experiments <- 100
# Initialize the ageing-error matrix
ageing_error_matrix <- matrix(0, nrow = n_bins, ncol = n_bins)

# Fill the ageing-error matrix
for (i in 1:n_bins) {
  # Generate a normally-distributed probability vector
  prob_vector <- dnorm(1:n_bins, mean = i, sd = 5)
  
  # Normalize the probability vector so it sums to 1
  prob_vector <- prob_vector / sum(prob_vector)

  # Assign the probability vector to the corresponding column in the ageing-error matrix
  ageing_error_matrix[, i] <- prob_vector
}

# Print the ageing-error matrix
#print(ageing_error_matrix)

# Plot a column of the ageing-error matrix
#plot(ageing_error_matrix[,10])

## put this in a loop over n_r
prop_c = exp(rnorm(n_bins)) ## randomly generate a vector of lognormal proportions for this experiment; length equal to number of bins 
prop_c = prop_c / sum(prop_c) ## rescale to get "true" population comp in bin (sums to 1)

# sample from multinomial; the prob of observing a given bin is given by prop_c
experiment_results <- t(rmultinom(n = n_experiments, size = n_samples, prob = prop_c)) 

## normalize as computation happens
routine1_results <- (experiment_results/rowSums(experiment_results)) %*% ageing_error_matrix

# Initialize the results matrix
routine2_results <- matrix(0, ncol = n_bins, nrow =n_experiments)
# Perform the second routine for each experiment
for (i in 1:n_experiments) {
    ## make a placeholder with a row for each age candidate
    temp <- matrix(0, nrow = n_bins, ncol = n_bins)
    for (j in 1:n_bins) {
    # Sample from the observed ages using the ageing-error matrix as the probabilities
    temp[j,] <- rmultinom(n=1, ## once per experiment
    size = experiment_results[i, j],   ## as if reading this many otoliths
    prob = ageing_error_matrix[, j])   ## the probability of observing a given bin is given by the ageing error matrix
  }
   routine2_results[i, ] <- colSums(temp)
}
#routine2_results[routine2_results==0] <- 1e-5 ## avoid division by zero
routine2_results <- routine2_results/(rowSums(routine2_results)) ## normalize

# Convert the matrices to data frames
routine1_df <- as.data.frame(routine1_results)
routine2_df <- as.data.frame(routine2_results)

# Add a column for the experiment number
routine1_df$experiment <- 1:nrow(routine1_df)
routine2_df$experiment <- 1:nrow(routine2_df)

# Reshape the data frames to a long format
routine1_long <- reshape2::melt(routine1_df, id.vars = "experiment", variable.name = "age_bin", value.name = "count")
routine2_long <- reshape2::melt(routine2_df, id.vars = "experiment", variable.name = "age_bin", value.name = "count")

# Add a column for the routine number
routine1_long$routine <- "Routine 1"
routine2_long$routine <- "Routine 2"

# Combine the data frames
combined_df <- rbind(routine1_long, routine2_long) 

# Convert the age_bin column to numeric
combined_df$age_bin <- as.numeric(gsub("V", "", combined_df$age_bin))

## average across experiments
#combined_df <- combined_df %>%
#group_by(age_bin, routine) %>%
#summarise(mean_count = mean(count)) #%>%
  #group_modify(~as.data.frame(t(quantile(.$count))))

#names(combined_df)[3:7] <- c("Q_0", "Q_25", "Q_50", "Q_75", "Q_100")


# Create the density plot
ggplot(combined_df, aes(x = age_bin, y  = count, color = routine, fill = routine, group = experiment)) +
    geom_line(data = data.frame(age_bin = 1:n_bins, count = prop_c,
     routine = 'true population', experiment = 1))+
    geom_line() +
  facet_wrap(~routine)+
 #geom_ribbon(aes(ymin = Q_25, ymax =Q_75), alpha = 0.3) +
  theme_minimal() +
  scale_fill_manual(values = c("Routine 1" = "red", "Routine 2" = "blue", "true population" = "black")) +
  scale_color_manual(values = c("Routine 1" = "red", "Routine 2" = "blue", "true population" = "black")) +
  labs(x = "Age Bin", y = "Frequency", color = "",
  title = "Comparison of Age Composition Sampling Routines") 


