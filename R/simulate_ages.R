require(dplyr)

n_c = 50 ## number of bins
n_s = 20 ## number of samples to take from each bin  
n_r = 1000 ## number of times to repeat the experiment
Size = 1:n_s
Results_rcz = array(NA, dim=c(n_r,n_c,2))


getP_aa <- function(theta, aa, sigma){
  pal = NA
  if(theta == 1){
    pal <- pnorm(theta, aa, sigma)
  } else if (theta > 1 & theta < n_c){
    pal <- pnorm(theta+1, aa, sigma) - pnorm(theta, aa, sigma)
  } else if (theta == n_c){
    pal <- 1-pnorm(theta, aa, sigma)
  }
  pal
  
}

P_aa <- matrix(NA, nrow = n_c, ncol = n_c)
for(a_obs in 1:n_c){
  P_aa[,a_obs] <- sapply(1:n_c, FUN = getP_aa, aa = a_obs, sigma = rnorm(1,3,0.2))
}

P_aa <- as.matrix(P_aa)
## put this in a loop over n_r
prop_c = exp(rnorm(n_c)) ## randomly generate a vector of lognormal proportions for this experiment; length equal to number of bins 
prop_c = prop_c / sum(prop_c) ## rescale to get "true" population comp in bin (sums to 1)

# sample from multinomial; the prob of observing a given bin is given by prop_c
Y_sc = t(sapply( Size, FUN=rmultinom, n=1, prob=prop_c))

# Convert to data frame
## These are the sampled expected ages a-tilde
Y_iz = data.frame(
  y = as.vector(Y_sc), ## number observed
  c = as.vector(col(Y_sc)), ## bin (e.g. age)
  s = as.vector(row(Y_sc)) ## sample id (akin to "haul", always max n_s)
)

## typical comps setup
expected_1 <- tidyr::pivot_wider(Y_iz, values_from = y, names_from = c, id_cols = s) %>% select(-s) %>%
  as.matrix(.)
## Approach 1 (SS approach)
## Multiply the expected ages, given above by Y_iz, by the ageing error matrix
# approach_1 <- as.matrix(expected_1)  %*% as.matrix(P_aa)
  # 
# sum(expected_1[1,]* P_aa[,1])
caa <- matrix(NA, nrow = n_s, ncol = n_c)
for(s in 1:nrow(expected_1)){
  for(a in 1:ncol(expected_1)){
    caa[s,a] <- sum(P_aa[,a] * expected_1[s,a])/sum(P_aa*n_s)
  }
}
rowSums(caa)
colSums(caa)
