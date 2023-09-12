## Annotated/Modified by M Kapur, Oct 2021

## Thorson's Original Version ----
# Based on:
#  https://stats.stackexchange.com/questions/24705/can-i-use-glm-algorithms-to-do-a-multinomial-logistic-regression

# Developed from:
#  Conversations with Brandon Chasco, Dave Warton, Joanna Mills-Flemming

# Please do not re-use without discussions with me (Jim Thorson)
#  NOTE that I would be happy to collaborate on work along these lines

n_c = 10 ## number of bins 
n_s = 20 ## number of samples to take from each bin (sort of like 'years', but there is not carrythru)
n_r = 1000 ## number of times to repeat the experiment

Size = 1:n_s
Results_rcz = array(NA, dim=c(n_r,n_c,2))

for( rI in 1:n_r ){
  prop_c = exp(rnorm(n_c)) ## 10 lognormal proportions
  prop_c = prop_c / sum(prop_c) ## rescale to get "true" population comp in bin (sums to 1)
  
  # sample from multinomial
  ## the prob of observing a given bin is given by prop_c
  Y_sc = t(sapply( Size, FUN=rmultinom, n=1, prob=prop_c))
  
  # Convert to data frame
  Y_iz = data.frame(
    y = as.vector(Y_sc), ## number observed
    c = as.vector(col(Y_sc)), ## bin (e.g. age)
    s = as.vector(row(Y_sc)) ## sample id (akin to "haul", always 20)
  )
  
  # Fit Poisson-GLM
  ## independent effects of haul and bin; based on your observations
  ## can be used to predict more obs...
  G = glm( y ~ 0 + factor(c) + factor(s), data=Y_iz, family=poisson() )
  
  # Extract effects (with data, try to estimate true pop comp (prop_c))
  prophat_c = summary(G)$coef[paste0("factor(c)",1:n_c),'Estimate']
  prophat_c = exp(prophat_c) / sum(exp(prophat_c)) ## rescale, again
  
  # Record results
  Results_rcz[rI,,1] = prop_c ## true observed props
  Results_rcz[rI,,2] = prophat_c ## simulated and or estimated props
}

mn = Results_rcz; rm(Results_rcz)

## shows that we get some error around observations
plot( x=mn[,,1], y=mn[,,2], log="",
      xlab = 'True',ylab = 'Simulated')
abline(a=0, b=1, col="blue")

## plot it in a clearer way
# Results_rcz has nrow = experiment, ncol = bins, arrays = true, obs
require(ggplot2); require(dplyr); require(reshape2)
plotme <-  data.frame(mn[,,1]) %>%
  mutate(haulid = 1:1000) %>%
  melt(id = 'haulid') %>%
  mutate(src = 'obs', 
         age = as.numeric(substr(variable,2,2))) %>%
  rbind( data.frame(mn[,,2]) %>%
           mutate(haulid = 1:1000) %>%
           melt(id = 'haulid') %>%
           mutate(src = 'pred', 
                  age = as.numeric(substr(variable,2,2)))) %>%
  select(-variable) %>%
  filter(haulid < 21) ## for plotting

mn2 <- ggplot(NULL, aes(x = age, y = value, fill = src, color = src)) +
  geom_histogram(data = subset(plotme, src == 'obs', fill = 'obs'),stat = 'identity' ) +
  geom_line(data = subset(plotme, src == 'pred', color = 'pred'), lwd = 1.1) +
  scale_fill_manual(values = c('grey55',NA), labels = c('obs','simulated'))+
  scale_color_manual(values = c(NA, 'blue'), labels = c('obs','simulated'))+
  labs(x = 'age or length', y = 'freq', fill = '',
       color='', title = 'Multinomial Simulator, 20 hauls')+
  theme_minimal()+
  facet_wrap(~haulid)


## Maia's Update: Use Dirichlet Multinomial ----

# help from: https://blog.byronjsmith.com/dirichlet-multinomial-example.html
# also check out https://towardsdatascience.com/adjust-for-overdispersion-in-poisson-regression-4b1f52baa2f1

## this package has the Dirichlet distribution, used to generate sims
# install.packages('dirmult')
require(dirmult)

## for this to work we need two values, the expected fraction (like prop_c) above,
## and  the "concentration" or what we use to control the overdispersion of our data.
## in Thorson's paper this is theta. Larger theta values yield results that look more like
## the multinomial.

theta = 0.05 ## not logged,  normally ranges in log space from -20 to 7
## could have a unique theta for each bin and/or experiment

n_c = 10 ## number of bins 
n_s = 20 ## number of samples to take from each bin
n_r = 1000 ## number of times to repeat the experiment

Size = 1:n_s
Results_rcz = array(NA, dim=c(n_r,n_c,2))

for( rI in 1:n_r ){
  prop_c = exp(rnorm(n_c)) ## 10 lognormal proportions
  prop_c = prop_c / sum(prop_c) ## rescale to get "true" population comp in bin (sums to 1)
  
  ## can't use Sapply because of the class returned by rdirichlet
  Y_sc = matrix(NA, nrow = n_s, ncol = n_c)
  for(i in Size){
    Y_sc[i,] <-  c(rdirichlet(1, alpha=prop_c*theta) ) 
  }
  
  # Convert to data frame
  Y_iz = data.frame(
    y = as.vector(Y_sc), ## number observed
    c = as.vector(col(Y_sc)), ## bin (e.g. age)
    s = as.vector(row(Y_sc)) ## sample id (akin to "haul", always 20)
  )
  
  # Fit Poisson-GLM 
  ## independent effects of haul and bin; based on your observations
  G = glm( y ~ 0 + factor(c) + factor(s), data=Y_iz, family=quasipoisson() )
  
  ## can hard calculate expected dispersion via the following, 
  ## or look within summary(G). dp's > 1 mean u are overdispersed
  dp = sum(residuals(G,type ="pearson")^2)/G$df.residual
  
  # Extract effects (with data, try to estimate true pop comp (prop_c))
  prophat_c = summary(G)$coef[paste0("factor(c)",1:n_c),'Estimate']
  prophat_c = exp(prophat_c) / sum(exp(prophat_c)) ## rescale, again
  
  # Record results
  Results_rcz[rI,,1] = prop_c ## true observed props
  Results_rcz[rI,,2] = prophat_c ## simulated/estimated props
}

## shows that we get even more error around observations
plot( x=Results_rcz[,,1], y=Results_rcz[,,2], log="",
      xlab = 'True',ylab = 'Simulated', main = 'Dirichlet with Theta = 0.05')
abline(a=0, b=1, col="blue")

## plot it in a clearer way
# Results_rcz has nrow = experiment, ncol = bins, arrays = true, obs
require(ggplot2); require(dplyr); require(reshape2)
plotme <-  data.frame(Results_rcz[,,1]) %>%
  mutate(haulid = 1:1000) %>%
  melt(id = 'haulid') %>%
  mutate(src = 'obs', 
         age = as.numeric(substr(variable,2,2))) %>%
  rbind( data.frame(Results_rcz[,,2]) %>%
           mutate(haulid = 1:1000) %>%
           melt(id = 'haulid') %>%
           mutate(src = 'pred', 
                  age = as.numeric(substr(variable,2,2)))) %>%
  select(-variable) %>%
  filter(haulid < 21) ## for plotting

dirich <- ggplot(NULL, aes(x = age, y = value, fill = src, color = src)) +
  geom_histogram(data = subset(plotme, src == 'obs', fill = 'obs'),stat = 'identity' ) +
  geom_line(data = subset(plotme, src == 'pred', color = 'pred'), lwd = 1.1) +
  scale_fill_manual(values = c('grey55',NA), labels = c('obs','simulated'))+
  scale_color_manual(values = c(NA, 'blue'), labels = c('obs','simulated'))+
  labs(x = 'age or length', y = 'freq', fill = '',
       color='', title = 'Dirichlet Simulator, 20 years of data')+
  theme_minimal()+
  facet_wrap(~haulid)

## comparison plots ----
par(mfrow = c(1,2))
plot( x=mn[,,1], y=mn[,,2], log="", cex.main = 0.9,
      xlab = 'True',ylab = 'Simulated', main = "Multinomial (Thorson's Code)")
abline(a=0, b=1, col="blue")
plot( x=Results_rcz[,,1], y=Results_rcz[,,2], log="", cex.main = 0.9,
      xlab = 'True',ylab = 'Simulated', main = 'Dirichlet with Theta (overdispersion) = 0.05')
abline(a=0, b=1, col="blue")


require(patchwork)
mn2  | dirich
