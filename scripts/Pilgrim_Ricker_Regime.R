# Pilgrim River sockeye salmon SR analysis

# Author: Adam Reimer
# Version: 2023-09-11

# Packages
packs <- c("tidyverse", "jagsUI", "mgcv")
lapply(packs, require, character.only = TRUE)
  
# Notes


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read data ---------------------------------------------------------------
# read in dat_chinBEG for profile comparisons
brood <- readxl::read_xlsx(".\\data\\brood_table.xlsx", sheet = 2)

#Jags data
data <- 
  list(
    S = brood$S,
    R = brood$R, 
    ar1 = 0,
    kf = 0,
    nyrs = length(brood$S),
    d = 4
)

# #   Simple Ricker Model ---------------------------------------------------------------
# cat("
#   model {
#     for (y in 1:nyrs) {
#         R[y] ~ dlnorm(mu[y], tau)
#         fit[y] <- log(S[y]) + lnalpha - (beta * S[y]/10^d)
#         e[y] <- log(R[y]) - fit[y]
#     }
#     mu[1] <- fit[1] + ar1 * phi * e0
#     cw[1] ~ dnorm(0, tau.w)
#     for (y in 2:nyrs) {
#         cw[y] ~ dnorm(cw[y - 1], tau.w)
#         mu[y] <- fit[y] + kf * cw[y] + ar1 * phi * e[y - 1]
#     }
#     lnalpha ~ dunif(0, 5)
#     beta ~ dunif(0, 5)
#     sigma ~ dunif(0, 2)
#     sigma.w ~ dunif(0, 2)
#     phi ~ dnorm(0, 4)T(-1, 1)
#     e0 ~ dnorm(0, 0.001)
#     tau <- 1/(sigma * sigma)
#     tau.w <- 1/(sigma.w * sigma.w)
#     for (y in 1:nyrs) {
#         lnalphai[y] <- lnalpha + cw[y]
#     }
#   }
#       ", file="Ricker_model.jag")
# 
# 
# #   * Chain params ------------------------------------------------------------
# #specify chains here because of the way we are doing ints. Better to write a function which automatically generates inits
# nchains <- 5
# nb <- 5000
# nt <- 5
# ns <- 10000
#  
# #   * Inits ---------------------------------------------------------------
# get_inits <- function(){
#   list(
#     beta = rlnorm(1, log(2e-5), 0.4),
#     lnalpha = rlnorm(1, log(1.6), 0.4),
#     e0 = rnorm(0, 0, 1),
#     phi = runif(1, 0.25, 0.75),
#     sigma = runif(1, 0, 1),
#     sigma.w = runif(1, 0, 1)
#   )
# }
# 
# 
# #   * Run jags ------------------------------------------------------------
# post <- jags(data = data,
#              parameters.to.save = c("lnalpha", "beta", "sigma"),
#                  inits = get_inits,
#                  model.file = "Ricker_model.jag",
#                  n.chains = nchains,
#                  n.iter = ns,
#                  n.burnin = nb,
#                  n.thin = nt,
#                  parallel = TRUE,
#                  store.data = TRUE
# )
# post


# Regime transition model -------------------------------------------------
cat("
  model {
    for (y in 1:nyrs) {
        R[y] ~ dlnorm(mu[y], tau[lambda[y]])
        fit[y] <- log(S[y]) + lnalpha[lambda[y]] - (beta[lambda[y]] * S[y]/10^d)
        e[y] <- log(R[y]) - fit[y]
        lambda[y] ~ dcat(gamma[y,1:2])
        lambda_1[y] <- lambda[y] - 1
    }
    mu[1] <- fit[1] + ar1 * phi * e0
    for (y in 2:nyrs) {
        mu[y] <- fit[y] + ar1 * phi * e[y - 1]
        for(r in 1:2){
          gamma[y, r] <- pi[lambda[y - 1], r]
        }
    }
    for (r in 1:2){
    #   lnalpha[r] ~ dunif(0, 5)
    #   beta[r] ~ dunif(0, 5)
    #  sigma[r] ~ dgamma(5, 15)
     tau[r] <- 1/(sigma[r] * sigma[r])
    }
    lnalpha[1] ~ dnorm(0, 1)T(lnalpha[2], )
    lnalpha[2] ~ dnorm(0, 1)T(0, )
    beta[1] ~ dnorm(0, 1)T(0, )
    beta[2] ~ dnorm(0, 1)T(beta[1], )
    sigma[1] ~ dgamma(5, 15)
    sigma[2] ~ dnorm(0, 1)T(sigma[1],)
    gamma[1, 1] <- 13/17
    gamma[1, 2] <- 4/17
    pi[1,1] <- p_1[1]
    pi[1,2] <- p_1[2]
    pi[2,1] <- p_2[1]
    pi[2,2] <- p_2[2]
    p_1[1:2] ~ ddirch(c(13, 2))
    p_2[1:2] ~ ddirch(c(2, 4))
    phi ~ dnorm(0, 4)T(-1, 1)
    e0 ~ dnorm(0, 0.001)
  }
      ", file="Ricker_hM_model.jag")


#   * Chain params ------------------------------------------------------------
#specify chains here because of the way ewe are doing intis. Better to write a function which automatically generates inits
nchains <- 5
nb <- 5000
nt <- 5
ns <- 20000


#   * Run jags ------------------------------------------------------------
post <- jags(data = data,
             parameters.to.save = c("lnalpha", "beta", "sigma", "pi", "lambda_1"),
             model.file = "Ricker_hM_model.jag",
             n.chains = nchains,
             n.iter = ns,
             n.burnin = nb,
             n.thin = nt,
             parallel = TRUE,
             store.data = TRUE
)


post
#saveRDS(post, ".\\post.rds")
#post <- readRDS(".\\post.rds")


#  * Ricker for each regime --------------------------------------------------
plot(1:90000, 
     1:90000*exp(post$q50$lnalpha[1])*exp(-post$q50$beta[1]/10^4*1:90000), 
     type = "l", 
     col = "red", 
     ylim = c(0, 7e4), 
     xlab = "S", 
     ylab = "R")
lines(1:90000, 1:90000*exp(post$q50$lnalpha[2])*exp(-post$q50$beta[2]/10^4*1:90000))
abline(0,1)
points(data$S, data$R)
text(data$S, data$R + 3000, round(post$mean$lambda_1, digits = 2))


#  * Sigma prior and posterior -----------------------------------------------
par(mfrow = c(3, 1))
hist(rgamma(1000, 5, 15), xlim = c(0, 2)) #prior
hist(post$sims.list$sigma[, 1], xlim = c(0, 2)) #post sigma 1
hist(post$sims.list$sigma[, 2], xlim = c(0, 2)) #post sigma 2
par(mfrow = c(1, 1))


#  * regime membership ----------------------------------------------------
post$sims.list$lambda_1
#distribution of broods in regime 1
membership <- 
  post$sims.list$lambda_1 %>% 
    rowSums() %>% 
    as.data.frame() %>% 
    setNames("R2") %>%
    mutate(R1 = 17 - R2) 
par(mfrow = c(2, 1))
hist(membership$R1, xlab = "# of broods from the high productivity regime", main = NULL)
hist(membership$R2, xlab = "# of broods from the low productivity regime", main = NULL)
par(mfrow = c(1, 1))
plot(post$sims.list$lnalpha[,1], post$sims.list$beta[,1], 
     col = ifelse(membership$R1 <= 1,'red','green'))
plot(post$sims.list$lnalpha[,2], post$sims.list$beta[,2], 
     col = ifelse(membership$R2 <= 1,'red','green'))

# Simulate yield across regimes ------------------------------------------------
#indices to use to draw 1000 samples from posterior
draw <- sample(1:dim(post$sims.list$lnalpha)[1],
               size = 1000,
               replace = FALSE)

#draw pi
sims_pi <- list(pi = post$sims.list$pi[draw, , ])
stationary_pi <- list()
stationary_pi$pi <- array(NA, dim(sims_pi$pi))
for(i in 1:dim(sims_pi$pi)[1]){
  stationary_pi$pi[i,,] <- 
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,] %*%
    sims_pi$pi[i,,]
}
head(stationary_pi$pi)
apply(stationary_pi$pi, MARGIN = c(2:3), median)
head(sims_pi$pi)

#function to simulate future regimes
get_regime <- function(sample, sims){
  regime <- 
    data.frame(
      year = 1:100,
      rand = runif(100, 0, 1),
      phi_11 = rep(diag(sims$pi[sample,,])[1], 100),
      phi_22 = rep(diag(sims$pi[sample,,])[2], 100),
      regime = c(sample(1:2, 1, prob = c(12/17, 5/17)), rep(NA, 99)))
  
  for(i in 2:100){
    if(regime$regime[i - 1] == 1){
      if(regime$rand[i] <= regime$phi_11[i]){regime$regime[i] = 1}
      else{regime$regime[i] = 2}
      }
    else{
      if(regime$rand[i] <= regime$phi_22[i]){regime$regime[i] = 2}
      else{regime$regime[i] = 1}
    }
  }
  as.integer(regime$regime)
}

#draw lnalpha, beta and sigma and format to join regime dataframe
sims_lnalpha <-
  data.frame(sample = 1:1000,
             lnalpha = post$sims.list$lnalpha[draw, ]) %>%
  pivot_longer(cols = starts_with("lnalpha"), 
               names_to = "regime", 
               names_prefix = "lnalpha.", 
               values_to = "lnalpha") %>%
  mutate(regime = as.integer(regime))

sims_beta <-
  data.frame(sample = 1:1000,
             beta = post$sims.list$beta[draw, ]) %>%
  pivot_longer(cols = starts_with("beta"), 
               names_to = "regime", 
               names_prefix = "beta.", 
               values_to = "beta") %>%
  mutate(regime = as.integer(regime))

sims_sigma <-
  data.frame(sample = 1:1000,
             sigma = post$sims.list$sigma[draw, ]) %>%
  pivot_longer(cols = starts_with("sigma"), 
               names_to = "regime", 
               names_prefix = "sigma.", 
               values_to = "sigma") %>%
  mutate(regime = as.integer(regime))

#create a dataframe with 100 years of future regimes for each 1000 samples
regimes <- 
  lapply(1:1000, function(x){
    data.frame(sample = x, year = 1:100, regime = get_regime(x, stationary_pi))
    }) %>%
  do.call(rbind, .) 

#join with appropriate SR parameter estimates to give SR parameter estimates 
# by sample and year after accounting for regime membership 
SRdata <- 
  regimes %>%
  left_join(sims_lnalpha, by = c("sample", "regime")) %>%
  left_join(sims_beta, by = c("sample", "regime")) %>%
  left_join(sims_sigma, by = c("sample", "regime")) %>%
  mutate(epsilon = exp(rnorm(n(), 0, sigma))) 

#Calculate yield and recruitment for a wide range of S across samples and years
#R_epsilon and Y_epsilon used to estimate EG range using simulation and GLMM fit.
#Rp, Yp, Smsy_p, MSY_p, OYP_90_p for mean based OYP
#R, Y, Smsy, MSY, OYP_90 for median based OYP
est_Y <-
  SRdata %>%
  slice(rep(1:n(), each = 100)) %>%
  mutate(lnalpha_p = lnalpha + sigma * sigma / 2,
         S = rep(seq(400, 40000, length.out = 100), times = 100 * 1000), #2),
         R = S * exp(lnalpha - beta / 10^data$d * S),
         Rp = S * exp(lnalpha_p - beta / 10^data$d * S),
         R_epsilon = R * epsilon,
         Y = R - S,
         Yp = Rp - S,
         Y_epsilon = R_epsilon - S,
         Smsy = lnalpha / (beta / 10^data$d) * (0.5 - 0.07 * lnalpha),
         Smsy_p = lnalpha_p / (beta / 10^data$d) * (0.5 - 0.07 * lnalpha_p),
         MSY = Smsy * exp(lnalpha - beta / 10^data$d * Smsy) - Smsy,
         MSY_p = Smsy_p * exp(lnalpha - beta / 10^data$d * Smsy_p) - Smsy_p,
         OYP_90 = ifelse(Y > 0.9 * MSY, 1, 0),
         OYP_90_p = ifelse(Yp > 0.9 * MSY_p, 1, 0),
         Smax = 1 / (beta / 10^data$d),
         MSR = Smax * exp(lnalpha - beta / 10^data$d * Smax),
         ORP_90 = ifelse(R > 0.9 * MSR, 1, 0))


# Simulation / GLMM eg ranges ---------------------------------------------
#build dataset
sim_Y <- 
  est_Y %>%
    group_by(S) %>%
    summarise(Y_mean = mean(Y_epsilon),
              Y_median = median(Y_epsilon),
              Y_10 = quantile(Y_epsilon, 0.1),
              Y_90 = quantile(Y_epsilon, 0.9),
              R_mean = mean(R_epsilon),
              R_median = median(R_epsilon))

#  * Smsy and eg range (mean) ------------------------------------------
mod_Ymean <- gam(Y_mean ~ s(S), data = sim_Y)
fit_Ymean <- predict(mod_Ymean, newdata = data.frame(S = 1:40000)) %>% round()
#Smsy
(index_Smsy <- fit_Ymean[which.max(fit_Ymean)])
#eg bounds based on 90% of MSY
(eg90 <- 
  fit_Ymean[fit_Ymean > round(max(fit_Ymean) * 0.9)] %>% 
  names() %>%
  as.numeric() %>%
  range())

#Sustained yield plot accounting for regime changes
plot(sim_Y$S, sim_Y$Y_mean, ylim = c(0, 100000))
lines(sim_Y$S, fit_Ymean[sim_Y$S]) #GLM fit
lines(sim_Y$S, sim_Y$Y_10, col = "red")
lines(sim_Y$S, sim_Y$Y_90, col = "red")
rect(eg90[1], 0, eg90[2], 1E5, #simulation/GLMM based goal (Y mean)
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))


#  * Smsy and eg range (median) ------------------------------------------
mod_Ymedian <- gam(Y_median ~ s(S), data = sim_Y)
fit_Ymedian <- predict(mod_Ymedian, newdata = data.frame(S = 1:40000)) %>% round()
#Smsy
(index_Smsy_median <- fit_Ymedian[which.max(fit_Ymedian)])
#eg bounds based on 90% of MSY
(eg90_median <-
   fit_Ymedian[fit_Ymedian > round(max(fit_Ymedian) * 0.9)] %>%
  names() %>%
  as.numeric() %>%
  range())

#Sustained yield plot accounting for regime changes
plot(sim_Y$S, sim_Y$Y_median, ylim = c(0, 100000))
lines(sim_Y$S, sim_Y$Y_median)
lines(sim_Y$S, sim_Y$Y_10, col = "red")
lines(sim_Y$S, sim_Y$Y_90, col = "red")
rect(eg90_median[1], 0, eg90_median[2], 1E5, #simulation/GLMM based goal (Y median)
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
rect(6400, 0, 11300, 1E5, #SR shiny app goal (Y median)
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))

#  * Smax and eg range (mean) ------------------------------------------
mod_Rmean <- gam(R_mean ~ s(S), data = sim_Y)
fit_Rmean <- predict(mod_Rmean, newdata = data.frame(S = 1:40000)) %>% round()
#Smax
fit_Rmean[which.max(fit_Rmean)]
#eg bounds based on 90% of Rmax
fit_Rmean[fit_Rmean > round(max(fit_Rmean) * 0.9)] %>% 
  names() %>%
  as.numeric() %>%
  range()
#Ricker accounting for regime changes
plot(sim_Y$S, sim_Y$R_mean)
lines(sim_Y$S, sim_Y$R_median)
abline(0,1)
#with yield based eg range
rect(eg90[1], 0, eg90[2], 1E5, #simulation/GLMM based goal (Y mean)
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
rect(6400, 0, 11300, 5E5, #SR shiny app goal (Y median)
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))

#  * Smax and eg range (median) ------------------------------------------
# mod_Rmedian <- gam(R_median ~ s(S), data = sim_Y)
# fit_Rmedian <- predict(mod_Rmedian, newdata = data.frame(S = 1:40000)) %>% round()
# fit_Rmedian[which.max(fit_Rmedian)]
# fit_Rmedian[fit_Rmedian > round(max(fit_Rmedian) * 0.9)] %>% 
#   names() %>%
#   as.numeric() %>%
#   range()



# OYP ---------------------------------------------------------------------
#Smsy quantiles
est_Y %>%
  group_by(regime) %>%
  summarize(min_Smsy = min(Smsy),
            q10_Smsy = quantile(Smsy, 0.05),
            q50_Smsy = quantile(Smsy, 0.5),
            q90_Smsy = quantile(Smsy, 0.95),
            max_Smsy = max(Smsy))

est_OYP_regime <- 
  est_Y %>%
  group_by(regime, S) %>%
  summarise(Y_mean = mean(Y),
            Yp_mean = mean(Yp),
            Y_mean_epsilon = mean(Y_epsilon),
            Y_median = median(Y),
            Y_median_epsilon = median(Y_epsilon),
            Y_10 = quantile(Y, 0.1),
            Y_90 = quantile(Y, 0.9),
            R_mean = mean(R),
            R_median = median(R),
            OYP_90 = mean(OYP_90),
            OYP_90_p = mean(OYP_90_p),
            ORP_90 = mean(ORP_90))

est_OYP <- 
  est_Y %>%
  group_by(S) %>%
  summarise(Y_mean = mean(Y),
            Yp_mean = mean(Yp),
            Y_mean_epsilon = mean(Y_epsilon),
            Y_median = median(Y),
            Y_median_epsilon = median(Y_epsilon),
            Y_10 = quantile(Y, 0.1),
            Y_90 = quantile(Y, 0.9),
            R_mean = mean(R),
            R_median = median(R),
            OYP_90 = mean(OYP_90),
            OYP_90_p = mean(OYP_90_p),
            ORP_90 = mean(ORP_90))


# OYP median --------------------------------------------------------------
plot(est_OYP$S, est_OYP$OYP_90, type = "l", ylim = c(0, 1))
abline(h = 0.9, col = "red")
#simulation/GLMM based goal (Y median)
rect(eg90_median[1], 0, eg90_median[2], 1,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
#SR shiny app goal (Y median)
rect(6400, 0, 11300, 1,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))
lines(est_OYP_regime$S[est_OYP_regime$regime == 1], 
      est_OYP_regime$OYP_90[est_OYP_regime$regime == 1], col = "green")
lines(est_OYP_regime$S[est_OYP_regime$regime == 2], 
      est_OYP_regime$OYP_90[est_OYP_regime$regime == 2], col = "orange")
#Max S from OYP for combined regimes
est_OYP[which.max(est_OYP$OYP_90), ]
#OYP prob at shiny app lb
est_OYP[est_OYP$S == 6400, ]
est_OYP_regime[est_OYP_regime$S == 6400, ]
#S where regime 1 OYP_90 == 90
est_OYP_regime[est_OYP_regime$OYP_90 > 0.85 & est_OYP_regime$regime == 1, ]
#regime 2 at same S
est_OYP_regime[est_OYP_regime$S %in% c(7200, 7600), ]


# OYP prime (mean) --------------------------------------------------------
plot(est_OYP$S, est_OYP$OYP_90_p, type = "l", ylim = c(0, 1))
abline(h = 0.9, col = "red")
rect(eg90[1], 0, eg90[2], 1,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
rect(6400, 0, 11300, 1,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))
lines(est_OYP_regime$S[est_OYP_regime$regime == 1], 
      est_OYP_regime$OYP_90_p[est_OYP_regime$regime == 1], col = "green")
lines(est_OYP_regime$S[est_OYP_regime$regime == 2], 
      est_OYP_regime$OYP_90_p[est_OYP_regime$regime == 2], col = "orange")

#Sustained yield plot each analysis type
#Note Y gives mean, Yp gives median as verified by simulation
#mean or median of Y negligible difference
plot(est_OYP$S, 
     est_OYP$Y_mean, 
     ylim = c(0, 60000),
     type = "l")
lines(est_OYP$S, #similar to mean
     est_OYP$Y_median, 
     ylim = c(0, 60000))
lines(est_OYP$S, 
      est_OYP$Yp_mean,  
      ylim = c(0, 60000),
      lt = 3)
lines(est_OYP$S, 
      est_OYP$Y_mean_epsilon,  
      ylim = c(0, 60000),
      lt = 2, 
      col = "purple")
lines(est_OYP$S, 
      est_OYP$Y_median_epsilon,  
      ylim = c(0, 60000),
      lt = 2, 
      col = "purple")

# Note mean median diff for regime 1 is small. 
# regime 2 mean is near regime 1
plot(est_OYP_regime$S[est_OYP_regime$regime == 1], 
     est_OYP_regime$Y_mean[est_OYP_regime$regime == 1], 
     ylim = c(0, 60000),
     type = "l",
     col = "green")
lines(est_OYP_regime$S[est_OYP_regime$regime == 1], 
      est_OYP_regime$Yp_mean[est_OYP_regime$regime == 1], 
      ylim = c(0, 60000),
      col = "green", 
      lt = 2)
lines(est_OYP_regime$S[est_OYP_regime$regime == 2], 
      est_OYP_regime$Y_mean[est_OYP_regime$regime == 2], 
      ylim = c(0, 60000),
      col = "orange")
lines(est_OYP_regime$S[est_OYP_regime$regime == 2], 
      est_OYP_regime$Yp_mean[est_OYP_regime$regime == 2], 
      ylim = c(0, 60000),
      col = "orange", 
      lt = 2)


#Yield Horsetail (median)
plot(NULL, xlim=c(0,50000), ylim=c(0,75000), ylab="R", xlab="S")
for(i in 1:30){
  lines(1:50000, 
        1:50000*exp(sims_lnalpha$lnalpha[sims_lnalpha$regime == 1 & sims_lnalpha$sample == i])*
          exp(-sims_beta$beta[sims_beta$regime == 1 & sims_beta$sample == i]/10^4*1:50000) - 1:50000,
        col = "green")
  lines(1:50000, 
        1:50000*exp(sims_lnalpha$lnalpha[sims_lnalpha$regime == 2 & sims_lnalpha$sample == i])*
          exp(-sims_beta$beta[sims_beta$regime == 2 & sims_beta$sample == i]/10^4*1:50000) - 1:50000,
        col = "orange")
}
lines(1:50000, 
      1:50000*exp(post$q50$lnalpha[1])*exp(-post$q50$beta[1]/10^4*1:50000) - 1:50000, 
      lwd = 4, 
      ylim = c(0, 7e4))
lines(1:50000, 
      1:50000*exp(post$q50$lnalpha[2])*exp(-post$q50$beta[2]/10^4*1:50000) - 1:50000, 
      lwd = 4)
rect(eg90[1], 0, eg90[2], 75000,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
rect(6400, 0, 11300, 7.5E4,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))

#Yield Horsetail (mean)
SRdata_p <- SRdata %>% mutate(lnalpha_p = lnalpha + sigma * sigma / 2)
plot(NULL, xlim=c(0,50000), ylim=c(0,75000), ylab="R", xlab="S")
for(i in 1:30){
  lines(1:50000, 
        1:50000*exp(SRdata_p$lnalpha_p[SRdata_p$regime == 1 & SRdata_p$sample == i])*
          exp(-SRdata_p$beta[SRdata_p$regime == 1 & SRdata_p$sample == i]/10^4*1:50000) - 1:50000,
        col = "green")
  lines(1:50000, 
        1:50000*exp(SRdata_p$lnalpha_p[SRdata_p$regime == 2 & SRdata_p$sample == i])*
          exp(-SRdata_p$beta[SRdata_p$regime == 2 & SRdata_p$sample == i]/10^4*1:50000) - 1:50000,
        col = "orange")
}
lnalpha_p <- post$q50$lnalpha + post$q50$sigma * post$q50$sigma / 2
lines(1:50000, 
      1:50000*exp(lnalpha_p[1])*exp(-post$q50$beta[1]/10^4*1:50000) - 1:50000, 
      lwd = 4, 
      ylim = c(0, 7e4))
lines(1:50000, 
      1:50000*exp(lnalpha_p[2])*exp(-post$q50$beta[2]/10^4*1:50000) - 1:50000, 
      lwd = 4)
rect(eg90[1], 0, eg90[2], 75000,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
rect(6400, 0, 11300, 7.5E4,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))

#ORP
plot(est_OYP$S, est_OYP$ORP_90, type = "l", ylim = c(0, 1))
abline(h = 0.9, col = "red")
#simulation/GLMM based goal (Y median)
rect(eg90_median[1], 0, eg90_median[2], 1,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
#SR shiny app goal (Y median)
rect(6400, 0, 11300, 1,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))
lines(est_OYP_regime$S[est_OYP_regime$regime == 1], 
      est_OYP_regime$ORP_90[est_OYP_regime$regime == 1], col = "green")
lines(est_OYP_regime$S[est_OYP_regime$regime == 2], 
      est_OYP_regime$ORP_90[est_OYP_regime$regime == 2], col = "orange")

#Horsetail
plot(NULL, xlim=c(0,90000), ylim=c(0,70000), ylab="R", xlab="S")
for(i in 1:30){
  lines(1:90000, 
        1:90000*exp(sims_lnalpha$lnalpha[sims_lnalpha$regime == 1 & sims_lnalpha$sample == i])*
          exp(-sims_beta$beta[sims_beta$regime == 1 & sims_beta$sample == i]/10^4*1:90000),
        col = "green")
  lines(1:90000, 
        1:90000*exp(sims_lnalpha$lnalpha[sims_lnalpha$regime == 2 & sims_lnalpha$sample == i])*
          exp(-sims_beta$beta[sims_beta$regime == 2 & sims_beta$sample == i]/10^4*1:90000),
        col = "orange")
}
lines(1:90000, 
      1:90000*exp(post$q50$lnalpha[1])*exp(-post$q50$beta[1]/10^4*1:90000), 
      lwd = 4, 
      ylim = c(0, 7e4))
lines(1:90000, 
      1:90000*exp(post$q50$lnalpha[2])*exp(-post$q50$beta[2]/10^4*1:90000), 
      lwd = 4)
abline(0, 1)
#simulation/GLMM based goal (Y median)
rect(eg90_median[1], 0, eg90_median[2], 75000,
     border = NA,
     col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.2))
#SR shiny app goal (Y median)
rect(6400, 0, 11300, 7.5E4,
     border = NA,
     col = rgb(red = 0, green = 0.5, blue = 0, alpha = 0.2))

