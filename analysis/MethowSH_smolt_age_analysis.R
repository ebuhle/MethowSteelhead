#===============================================================
# METHOW HATCHERY STEELHEAD ANALYSIS
# 
# To do:
# - RMarkdown vignette figure numbers & captions (bookdown?)
# - CJS PPD
# - Fix p_TDAWEA = 1 for *all* capture histories
# - (Put this on GitHub?)
#
# Issues:
# - Size data contain one tag (3D9.1C2CA3C9DD) that is not in
#   the capture history data. (Many tags in the CH data are
#   not in the size data; presumably these weren't measured?)
# - Should CJS allow time-varying covariate effects?
#   => No for now; see if size explains some interannual var.
# - Should we consider alternative groupings for the
#   time-varying effects in CJS? E.g., brood year, or upstream
#   survivals grouped by release year?
#   => No
# - How about effects of smolt age on upstream survival?
#   => No
# - Or effects of smolt age on downstream detection?
#   => No, S1 & S2 probably have similar detection probs
# - Why aren't we including smolt length, given its crucial
#   role in the previous analysis?
#   => We are now
# - SAR should depend on ocean age somehow, but oa is undefined
#   (not to mention unobserved) for non-survivors...? IPM!
#   => In full conditional posterior, factors for age and
#      survival are independent (?). "Dis-integrated model"?
#===============================================================

setwd(file.path("~","MethowSteelhead","analysis"))
options(device = windows)
library(lubridate)
library(gtools)
library(Hmisc)
library(vioplot)
source("vioplot2.R")
library(yarrr)
library(rstan)
library(rstanarm)
library(loo)
library(shinystan)
source("extract1.R")
source("stan_mean.R")
source("sim_phi_tot.R")
if(file.exists("MethowSH_stanfit.RData")) load("MethowSH_stanfit.RData")


#-----------------
# DATA
#-----------------

# Read in capture history data (observations are dates)
methowSH <- read.csv("methowSH.csv", header = T)
for(i in 8:29)
  methowSH[,i] <- as.Date(as.character(methowSH[,i]), "%m/%d/%Y")

# Read in size data, change unobserved values and "outliers" to NA,
# and merge into capture history data by tag ID
methowSHsize <- read.csv("methowSHsize.csv", header = T)
methowSHsize$length_rel[methowSHsize$error==1 | methowSHsize$type != 0] <- NA
methowSHsize[,c("length_rel","length_tag")] <- methowSHsize[,c("length_rel","length_tag")]/10 # convert to cm
methowSH <- data.frame(methowSH[,1:4], 
                       methowSHsize[match(methowSH$tag, methowSHsize$tag), c(2,6)],
                       methowSH[,-(1:4)])

# "Fix" return year and age at return

# (1) Some tags were detected in adult ladders the same year they were released
# (usually in April). Change those detections to NA.
rty_test <- sweep(apply(methowSH[,15:31], 2, year), 1, methowSH$release_year, "==") 
rty_indx <- which(rty_test & !is.na(rty_test), arr.ind = T)
rty_indx[,"col"] <- rty_indx[,"col"] + 14
methowSH[rty_indx] <- NA

# (2) return_year and adult_age are missing for some  fish that did, in fact, return as adults. 
# Fill in those missing values.
rty <- apply(methowSH[,15:31], 1, function(x) ifelse(all(is.na(x)), NA, min(year(x), na.rm = T)))
rty_NA <- !is.na(rty) & (is.na(methowSH$return_year) | is.na(methowSH$adult_age))
methowSH$return_year[rty_NA] <- rty[rty_NA]
methowSH$adult_age[rty_NA] <- methowSH$return_year[rty_NA] - methowSH$brood_year[rty_NA]

# Convert dates to 0/1
methowSH[,10:31] <- as.numeric(!is.na(methowSH[,10:31]))

# Pool detections from MCJ:TWX (juvenile), TDA:WEA (adult), and MRC:BROOD (adult)
methowSH <- cbind(methowSH[,1:8],
                  ocean_age = factor(methowSH$return_year - methowSH$release_year),
                  methowSH[,9:14],
                  MCJTWX = as.numeric(apply(methowSH[,c("MCJ","JDJ","BON","TWX")] > 0, 1, any)),
                  methowSH[,15:21],
                  TDAWEA = as.numeric(apply(methowSH[,c("TDA","MCN","PRA","RIA","RRF","WEA")] > 0, 1, any)),
                  methowSH[,22:31],
                  MRCBRD = as.numeric(apply(methowSH[,c("MRC","MRT","MRW","SPRING","WFC","LOR","EWC","CRW","BROOD")] > 0, 1, any)))
levels(methowSH$ocean_age) <- c("1","2+","2+")  # very few 3-ocean; group with 2-ocean

# Convert crosstabs of ocean age by release year and smolt age to data frame
# for use in posterior predictive checking of binomial GLMMs
methowSHoa <- ftable(release_year = methowSH$release_year, 
                     smolt_age = methowSH$smolt_age, 
                     ocean_age = methowSH$ocean_age, 
                     row.vars = 1:2, col.vars = 3, exclude = NULL)
methowSHoa <- as.data.frame(as.matrix(methowSHoa))
methowSHoa <- cbind(t(sapply(strsplit(row.names(methowSHoa), "_"), identity)),
                         methowSHoa)
names(methowSHoa) <- c("release_year","smolt_age","oa1","oa2plus","mort")
methowSHoa$release_year <- as.numeric(as.character(methowSHoa$release_year))
methowSHoa <- data.frame(methowSHoa[,1:2], smolt_age_num = as.numeric(methowSHoa$smolt_age),
                         methowSHoa[,3:5])
row.names(methowSHoa) <- NULL

# Convert to m-array format, including the covariates (some redundant):
# smolt_age, brood_year, release_year, return_year, adult_age, adult_age_factor
# (which combines 4 and 5 and arbitrarily fills NAs)
methowSHm <- cbind(methowSH[,1:4], 
                   return_year = factor(methowSH$return_year, exclude = NULL),
                   adult_age = factor(methowSH$adult_age, exclude = NULL),
                   adult_age_factor = factor(methowSH$adult_age, exclude = NULL),
                   methowSH[,c("WNFH","RRJ","MCJTWX","BOA", "TDAWEA","LMR","MRCBRD")])
levels(methowSHm$adult_age_factor) <- c("2","3","4+","4+","2")
methowSHm <- aggregate(tag ~ smolt_age + brood_year + release_year + return_year + adult_age + adult_age_factor +
                         WNFH + RRJ + MCJTWX + BOA + TDAWEA + LMR + MRCBRD,
                       data = methowSHm, length)
names(methowSHm)[names(methowSHm)=="tag"] <- "n"
methowSHm$return_year <- as.numeric(as.character(methowSHm$return_year))
methowSHm$adult_age <- as.numeric(as.character(methowSHm$adult_age))

# Year indices as grouping variables for phi and p random effects
# Replace NA return years with arbitrary index that will not enter into the likelihood
# (these capture histories contain no information on upstream survival and detection probs)
release_year <- as.numeric(factor(methowSHm$release_year))
return_year <- ifelse(!is.na(methowSHm$return_year),
                      as.numeric(factor(methowSHm$return_year)),
                      0) # code NA as zero

# phi: downstream, SAR grouped by release year; upstream grouped by return year
group_phi <- cbind(matrix(rep(release_year, 3), ncol = 3), matrix(rep(return_year, 3), ncol = 3))

# p: downstream grouped by release year; upstream grouped by return year
group_p <- cbind(matrix(rep(release_year, 3), ncol = 3), matrix(rep(return_year, 4), ncol = 4))

# Centered predictors for CJS models
smolt_age <- scale(as.numeric(methowSHm$smolt_age), scale=F)
adult_age <- na.replace(scale(methowSHm$adult_age, scale = F), 0)  # arbitrary NA value


#--------------------------------------------------------------
# BINOMIAL (BERNOULLI) GLMMs FOR OCEAN AGE
# Condition on survival to adulthood (and detection)
# Predict ocean age given smolt age, grouped by release year
#--------------------------------------------------------------

# Null model
# Intercept varies by release year
fit_oa0 <- stan_glmer(ocean_age ~ 1 | release_year, 
                      data = methowSH,
                      family = binomial("logit"),
                      prior_intercept = normal(0,1.5), 
                      prior_covariance = decov(),
                      chains = 4, iter = 2000, warmup = 1000, cores = 4)

summary(fit_oa0)


# ocean_age ~ smolt_age
# Intercept varies by release year
fit_oa1a <- stan_glmer(ocean_age ~ smolt_age + (1 | release_year), 
                      data = methowSH,
                      family = binomial("logit"),
                      prior = normal(0,3),
                      prior_intercept = normal(0,1.5), 
                      prior_covariance = decov(),
                      chains = 4, iter = 2000, warmup = 1000, cores = 4)

summary(fit_oa1a)


# ocean_age ~ smolt_age
# Intercept and smolt_age effect both vary by release year, but are uncorrelated
# (Note that the uncorrelated specification entails tricking the formula parser
# into treating smolt_age as a binary 0/1 numeric variable)
fit_oa1b <- stan_glmer(ocean_age ~ smolt_age + (smolt_age_num || release_year), 
                       data = cbind(methowSH, smolt_age_num = as.numeric(methowSH$smolt_age)),
                       family = binomial("logit"),
                       prior = normal(0,3),
                       prior_intercept = normal(0,1.5), 
                       prior_covariance = decov(),
                       chains = 4, iter = 2000, warmup = 1000, cores = 4)

summary(fit_oa1b)


# Model selection
loo_oa <- list(fit_oa0 = loo(fit_oa0), fit_oa1a = loo(fit_oa1a), fit_oa1b = loo(fit_oa1b))
compare_oa <- compare_models(loos = loo_oa)
compare_oa[,"looic"] <- compare_oa[,"looic"] - min(compare_oa[,"looic"])
compare_oa[,c("looic","se_looic","p_loo","se_p_loo")]
# pairwise is often more discriminating
compair_oa <- vector("list", length(loo_oa) -1 ) 
for(i in 1:length(compair_oa))
{
  names(compair_oa)[i] <- paste(row.names(compare_oa)[c(i+1,i)], collapse = " vs. ")
  compair_oa[[i]] <- 2*compare_models(loos = loo_oa[row.names(compare_oa)[c(i+1,i)]])
  names(compair_oa[[i]])[1] <- "looic_diff"
}
compair_oa


#-----------------------------------------------------------------------------------
# CJS MODELS
# (WNFH) - RRJ - MCJTWX - BOA - TDAWEA - LMR - MRCBRD
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# All years pooled, no random effects or covariates
#-----------------------------------------------------------------------------------

cjs_all <- stan(file = "CJS-marray.stan",
                data = list(T = 7, M = nrow(methowSHm), 
                            y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                            n = methowSHm$n),
                pars = c("phi","p","lambda","LL"),
                chains = 3, cores = 3, iter = 2000, warmup = 1000,
                control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_all, pars = "LL", include = FALSE)

#-----------------------------------------------------------------------------------
# Random effects on capture and detection probs
# phi: downstream, SAR grouped by release year; upstream grouped by return year
# p: downstream grouped by release year; upstream grouped by return year
# If not detected as adult (group code 0) fix phi[,3:6] = 0 and p[,4:7] = 1
#-----------------------------------------------------------------------------------

### No covariates
cjs_fixNA0 <- stan(file = "CJS-marray-phiXRE-pXRE-fixNA.stan",
                   data = list(T = 7, M = nrow(methowSHm), K = 1,
                               X = matrix(1, nrow(methowSHm), 1),
                               indX_phi = matrix(1,1,6),
                               group_phi = group_phi, 
                               indX_p = matrix(1,1,7),
                               group_p = group_p,
                               y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                               n = methowSHm$n),
                   pars = c("beta","sigma","b","s","LL"),
                   chains = 3, cores = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_fixNA0, pars = "LL", include = FALSE)                             


### phi[1:3] ~ smolt_age
cjs_fixNA1 <- stan(file = "CJS-marray-phiXRE-pXRE-fixNA.stan",
                   data = list(T = 7, M = nrow(methowSHm), K = 2,
                               X = model.matrix(~ smolt_age),
                               indX_phi = rbind(1, c(1,1,1,0,0,0)),
                               group_phi = group_phi, 
                               indX_p = rbind(rep(1,7), 0),
                               group_p = group_p,
                               y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                               n = methowSHm$n),
                   pars = c("beta","sigma","b","s","LL"),
                   chains = 3, cores = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_fixNA1, pars = "LL", include = FALSE)                             


### phi[4:5] ~ adult_age (categorical: 2, 3, 4+; NA = 2)
cjs_fixNA2a <- stan(file = "CJS-marray-phiXRE-pXRE-fixNA.stan",
                   data = list(T = 7, M = nrow(methowSHm), K = 3,
                               X = model.matrix(~ adult_age_factor, data = methowSHm),
                               indX_phi = rbind(1, c(0,0,0,1,1,0), c(0,0,0,1,1,0)),
                               group_phi = group_phi, 
                               indX_p = rbind(rep(1,7), 0, 0),
                               group_p = group_p,
                               y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                               n = methowSHm$n),
                   pars = c("beta","sigma","b","s","LL"),
                   chains = 3, cores = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_fixNA2a, pars = "LL", include = FALSE)                             


### phi[4:5] ~ adult_age (continuous; NA = 2)
cjs_fixNA2b <- stan(file = "CJS-marray-phiXRE-pXRE-fixNA.stan",
                   data = list(T = 7, M = nrow(methowSHm), K = 2,
                               X = model.matrix(~ adult_age),
                               indX_phi = rbind(1, c(0,0,0,1,1,0)),
                               group_phi = group_phi, 
                               indX_p = rbind(rep(1,7), 0),
                               group_p = group_p,
                               y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                               n = methowSHm$n),
                   pars = c("beta","sigma","b","s","LL"),
                   chains = 3, cores = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_fixNA2b, pars = "LL", include = FALSE)                             


### phi[1:3] ~ smolt_age, phi[4:5] ~ adult_age (continuous; NA = 2)
cjs_fixNA3 <- stan(file = "CJS-marray-phiXRE-pXRE-fixNA.stan",
                    data = list(T = 7, M = nrow(methowSHm), K = 3,
                                X = model.matrix(~ smolt_age + adult_age),
                                indX_phi = rbind(1, c(1,1,1,0,0,0), c(0,0,0,1,1,0)),
                                group_phi = group_phi, 
                                indX_p = rbind(rep(1,7), 0, 0),
                                group_p = group_p,
                                y = methowSHm[,c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")],
                                n = methowSHm$n),
                    pars = c("beta","sigma","epsilon_z","b","s","e_z","LL"),
                    chains = 3, cores = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.95, max_treedepth = 12))

print(cjs_fixNA3, pars = c("epsilon_z","e_z","LL"), include = FALSE)                             


### Model selection

# Compare LOO for individual capture histories
# Because all histories within a cell of the m-array have the same likelihood,
# we only need to call loo() on the posterior samples of the unique log-likelihoods.
# This is waaaaay faster than calling it on the full individual-level matrix. 
# We can then edit the loo object so it corresponds to the individual history data.
loo_cjs <- lapply(list(cjs_all = cjs_all,
                       cjs_fixNA0 = cjs_fixNA0, 
                       cjs_fixNA1 = cjs_fixNA1, 
                       cjs_fixNA2a = cjs_fixNA2a, 
                       cjs_fixNA2b = cjs_fixNA2b, 
                       cjs_fixNA3 = cjs_fixNA3),
                  function(obj) {
                    LL <- sweep(as.matrix(obj,"LL"), 2, methowSHm$n, "/") # LL of unique CH
                    LOO <- loo(LL)
                    LOO$pointwise <- LOO$pointwise[rep(1:ncol(LL), times = methowSHm$n),]
                    LOO$pareto_k <- LOO$pareto_k[rep(1:ncol(LL), times = methowSHm$n)]
                    LOO[c("elpd_loo","p_loo","looic")] <- colSums(LOO$pointwise)
                    LOO[c("se_elpd_loo","se_p_loo","se_looic")] <- 
                      apply(LOO$pointwise,2,sd)*sqrt(nrow(LOO$pointwise)) # SE of *sum*
                    return(LOO)
                  })

compare_cjs <- compare(x = loo_cjs)
# deltas
compare_cjs[,"looic"] <- compare_cjs[,"looic"] - min(compare_cjs[,"looic"])
compare_cjs[,c("looic","se_looic","p_loo","se_p_loo")]
# pairwise comparisons are often more discriminating
compair_cjs <- vector("list", length(loo_cjs) - 1) 
for(i in 1:length(compair_cjs))
{
  names(compair_cjs)[i] <- paste(row.names(compare_cjs)[c(i+1,i)], collapse = " vs. ")
  compair_cjs[[i]] <- 2*compare(x = loo_cjs[row.names(compare_cjs)[c(i+1,i)]])
  names(compair_cjs[[i]])[1] <- "looic_diff"
}
compair_cjs


#---------------------------------------------
# Save workspace (only the stanfit objects)
#---------------------------------------------

save(list = ls()[substring(ls(),1,4) %in% c("cjs_","fit_")], file = "MethowSH_stanfit.RData")


#-----------------
# TABLES
#-----------------

# Crosstabs of smolt, adult and ocean age by release or return year
table(release_year = methowSH$release_year, smolt_age = methowSH$smolt_age)
addmargins(table(release_year = methowSH$release_year, adult_age = methowSH$adult_age), 2)
addmargins(table(return_year = methowSH$return_year, adult_age = methowSH$adult_age), 2)
table(release_year = methowSH$release_year, ocean_age = methowSH$return_year - methowSH$release_year)

# Crosstabs of smolt and adult age by reach and release or return year respectively
tab1 <- aggregate(cbind(WNFH, RRJ, MCJTWX, BOA, TDAWEA, LMR) ~ smolt_age + release_year, data = methowSH, sum)
tab1[,-(1:2)] <- round(sweep(tab1[,-(1:2)], 1, tab1$WNFH, "/"), 3)
print(tab1, digits = 2)

tab2 <- aggregate(cbind(BOA, TDAWEA, LMR) ~ adult_age + return_year, data = methowSH, sum)
tab2[,-(1:2)] <- round(sweep(tab2[,-(1:2)], 1, tab2$BOA, "/"), 3)
print(tab2, digits = 2)


#-------------------------------------------------
# FIGURES
#-------------------------------------------------

#------------------------
# GLMMs for ocean age
#------------------------

# Posterior predictive distribution and data (proportion OA 2+) by year

p_hat <- posterior_linpred(fit_oa1b, newdata = methowSHoa, transform = TRUE)
p_hyper <- posterior_linpred(fit_oa1b, transform = TRUE, re.form = ~ 0,
                             newdata = data.frame(release_year = rep(0,2),
                                                  smolt_age = c("S1","S2"),
                                                  smolt_age_num = c(1,2)))
n <- rowSums(methowSHoa[,c("oa1","oa2plus")])
p_samp <- rbinom(length(p_hat), prob = p_hat, size = rep(n, each = nrow(p_hat)))
p_samp <- sweep(matrix(p_samp, nrow(p_hat), ncol(p_hat)), 2, n, "/")
ry <- methowSHoa$release_year + ifelse(methowSHoa$smolt_age=="S1", -0.15, 0.15) 

dev.new()
# png(filename="pp_ocean_age.png", width=7, height=7, units="in", res=300, type="cairo-png")
plot(ry, methowSHoa$oa2plus/n, pch = "", ylim = c(0,1), yaxs = "i",
     xlab = "Release year", ylab = "Proportion returning at ocean age 2+",
     main = "Data and posterior predictive distribution",
     cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, las = 1)
vioplot2(p_samp, at = ry, quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, 
         add = TRUE, border = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"),
         lwd.quantile = c(1,2,1))
points(ry, methowSHoa$oa2plus/n, pch = 16, cex = 1.2, 
       col = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"))
text(unique(round(ry)), 0.01, adj = c(0.5,0), 
     labels = rowSums(matrix(n, ncol = 2, byrow = TRUE)))
legend("topleft", c("1","2"), title = "Smolt age", pch = 16, cex = 1.2, lty = 1,
       col = c("darkgray","black"))
# dev.off()


# Fitted values and sample estimates (proportion OA 2+) by year

dev.new()
# png(filename="fitted_ocean_age.png", width=7, height=7, units="in", res=300, type="cairo-png")
plot(ry, methowSHoa$oa2plus/n, pch = "", xlim = c(min(ry), max(ry) + 1), ylim = c(0,1), 
     xaxt = "n", yaxs = "i", xlab = "Release year", ylab = "Proportion returning at ocean age 2+",
     main = "Sample proportions and fitted values", las = 1,
     cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5)
axis(1, at = c(unique(round(ry))), labels = unique(round(ry)))
axis(1, at = max(round(ry)) + 1, tick = FALSE, labels = "hyper-\nmean")
vioplot2(cbind(p_hat, p_hyper), at = c(ry, rep(max(ry),2) + c(0.7,1)), 
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = c(ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"), "darkgray", "black"),
         lwd.quantile = c(1,2,1))
abline(v = max(ry) + 0.35, lty = 3)
points(ry, methowSHoa$oa2plus/n, pch = 16, cex = 1.2, 
       col = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"))
segments(ry, y0 = binconf(methowSHoa$oa2plus, n, alpha = 0.1)[,"Lower"], 
         y1 = binconf(methowSHoa$oa2plus, n, alpha = 0.1)[,"Upper"],
         col = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"))
text(unique(round(ry)), 0.01, adj = c(0.5,0), 
     labels = rowSums(matrix(n, ncol = 2, byrow = TRUE)))
legend("topleft", c("1","2"), title = "Smolt age", pch = 16, cex = 1.2, lty = 1,
       col = c("darkgray","black"))
# dev.off()

# (Marginal) adult age distribution of S1 vs. S2: observed and predicted 

dev.new()
# png(filename="adult_age.png", width=7, height=7, units="in", res=300, type="cairo-png")
par(mar = c(5,5,4,1))
histS1S2 <- table(smolt_age = methowSH$smolt_age, adult_age = methowSH$adult_age)
histS1S2 <- sweep(histS1S2, 1,rowSums(histS1S2), "/")
bp <- barplot(histS1S2, beside = TRUE, space = c(0.05,1), names = dimnames(histS1S2)$adult_age,
              ylim = c(0,1), cex.names = 1.4, cex.axis = 1.4, cex.lab = 1.8,
              border = c("darkgray","black"), col = c("darkgray","black"), 
              las = 1, axis.lty = 1, xpd = NA, density = 30, 
              xlab = "Adult age", ylab = "Proportion")
box()
mtext("Sample proportions and fitted values", side = 3, line = 1, 
      font = 2, cex = par("cex")*1.5)
vioplot2(cbind(1 - p_hyper, p_hyper), at = c(bp[1,1], bp[2,2], bp[1,2], bp[2,3]),
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.5, drawRect = FALSE, add = TRUE,
         border = rep(c("darkgray", "black"), 2), lwd.quantile = c(1,2,1))
legend("topright", c("1","2"), title = "Smolt age", pch = "", lty = 1, cex = 1.4,
       col = c("darkgray","black"))
# dev.off()


#------------------------
# CJS MODELS
#------------------------

# Posterior predictive checking for CJS models based on m-array cell proportions

cjs_ppd_iter <- 100
system.time(
  cjs_ppd <- sim_phi_tot(cjs = cjs_fixNA3, fit_oa = fit_oa1b, iter = cjs_ppd_iter, 
                         adult_age = adult_age,release_years = unique(methowSHm$release_year),
                         return_years = unique(na.omit(methowSHm$return_year)),
                         CH = methowSH[,c("tag","release_year","smolt_age","return_year","adult_age",
                                          "WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR","MRCBRD")])
)

nm <- merge(cjs_ppd$CHm_sim, methowSHm, all.x = TRUE)
nm$n[is.na(nm$n)] <- 0
mean(nm$n > apply(nm[,as.character(1:cjs_ppd_iter)], 1, quantile, 0.05) &
       nm$n < apply(nm[,as.character(1:cjs_ppd_iter)], 1, quantile, 0.95))

dev.new(width = 7, height = 7)
# png(filename="pp_m-array.png", width=7, height=7, units="in", res=300, type="cairo-png")
plot(apply(nm[,as.character(1:cjs_ppd_iter)], 1, median)/sum(nm$n), nm$n/sum(nm$n), 
     pch = "", cex = 1.2, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, las = 1, 
     xlim = range(apply(nm[,as.character(1:cjs_ppd_iter)], 1, median), nm$n)/sum(nm$n),
     ylim = range(apply(nm[,as.character(1:cjs_ppd_iter)], 1, median), nm$n)/sum(nm$n), 
     xlab = "Predicted m-array cell proportion", ylab = "Observed m-array cell proportion",
     main = "Data and posterior predictive distribution")
abline(0,1)
points(apply(nm[,as.character(1:cjs_ppd_iter)], 1, median)/sum(nm$n), nm$n/sum(nm$n), 
       pch = 16, cex = 1.5, col = transparent("darkgray", 0.5))
segments(x0 = apply(nm[,as.character(1:cjs_ppd_iter)], 1, quantile, 0.05)/sum(nm$n),
         x1 = apply(nm[,as.character(1:cjs_ppd_iter)], 1, quantile, 0.95)/sum(nm$n),
         y0 = nm$n/sum(nm$n), col = transparent("darkgray", 0.5))
# dev.off()


# Marginal posteriors (and priors) of survival parameters in "full" model

beta <- as.matrix(cjs_fixNA3,"beta")
beta0 <- beta[, paste0("beta[1,", 1:6, "]")]
beta_smolt_age <- beta[, paste0("beta[2,", 1:6, "]")]
beta_adult_age <- beta[, paste0("beta[3,", 1:6, "]")]
sigma <- as.matrix(cjs_fixNA3,"sigma")
epsilon_z <- matrix(get_posterior_mean(cjs_fixNA3,"epsilon_z")[,"mean-all chains"],
                      ncol = 6, byrow = T)
occ_names <- c("WNFH","RRJ","MCJTWX","BOA","TDAWEA","LMR")

dev.new(width = 4, height = 8)
# png(filename="posteriors_phi.png", width=4, height=8, units="in", res=300, type="cairo-png")
par(mfrow = c(4,1), mar = c(3,5,1,1), oma = c(2,0,0,0))
# phi intercepts
plot(1:6, rep(0.5,6), pch = "", xaxt = "n", ylim = c(0,1), yaxs = "i", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xlab = "", 
     ylab = "Baseline survival")
axis(1, at = 1:6, labels = occ_names, cex.axis = 1.2)
vioplot2(plogis(beta0[,1:5]), at = 1:5 + 0.5, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.5, drawRect = FALSE, add = TRUE)
for(j in 1:5)
  for(i in unique(group_phi[,j]))
    points(j + 0.5, median(plogis(beta0[,j] + sigma[,j] * epsilon_z[i,j])),
           pch = 16, cex = 1.2, col = transparent("darkgray", 0.5))
# phi interannual SDs
plot(1:6, rep(0,6), ylim = range(0, sigma[,1:5]),  pch = "", xaxt = "n", yaxs = "i",
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xlab = "", 
     ylab = bquote("Survival SD (" ~ sigma * ")"))
axis(1, at = 1:6, labels = occ_names, cex.axis = 1.2)
vioplot2(matrix(qnorm(seq(0.001, 0.999, length = 1000),0,3), 1000, 5), at = 1:5 + 0.5, 
         quantiles = c(0.75,0.975), lwd.quantile = c(1,2,1), col = NA, 
         border = transparent("darkgray", 0.3), wex = 0.5, drawRect = FALSE, add = TRUE)
vioplot2(sigma[,1:5], at = 1:5 + 0.5, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.5, drawRect = FALSE, add = TRUE)
# phi regression coefs for smolt age
plot(1:6, rep(0,6), ylim = range(beta_smolt_age[,1:5]),  
     pch = "", xaxt = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     xlab = "", ylab = bquote("Smolt age effect (" ~ beta[2] * ")"))
axis(1, at = 1:6, labels = occ_names, cex.axis = 1.2)
abline(h = 0, lty = 2)
vioplot2(matrix(qnorm(seq(0.001, 0.999, length = 1000),0,3), 1000, 3), at = 1:3 + 0.5, 
         quantiles = c(0.75,0.975), lwd.quantile = c(1,2,1), col = NA, 
         border = transparent("darkgray", 0.3), wex = 0.5, drawRect = FALSE, add = TRUE)
vioplot2(beta_smolt_age[,1:3], at = 1:3 + 0.5, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.5, drawRect = FALSE, add = TRUE)
# phi regression coefs for adult age
plot(1:6, rep(0,6), ylim = range(beta_adult_age[,1:5]),  
     pch = "", xaxt = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5, xpd = NA,
     xlab = "Location", ylab = bquote("Adult age effect (" ~ beta[3] * ")"))
axis(1, at = 1:6, labels = occ_names, cex.axis = 1.2)
abline(h = 0, lty = 2)
vioplot2(matrix(qnorm(seq(0.001, 0.999, length = 1000),0,3), 1000, 2), at = 4:5 + 0.5, 
         quantiles = c(0.75,0.975), lwd.quantile = c(1,2,1), col = NA, 
         border = transparent("darkgray", 0.3), wex = 0.5, drawRect = FALSE, add = TRUE)
vioplot2(beta_adult_age[,4:5], at = 4:5 + 0.5, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.5, drawRect = FALSE, add = TRUE)
# dev.off()


# Marginal posteriors (and priors) of detection parameters in "full" model

b <- as.matrix(cjs_fixNA3,"b")
b0 <- b[, paste0("b[1,", 1:7, "]")]
s <- as.matrix(cjs_fixNA3,"s")
e_z <- matrix(get_posterior_mean(cjs_fixNA3,"e_z")[,"mean-all chains"],
                    ncol = 7, byrow = T)

dev.new(width = 7, height = 7)
# png(filename="posteriors_p.png", width=7, height=7, units="in", res=300, type="cairo-png")
par(mfrow = c(2,1), mar = c(3,5,1,1), oma = c(2,0,0,0))
# p intercepts
plot(2:6, rep(0.5,5), pch = "", xaxt = "n", xlim = c(1.5,6.5), ylim = c(0,1), 
     xaxs = "i", yaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     xlab = "", ylab = "P(detection)")
axis(1, at = 2:6, labels = occ_names[2:6], cex.axis = 1.2)
vioplot2(plogis(b0[,2:6]), at = 2:6, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE)
for(j in 2:6)
  for(i in unique(group_p[,j]))
    points(j, median(plogis(b0[,j] + s[,j] * e_z[i,j])),
           pch = 16, cex = 1.2, col = transparent("darkgray", 0.5))
# p interannual SDs
plot(2:6, rep(0,5), xlim = c(1.5,6.5), ylim = range(0, s[,2:6]),  pch = "", 
     xaxt = "n", xaxs = "i", yaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.5, xpd = NA, 
     xlab = "Location", ylab = bquote("Detection SD (" * italic(s) * ")"))
axis(1, at = 2:6, labels = occ_names[2:6], cex.axis = 1.2)
vioplot2(matrix(qnorm(seq(0.001, 0.999, length = 1000),0,3), 1000, 5), at = 2:6, 
         quantiles = c(0.75,0.975), lwd.quantile = c(1,2,1), col = NA, 
         border = transparent("darkgray", 0.3), wex = 0.3, drawRect = FALSE, add = TRUE)
vioplot2(s[,2:6], at = 2:6, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE)
# dev.off()


# Expected survival by reach for S1 vs. S2 smolts in an average year

dev.new(width = 7, height = 7)
# png(filename="reach_phi_S1vS2.png", width=7, height=7, units="in", res=300, type="cairo-png")
plot(1:4, rep(0.5,4), pch = "", xaxt = "n", ylim = c(0,1), yaxs = "i", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, xlab = "Location", ylab = "Survival")
axis(1, at = 1:4, labels = occ_names[1:4], cex.axis = 1.2)
# S1
vioplot2(plogis(beta0[,1:3]), at = 1:3 + 0.33, quantiles = c(0.05,0.5,0.95), 
         lwd.quantile = c(1,2,1), col = NA, border = "darkgray", 
         wex = 0.2, drawRect = FALSE, add = TRUE)
# S2
vioplot2(plogis(beta0[,1:3] + beta_smolt_age[,1:3]), at = 1:3 + 0.66, 
         quantiles = c(0.05,0.5,0.95), lwd.quantile = c(1,2,1), 
         col = NA, wex = 0.2, drawRect = FALSE, add = TRUE)
legend("topright", c("1","2"), title = "Smolt age", pch = "", cex = 1.2, lty = 1,
       col = c("darkgray","black"))
# dev.off()


# Expected reach-specific survival vs. adult age in an average year

dev.new(width = 7, height = 7)
# png(filename="reach_phi_adult_age.png", width=7, height=7, units="in", res=300, type="cairo-png")
par(mfrow = c(2,1), mar = c(2,5,2,1), oma = c(2.1,0,0,0))
# BOA to TDAWEA
plot(2:5, rep(0.5,4), pch = "", xaxt = "n", xlim = c(1.8,5.2), ylim = c(0,1), yaxs = "i", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, 
     xlab = "", ylab = "Survival", main = bquote("BOA" %->% "  TDAWEA"))
axis(1, at = 2:5, cex.axis = 1.2)
vioplot2(plogis(beta0[,4] + beta_adult_age[,4] %*% t(2:5)), 
         at = 2:5, quantiles = c(0.05,0.5,0.95),  lwd.quantile = c(1,2,1), 
         col = NA, wex = 0.3, drawRect = FALSE, add = TRUE)
# TDAWEA to LMR
plot(2:5, rep(0.5,4), pch = "", xaxt = "n", xlim = c(1.8,5.2), ylim = c(0,1), yaxs = "i", 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, xpd = NA,
     xlab = "Adult age", ylab = "Survival", main = bquote("TDAWEA" %->% "  LMR"))
axis(1, at = 2:5, cex.axis = 1.2)
vioplot2(plogis(beta0[,5] + beta_adult_age[,5] %*% t(2:5)), 
         at = 2:5, quantiles = c(0.05,0.5,0.95),  lwd.quantile = c(1,2,1), 
         col = NA, wex = 0.3, drawRect = FALSE, add = TRUE)
# dev.off()



# Time series of life-cycle survival probabilities for S1 and S2 smolts

phi_r2a <- sim_phi_tot(cjs = cjs_fixNA3, fit_oa = fit_oa1b, adult_age = adult_age,
                        release_years = unique(methowSHm$release_year),
                        return_years = unique(na.omit(methowSHm$return_year)), Tmax = 3)
phi_tot <- sim_phi_tot(cjs = cjs_fixNA3, fit_oa = fit_oa1b, adult_age = adult_age,
                        release_years = unique(methowSHm$release_year),
                        return_years = unique(na.omit(methowSHm$return_year)))
n_rel <- rowSums(methowSHoa[,c("oa1","oa2plus","mort")])

dev.new(width = 7, height = 7)
# png(filename="phi_tot.png", width=7, height=7, units="in", res=300, type="cairo-png")
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths = c(4,1))

# Annual phi_WNFH-BOA
par(mar = c(4.1,6,2,1))
plot(ry, rep(0.5,length(ry)), pch = "", xlim = c(min(ry), max(ry) + 1), ylim = range(unlist(phi_r2a)), 
     xaxt = "n", yaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.8, cex.main = 1.8,
     xlab = "", ylab = "", main = bquote("WNFH" %->% " BOA"))
mtext(bquote("Survival"), side = 2, line = 4, cex = par("cex")*1.8)
axis(1, at = c(unique(round(ry))), labels = unique(round(ry)), cex.axis = 1.2)
axis(1, at = max(round(ry)) + 1, tick = FALSE, labels = "hyper-\nmean", 
     padj = 0.5, cex.axis = 1.2)
vioplot2(cbind(phi_r2a$S1, rowMeans(phi_r2a$S1)), at = c(ry[seq(1, length(ry), 2)], max(ry) + 0.7), 
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "darkgray", lwd.quantile = c(1,2,1))
vioplot2(cbind(phi_r2a$S2, rowMeans(phi_r2a$S2)), at = c(ry[seq(2, length(ry), 2)], max(ry) + 1), 
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "black", lwd.quantile = c(1,2,1))
points(ry, n/n_rel, pch = 16, cex = 1.2, 
       col = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"))
segments(ry, y0 = binconf(n, n_rel, alpha = 0.1)[,"Lower"], 
         y1 = binconf(n, n_rel, alpha = 0.1)[,"Upper"],
         col = ifelse(methowSHoa$smolt_age=="S1", "darkgray", "black"))
abline(v = max(ry) + 0.35, lty = 3)
legend("top", c("1","2"), title = "Smolt age", horiz = TRUE, pch = 16, cex = 1.2, lty = 1,
       col = c("darkgray","black"))
# Hyper-mean logit-scale smolt age effect on phi_WNFH-BOA
par(mar = c(4.1,0,2,5))
plot(1, 0, pch = "", ylim = range(qlogis(rowMeans(phi_r2a$S2)) - qlogis(rowMeans(phi_r2a$S1))),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4, cex.axis = 1.2, las = 1)
axis(side = 1, at = 1, tick = FALSE, labels = "hyper-\nmean", padj = 0.5, cex.axis = 1.2)
text(par("usr")[2] + 1, mean(par("usr")[3:4]), labels = "Smolt age effect", 
     srt = -90, xpd = NA, cex = 1.8)
abline(h = 0, lty = 2)
vioplot2(qlogis(rowMeans(phi_r2a$S2)) - qlogis(rowMeans(phi_r2a$S1)), at = 1,
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "black", lwd.quantile = c(1,2,1))
# Annual phi_tot
par(mar = c(4.1,6,2,1))
plot(ry, rep(0.5,length(ry)), pch = "", xlim = c(min(ry), max(ry) + 1), ylim = range(unlist(phi_tot)), 
     xaxt = "n", yaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.8, cex.main = 1.8,
     xlab = "Release year", ylab = "", main = bquote("WNFH" %->% " LMR"))
mtext(bquote("Survival"), side = 2, line = 4, cex = par("cex")*1.8)
axis(1, at = c(unique(round(ry))), labels = unique(round(ry)), cex.axis = 1.2)
axis(1, at = max(round(ry)) + 1, tick = FALSE, labels = "hyper-\nmean", 
     padj = 0.5, cex.axis = 1.2)
vioplot2(cbind(phi_tot$S1, rowMeans(phi_tot$S1)), at = c(ry[seq(1, length(ry), 2)], max(ry) + 0.7), 
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "darkgray", lwd.quantile = c(1,2,1))
vioplot2(cbind(phi_tot$S2, rowMeans(phi_tot$S2)), at = c(ry[seq(2, length(ry), 2)], max(ry) + 1), 
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "black", lwd.quantile = c(1,2,1))
abline(v = max(ry) + 0.35, lty = 3)
# Hyper-mean logit-scale smolt age effect on phi_tot
par(mar = c(4.1,0,2,5))
plot(1, 0, pch = "", ylim = range(qlogis(rowMeans(phi_tot$S2)) - qlogis(rowMeans(phi_tot$S1))),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4, cex.axis = 1.2, las = 1)
axis(side = 1, at = 1, tick = FALSE, labels = "hyper-\nmean", padj = 0.5, cex.axis = 1.2)
text(par("usr")[2] + 1, mean(par("usr")[3:4]), labels = "Smolt age effect", 
     srt = -90, xpd = NA, cex = 1.8)
abline(h = 0, lty = 2)
vioplot2(qlogis(rowMeans(phi_tot$S2)) - qlogis(rowMeans(phi_tot$S1)), at = 1,
         quantiles = c(0.05,0.5,0.95), col = NA, wex = 0.3, drawRect = FALSE, add = TRUE, 
         border = "black", lwd.quantile = c(1,2,1))
# dev.off()








