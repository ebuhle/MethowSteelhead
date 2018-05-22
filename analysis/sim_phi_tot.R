sim_phi_tot <- function(cjs, fit_oa, iter = NULL, adult_age, release_years, return_years, 
                        Tmax = NULL, CH = NULL) 
{
  # Actual unique years, not factor levels
  rlsy <- sort(unique(release_years))
  rtny <- sort(unique(return_years))
  # adult_age is the centered covariate used in CJS models
  mean_adult_age <- attributes(adult_age)[["scaled:center"]]
  aa <- 2:4 - mean_adult_age
  names(aa) <- 2:4 # ignore 5-year-olds
  
  # CJS
  TT <- 7  # nr. occasions
  if(is.null(Tmax)) Tmax <- TT - 1      # calculate survival for occasions 1:Tmax
  if(is.null(iter)) iter <- length(as.matrix(cjs,"lp__")) # nr. MCMC iters 
  J <- max(length(rlsy), length(rtny))  # max nr. groups 
  K <- 3                                # nr. covariates
  mbeta <- as.matrix(cjs, "beta")[1:iter,]
  beta <- array(NA, c(iter, TT - 1, K))
  for(k in 1:K)
    beta[,,k] <- mbeta[, paste0("beta[", k, ",", 1:(TT-1), "]")]
  mepsilon_z <- as.matrix(cjs, "epsilon_z")[1:iter,]
  epsilon_z <- array(NA, c(iter, TT - 1, J))
  for(j in 1:J)
    epsilon_z[,,j] <- mepsilon_z[, paste0("epsilon_z[", j, ",", 1:(TT-1), "]")]
  sigma <- as.matrix(cjs, "sigma")[1:iter,]  # iter * TT - 1
  epsilon <- sweep(epsilon_z, 1:2, sigma, "*")  # iter * TT - 1 * J
  b <- as.matrix(cjs, "b")[1:iter, paste0("b[", 1, ",", 1:TT, "]")]  # iter * TT
  me_z <- as.matrix(cjs, "e_z")[1:iter,]
  e_z <- array(NA, c(iter, TT, J))
  for(j in 1:J)
    e_z[,,j] <- me_z[, paste0("e_z[", j, ",", 1:TT, "]")]
  s <- as.matrix(cjs, "s")[1:iter,]  # iter * TT
  e <- sweep(e_z, 1:2, s, "*")

  # OA GLMM
  dat <- data.frame(release_year = rep(rlsy, each = 2), smolt_age = c("S1","S2"), smolt_age_num = 1:2)
  prob_oa2 <- posterior_linpred(fit_oa, transform = TRUE, newdata = dat)
  prob_oa2 <- prob_oa2[1:iter,]  # GLMM iter > CJS iter
                             
  # Calculate cumulative survival (Gilbert & Rich age notation)
  phi_1.1 <- array(NA, c(iter, TT - 1, length(rlsy)))
  dimnames(phi_1.1)[[3]] <- rlsy
  phi_1.2 <- phi_1.1
  phi_2.1 <- phi_1.1
  phi_2.2 <- phi_1.1
  phi_tot_S1 <- matrix(NA, iter, length(rlsy))
  colnames(phi_tot_S1) <- rlsy
  phi_tot_S2 <- phi_tot_S1
  p_x.1 <- array(NA, c(iter, TT, length(rlsy)))
  dimnames(p_x.1)[[3]] <- rlsy
  p_x.2 <- p_x.1
  for(j in 1:length(rlsy))
  {
    phi_1.1[,,j] <- plogis(beta[,,1] + aa["2"] * beta[,,3] + 
                             cbind(epsilon[,1:3,j], epsilon[,4:6,j]))
    phi_1.2[,,j] <- plogis(beta[,,1] + aa["3"] * beta[,,3] + 
                             cbind(epsilon[,1:3,j], epsilon[,4:6,j+1]))
    prob_1.2 <- prob_oa2[, dat$release_year==rlsy[j] & dat$smolt_age=="S1"]
    phi_tot_S1[,j] <- (1 - prob_1.2) * apply(phi_1.1[,1:Tmax,j], 1, prod) + 
      prob_1.2 * apply(phi_1.2[,1:Tmax,j], 1, prod)
    phi_2.1[,,j] <- plogis(beta[,,1] + beta[,,2] + aa["3"] * beta[,,3] + 
                             cbind(epsilon[,1:3,j], epsilon[,4:6,j]))
    phi_2.2[,,j] <- plogis(beta[,,1] + beta[,,2] + aa["4"] * beta[,,3] + 
                             cbind(epsilon[,1:3,j], epsilon[,4:6,j+1]))
    prob_2.2 <- prob_oa2[, dat$release_year==rlsy[j] & dat$smolt_age=="S2"]
    phi_tot_S2[,j] <- (1 - prob_2.2) * apply(phi_2.1[,1:Tmax,j], 1, prod) + 
      prob_2.2 * apply(phi_2.2[,1:Tmax,j], 1, prod)
    p_x.1[,,j] <- plogis(b + cbind(e[,1:3,j], e[,4:7,j]))
    p_x.2[,,j] <- plogis(b + cbind(e[,1:3,j], e[,4:7,j+1]))
  }
  
  # Generate replicate data from PPD using initial conditions in CH
  # CH is data frame with ID, release_year, smolt_age, adult_age, return_year 
  # (latter two can be NA), and capture histories for occasions 1:TT (rows are individuals)
  if(!is.null(CH))   
  {
    N <- nrow(CH)
    CH_sim <- CH
    CHm_sim <- matrix(NA, 2^(TT - 1) * length(rlsy) * 2 * 4, iter)  # simulated m-arrays
    s_a <- CH_sim$smolt_age
    r_y <- as.character(CH_sim$release_year)
    for(i in 1:iter)
    {
      o_a <- 1 + rbinom(N, size = 1, 
                        prob = prob_oa2[i, match(paste(r_y, s_a), 
                                                 paste(dat$release_year, dat$smolt_age))])
      CH_sim$adult_age <- ifelse(s_a=="S1", 1, 2) + o_a
      CH_sim$return_year <- CH_sim$release_year + o_a
      phi <- matrix(NA, N, TT - 1)
      phi[s_a=="S1" & o_a==1,] <- t(phi_1.1[i,,r_y[s_a=="S1" & o_a==1]])
      phi[s_a=="S1" & o_a==2,] <- t(phi_1.2[i,,r_y[s_a=="S1" & o_a==2]])
      phi[s_a=="S2" & o_a==1,] <- t(phi_2.1[i,,r_y[s_a=="S2" & o_a==1]])
      phi[s_a=="S2" & o_a==2,] <- t(phi_2.2[i,,r_y[s_a=="S2" & o_a==2]])
      p <- matrix(NA, N, TT)
      p[o_a==1,] <- t(p_x.1[i,,r_y[o_a==1]])
      p[o_a==2,] <- t(p_x.2[i,,r_y[o_a==2]])
      z <- cbind(1, matrix(NA, N, TT - 1))  # latent states: alive or dead
      for(tt in 2:TT)
        z[,tt] <- z[,tt-1] * rbinom(N, size = 1, prob = phi[,tt - 1])
      # observations
      y <- cbind(1, matrix(rbinom(N * (TT - 1), size = 1, 
                                  prob = z[,-1] * p[,-1]), N, TT - 1))
      colnames(y) <- names(CH_sim)[!names(CH_sim) %in% c("tag","release_year","smolt_age",
                                                         "adult_age","return_year")]
      CH_sim$adult_age <- factor(replace(CH_sim$adult_age, CH_sim$TDAWEA==0, NA), exclude = "")
      CHm_sim_i <- aggregate(WNFH ~ .,
                             data = cbind(CH_sim[,c("release_year","smolt_age","adult_age")], y), 
                             length, drop = FALSE)
      CHm_sim[,i] <- CHm_sim_i$WNFH
    }
    CHm_sim <- cbind(CHm_sim_i[,colnames(CHm_sim_i) != "WNFH"], CHm_sim)
    CHm_sim$adult_age <- as.numeric(as.character(CHm_sim$adult_age))
  } else {
    CHm_sim <- NULL
  }
  return(list(S1 = phi_tot_S1, S2 = phi_tot_S2, CHm_sim = CHm_sim))
}



