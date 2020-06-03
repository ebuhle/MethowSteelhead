library(ggplot2)

J <- 6
NJ <- 200
N <- J*NJ
group <- rep(1:J, each = NJ)
X <- cbind(1, rep(rep(0:1, each = NJ/2), J))

L0 <- rlnorm(N, log(ifelse(X[,2], 10, 15)), 0.5)
dt <- rep(1,N)
beta_log_r <- log(c(5,10))
sigma_log_r <-  0.2
r <- exp(X %*% beta_log_r + rnorm(J, 0, sigma_log_r)[group])
beta_log_q <- c(-0.2,-0.1)
sigma_log_q <- 0.2
q <- exp(X %*% beta_log_q + rnorm(J, 0, sigma_log_q)[group])
sigma <- 0.1

Lt <- rlnorm(N, q * log(L0^(1/q) + r * dt), sigma)

windows()
ggplot(data.frame(L0 = L0, Lt = Lt, X = factor(X[,2]), group = group), 
       aes(x = L0, y = Lt, shape = X, color = X)) +
  geom_point(size = 2) + scale_shape_manual(values = c(1,16)) + 
  scale_color_manual(values = c("darkgray","black")) + labs(x = "L0", y = "Lt") + 
  theme(axis.title = element_text(size = rel(4)), axis.ticks = element_text(size = rel(3)),
        legend.text = element_text(size = rel(4)), strip.text = element_text(size = rel(10))) +
  theme_bw() + facet_wrap(~ group)

# xyplot(Lt ~ L0 | X[,2] + group)


# Fit model
pg_test <- stan(file = here("analysis","stan","parabolic-growth-rXRE-qXRE-NA.stan"),
                data = list(N = N, group = group, K = 2, X = X, L0 = L0, Lt = Lt, dt = dt),
                init = function() { list(beta_log_q = rnorm(2,0,0.5), log_q_z = rep(0,J)) },
                chains = 3, cores = 3, iter = 2000, warmup = 1000,
                control = list(adapt_delta = 0.95, max_treedepth = 12))

print(pg_test)
