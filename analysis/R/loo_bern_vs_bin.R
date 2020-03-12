bern <- data.frame(obs = rep(1:2, each=50), y = rbinom(100, size = 1, prob = 0.17))
bin <- data.frame(N = tapply(bern$obs, bern$obs, length), 
                  y = tapply(bern$y, bern$obs, sum))

glm_bern <- stan_glm(y ~ 1, data = bern, family = "binomial",
                     prior_intercept = normal(0,2),
                     chains = 3, iter = 2000, warmup = 1000, cores = 3)
summary(glm_bern)
loo(glm_bern)


glm_bin <- stan_glm(cbind(y, N - y) ~ 1, data = bin, family = "binomial",
                     prior_intercept = normal(0,2),
                     chains = 3, iter = 2000, warmup = 1000, cores = 3)
summary(glm_bin)
loo(glm_bin)


