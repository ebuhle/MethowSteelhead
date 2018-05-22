Ni <- 10
M <- 10
group <- rep(1:M, each = Ni)
n <- 50
x <- runif(Ni*M, -1, 1)
f <- factor(as.numeric(x > 0))
mu_beta <- c(0,0.5)
sigma_beta <- c(0.2,0.2)
beta <- sapply(1:2, function(i) rnorm(M, mu_beta[i], sigma_beta[i]))
y <- rbinom(Ni*M, size = n, prob = plogis(beta[group,1] + beta[group,2]*x))
dat_test <- data.frame(group = group, x = x, f = f, n = n, y = y)

test <- stan_glmer(cbind(y, n - y) ~ f + (as.numeric(f) || group), 
                   data = dat_test,
                   family = binomial("logit"),
                   prior = normal(0,3),
                   prior_intercept = normal(0,1.5), 
                   prior_covariance = decov(),
                   chains = 4, iter = 2000, warmup = 1000, cores = 4)
summary(test)

