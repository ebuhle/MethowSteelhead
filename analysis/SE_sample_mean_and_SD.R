# "Power analysis" for precision of mean and SD of normal distribution
# applied to hatchery smolt lengths

options(device = windows)

# Data
len <- unlist(read.csv("MethowSPCH-Lengths2017.csv", header = T))
log_len <- log(len)

# Confirm length distribution is close to lognormal
# windows()
png(filename="qq_log_length.png", width=7, height=7, units="in", res=300, type="cairo-png")
qqnorm(log_len, main = "Normal QQ plot of log(length)", cex.axis = 1.2, cex.lab = 1.5, las = 1)
qqline(log_len)
dev.off()

# Estimate "true" log-mean and log-SD from large sample
mu <- mean(log_len)
sigma <- sd(log_len)

# Functions to calculate SE of sample mean and SD given their true values and sample size
# The latter is from https://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf
# and is valid for N > 10
se_mean <- function(sigma, N) sigma/sqrt(N)
se_sd <- function(sigma, N) sigma/sqrt(2*(N - 1))
  
# Plot results
# windows(width = 7, height = 10)
png(filename="SE_mu_sigma_log_length.png", width=7, height=10, units="in", res=300, type="cairo-png")
par(mfrow=c(2,1), mar = c(5.1,2.1,0.1,1.1), oma = c(0,3.5,0,0))
curve(se_mean(sigma, x), from = 10, to = 2000, las = 1, xlab = "", ylab = "",
      cex.axis = 1.2, cex.lab = 1.5)
mtext("SE of log-mean", side = 2, line = 4, cex = 1.5)
curve(se_sd(sigma, x), from = 10, to = 2000, las = 1, xlab = "Sample size", ylab = "",
      cex.axis = 1.2, cex.lab = 1.5)
mtext("SE of log-SD", side = 2, line = 4, cex = 1.5)
dev.off()
