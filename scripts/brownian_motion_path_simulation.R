library(layeredBB)

##### BROWNIAN MOTION #####

set.seed(2021)
bm_paths <- lapply(1:3, function(i) Brownian_motion_path_sampler(0, seq(0, 1, 0.01)))
plot(x = bm_paths[[1]]['time',], y = bm_paths[[1]]['X',], type = 'l', ylim = c(-2, 2),
     xlab = '', ylab = '', lwd = 3, xaxt = 'n', yaxt = 'n')
lines(x = bm_paths[[2]]['time',], y = bm_paths[[2]]['X',], lwd = 3, lty = 2)
lines(x = bm_paths[[3]]['time',], y = bm_paths[[3]]['X',], lwd = 3, lty = 3)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(-2, 2, 1), labels=seq(-2, 2, 1), font = 2, cex = 1.5)
axis(2, at=seq(-2, 2, 0.5), labels=rep("", 9), lwd.ticks = 0.5)

# looking at the variance of the end time (X_{1}~N(0,1))
# set variables
x <- 0
s <- 0
t <- 2
replicates <- 10000
test <- list()
test[[1]] <- Brownian_motion_path_sampler(x, seq(s, t, 0.01))
plot(x = test[[1]]['time',], y = test[[1]]['X',], type = 'l', ylim = c(-4*sqrt(t-s), 4*sqrt(t-s)),
     xlab = 'Time', ylab = 'X', lwd = 0.5)
for (i in 2:replicates) {
  test[[i]] <- Brownian_motion_path_sampler(x, seq(s, t, 0.01))
  lines(x = test[[i]]['time',], y = test[[i]]['X',], lwd = 0.5)  
}

end_points <- sapply(1:replicates, function(i) test[[i]]['X', which(seq(s, t, 0.01)==t)])
curve(dnorm(x, 0, sqrt(t-s)), -4*sqrt(t-s), 4*sqrt(t-s))
lines(density(end_points), col = 'red')
