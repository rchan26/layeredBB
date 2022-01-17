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
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-2, 2, 1), labels=seq(-2, 2, 1), font = 2, cex = 1.5)
axis(2, at=seq(-2, 2, 0.5), labels=rep("", 9), lwd.ticks = 0.5)

# looking at the variance of the end time (X_{t}~N(0,t))
# set variables
x <- 0
s <- 0
t <- 2
replicates <- 10000
test <- lapply(1:replicates, function(i) Brownian_motion_path_sampler(x, seq(s, t, 0.01)))
plot(x = test[[1]]['time',], y = test[[1]]['X',], type = 'l', ylim = c(-4*sqrt(t-s), 4*sqrt(t-s)),
     xlab = '', ylab = '', lwd = 0.5, xaxt = 'n', yaxt = 'n')
for (i in 2:replicates) {
  # test[[i]] <- Brownian_motion_path_sampler(x, seq(s, t, 0.01))
  lines(x = test[[i]]['time',], y = test[[i]]['X',], lwd = 0.5)  
}
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.2), labels=seq(s, t, 0.2), font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.1), labels=rep("", length(seq(s, t, 0.1))), lwd.ticks = 0.5)
axis(2, at=seq(floor(min(sapply(test, function(path) path['X',]))),
               ceiling(max(sapply(test, function(path) path['X',]))), 1),
     labels=seq(floor(min(sapply(test, function(path) path['X',]))),
                ceiling(max(sapply(test, function(path) path['X',]))), 1),
     font = 2, cex = 1.5)

end_points <- sapply(1:replicates, function(i) test[[i]]['X', which(seq(s, t, 0.01)==t)])
curve(dnorm(x, 0, sqrt(t-s)), -4*sqrt(t-s), 4*sqrt(t-s), lwd = 2)
lines(density(end_points), col = 'red', lty = 2, lwd = 2)
legend('topleft',
       legend = c('theoretical density', 'simulated density'),
       col = c('black', 'red'),
       lty = c(1, 2),
       lwd = c(2, 2))
