library(layeredBB)

##### BROWNIAN BRIDGE #####

set.seed(2021)
bb_paths <- lapply(1:3, function(i) Brownian_bridge_path_sampler(0, 0, 0, 1, seq(0, 1, 0.01))$simulated_path)
plot(x = bb_paths[[1]]['time',], y = bb_paths[[1]]['X',], type = 'l', ylim = c(-2, 2),
     xlab = '', ylab = '', lwd = 3, xaxt = 'n', yaxt = 'n')
lines(x = bb_paths[[2]]['time',], y = bb_paths[[2]]['X',], lwd = 3, lty = 2)
lines(x = bb_paths[[3]]['time',], y = bb_paths[[3]]['X',], lwd = 3, lty = 3)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-2, 2, 1), labels=seq(-2, 2, 1), font = 2, cex = 1.5)
axis(2, at=seq(-2, 2, 0.5), labels=rep("", 9), lwd.ticks = 0.5)
points(x = c(0,1), y = c(0,0), pch = 19, cex = 1)

# looking at the variance of the mid point time
# set variables
x <- 0
y <- 0
s <- 0
t <- 1
q <- 0.5
replicates <- 10000
# repeatedly simulate Brownian bridges
paths <- lapply(1:replicates, function(i) {
  Brownian_bridge_path_sampler(x = x,
                               y = y,
                               s = s,
                               t = t,
                               times = seq(s, t, 0.01))
})
plot(x = paths[[1]]$full_path['time',], y = paths[[1]]$full_path['X',], type = 'l',
     ylim = c(-4*sqrt(t-s), 4*sqrt(t-s)),
     xlab = '', ylab = '', lwd = 0.5, xaxt = 'n', yaxt = 'n')
for (i in 2:replicates) {
  lines(x = paths[[i]]$full_path['time',], y = paths[[i]]$full_path['X',], lwd = 0.5)  
}
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.2), labels=seq(s, t, 0.2), font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.1), labels=rep("", length(seq(s, t, 0.1))), lwd.ticks = 0.5)
axis(2, at=seq(floor(min(sapply(paths, function(path) path$full_path['X',]))),
               ceiling(max(sapply(paths, function(path) path$full_path['X',]))), 1),
     labels=seq(floor(min(sapply(paths, function(path) path$full_path['X',]))),
                ceiling(max(sapply(paths, function(path) path$full_path['X',]))), 1),
     font = 2, cex = 1.5)
# select the points at the specified time q
index <- which(seq(s, t, 0.01)==q)
simulated_points <- sapply(1:replicates, function(i) paths[[i]]$full_path['X', index])
# calculate the theoretical mean and standard deviation of the simulated points at time q
theoretical_mean <- x + (q-s)*(y-x)/(t-s)
theoretical_sd <- sqrt((t-q)*(q-s)/(t-s))
# plot distribution of the simulated points and the theoretical distribution
curve(dnorm(x, theoretical_mean, theoretical_sd),
      theoretical_mean-4*theoretical_sd,
      theoretical_mean+4*theoretical_sd, lwd = 2)
lines(density(simulated_points), col = 'red', lty = 2, lwd = 2)
legend('topleft',
       legend = c('theoretical density', 'simulated density'),
       col = c('black', 'red'),
       lty = c(1, 2),
       lwd = c(2, 2))
