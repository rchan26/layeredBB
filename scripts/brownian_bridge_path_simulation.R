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
paths <- list()
# repeatedly simulate Brownian bridge 
for (i in 1:replicates) {
  paths[[i]] <- Brownian_bridge_path_sampler(x = x,
                                             y = y,
                                             s = s,
                                             t = t,
                                             times = seq(s, t, 0.01))
}
# select the points at the specified time q
index <- which(seq(s, t, 0.01)==q)
simulated_points <- sapply(1:replicates, function(i) paths[[i]]$full_path['X', index])
# calculate the theoretical mean and standard deviation of the simulated points at time q
theoretical_mean <- x + (q-s)*(y-x)/(t-s)
theoretical_sd <- sqrt((t-q)*(q-s)/(t-s))
# plot distribution of the simulated points and the theoretical distribution
plot(density(simulated_points))
curve(dnorm(x, theoretical_mean, theoretical_sd), add = T, col = 'red')
