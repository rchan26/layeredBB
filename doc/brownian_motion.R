## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(layeredBB)

## ---- fig.align = 'center', fig.dim=c(6,4)------------------------------------
set.seed(2021)
bm_paths <- lapply(1:3, function(i) Brownian_motion_path_sampler(0, seq(0, 1, 0.01)))
# plotting the simulated paths
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

## -----------------------------------------------------------------------------
x <- 0
s <- 0
t <- 2
replicates <- 10000
test <- lapply(1:replicates, function(i) Brownian_motion_path_sampler(x, seq(s, t, 0.01)))

## ---- fig.align = 'center', fig.dim=c(6,4)------------------------------------
plot(x = test[[1]]['time',], y = test[[1]]['X',], type = 'l',
     ylim = c(-4*sqrt(t-s), 4*sqrt(t-s)),
     xlab = '', ylab = '', lwd = 0.5, xaxt = 'n', yaxt = 'n')
for (i in 2:replicates) {
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

## ---- fig.align = 'center', fig.dim=c(6,5)------------------------------------
# select the points at the end time t
end_points <- sapply(1:replicates, function(i) test[[i]]['X', which(seq(s, t, 0.01)==t)])
# plot distribution of the simulated points and the theoretical distribution
curve(dnorm(x, 0, sqrt(t-s)), -4*sqrt(t-s), 4*sqrt(t-s), lwd = 2,
      ylab = 'density')
lines(density(end_points), col = 'red', lty = 2, lwd = 2)
legend('topleft',
       legend = c('theoretical density', 'simulated density'),
       col = c('black', 'red'),
       lty = c(1, 2),
       lwd = c(2, 2))

## ---- fig.align = 'center', fig.dim=c(6,4)------------------------------------
set.seed(2021)
bb_paths <- lapply(1:3, function(i) {
  Brownian_bridge_path_sampler(0, 0, 0, 1, seq(0, 1, 0.01))$full_path})
# plotting the simulated paths
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

## -----------------------------------------------------------------------------
x <- 0
y <- 0
s <- 0
t <- 1
q <- 0.5
replicates <- 10000
paths <- lapply(1:replicates, function(i) {
  Brownian_bridge_path_sampler(x = x,
                               y = y,
                               s = s,
                               t = t,
                               times = seq(s, t, 0.01))$full_path})

## ---- fig.align = 'center', fig.dim=c(6,5)------------------------------------
plot(x = paths[[1]]['time',], y = paths[[1]]['X',], type = 'l',
     ylim = c(-4*sqrt(t-s), 4*sqrt(t-s)),
     xlab = '', ylab = '', lwd = 0.5, xaxt = 'n', yaxt = 'n')
for (i in 2:replicates) {
  lines(x = paths[[i]]['time',], y = paths[[i]]['X',], lwd = 0.5)  
}
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.2), labels=seq(s, t, 0.2), font = 2, cex = 1.5)
axis(1, at=seq(s, t, 0.1), labels=rep("", length(seq(s, t, 0.1))), lwd.ticks = 0.5)
axis(2, at=seq(floor(min(sapply(paths, function(path) path['X',]))),
               ceiling(max(sapply(paths, function(path) path['X',]))), 1),
     labels=seq(floor(min(sapply(paths, function(path) path['X',]))),
                ceiling(max(sapply(paths, function(path) path['X',]))), 1),
     font = 2, cex = 1.5)

## ---- fig.align = 'center', fig.dim=c(6,5)------------------------------------
# select the points at the specified time q
index <- which(seq(s, t, 0.01)==q)
simulated_points <- sapply(1:replicates, function(i) paths[[i]]['X', index])
# calculate the theoretical mean and standard deviation of the simulated points at time q
theoretical_mean <- x + (q-s)*(y-x)/(t-s)
theoretical_sd <- sqrt((t-q)*(q-s)/(t-s))
# plot distribution of the simulated points and the theoretical distribution
curve(dnorm(x, theoretical_mean, theoretical_sd),
      theoretical_mean-4*theoretical_sd,
      theoretical_mean+4*theoretical_sd, lwd = 2,
      ylab = 'density')
lines(density(simulated_points), col = 'red', lty = 2, lwd = 2)
legend('topleft',
       legend = c('theoretical density', 'simulated density'),
       col = c('black', 'red'),
       lty = c(1, 2),
       lwd = c(2, 2))

## ---- fig.align = 'center', fig.dim=c(6,5)------------------------------------
set.seed(2021)
replicates <- 10000
x <- 0
y <- 0
s <- 0
t <- 1
min_low_bound <- min(x,y)-1000
min_up_bound <- min(x,y)
max_low_bound <- max(x,y)
max_up_bound <- max(x,y)+1000
sim_min <- sapply(1:replicates, function(i) min_sampler(x = x,
                                                        y = y,
                                                        s = s,
                                                        t = t,
                                                        low_bound = min_low_bound,
                                                        up_bound = min_up_bound))
sim_max <- sapply(1:replicates, function(i) max_sampler(x = x,
                                                        y = y,
                                                        s = s,
                                                        t = t,
                                                        low_bound = max_low_bound,
                                                        up_bound = max_up_bound))
# plotting minima and maxima
plot(x = sim_min[2,], y = sim_min[1,], ylim = c(-2.5, 2.5),
     pch = 20, lwd = 0.5, cex = 0.25,
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
points(x = sim_max[2,], y = sim_max[1,], pch = 20, lwd = 0.5, cex = 0.25)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)

## ---- fig.align = 'center', fig.dim=c(6,5)------------------------------------
set.seed(2021)
replicates <- 10000
x <- 0
y <- 0
s <- 0
t <- 1
min_low_bound <- min(x,y)-1.5
min_up_bound <- min(x,y)-0.5
max_low_bound <- max(x,y)+1
max_up_bound <- max(x,y)+1.5
sim_min_rest <- sapply(1:replicates, function(i) min_sampler(x = x,
                                                        y = y,
                                                        s = s,
                                                        t = t,
                                                        low_bound = min_low_bound,
                                                        up_bound = min_up_bound))
sim_max_rest <- sapply(1:replicates, function(i) max_sampler(x = x,
                                                        y = y,
                                                        s = s,
                                                        t = t,
                                                        low_bound = max_low_bound,
                                                        up_bound = max_up_bound))
# plotting minima and maxima
plot(x = sim_min_rest[2,], y = sim_min_rest[1,], ylim = c(-2.5, 2.5),
     pch = 20, lwd = 0.5, cex = 0.25,
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
points(x = sim_max_rest[2,], y = sim_max_rest[1,], pch = 20, lwd = 0.5, cex = 0.25)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)
abline(h=c(min_low_bound, min_up_bound, max_low_bound, max_up_bound), lwd = 2)

## ---- fig.align = 'center', fig.dim=c(6,4)------------------------------------
set.seed(2021)
x <- 0
y <- 0
s <- 0
t <- 1
min <- c('min' = -2, 'tau' = 0.5)
besb_paths <- lapply(1:3, function(i) {
  return(min_Bessel_bridge_path_sampler(x = x,
                                        y = y,
                                        s = s,
                                        t = t,
                                        m = min['min'],
                                        tau = min['tau'],
                                        times = seq(s, t, 0.01))$full_path)
})
# plotting simulated paths
plot(x = besb_paths[[1]]['time',], y = besb_paths[[1]]['X',], type = 'l', ylim = c(-3, 3),
     xlab = '', ylab = '', lwd = 3, xaxt = 'n', yaxt = 'n')
lines(x = besb_paths[[2]]['time',], y = besb_paths[[2]]['X',], lwd = 3, lty = 2)
lines(x = besb_paths[[3]]['time',], y = besb_paths[[3]]['X',], lwd = 3, lty = 3)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)
points(x = c(s,t), y = c(x,y), pch = 19, cex = 1)
points(x = min['tau'], y = min['min'], pch = 19, cex = 1.5)

## ---- fig.align = 'center', fig.dim=c(6,4)------------------------------------
set.seed(2021)
x <- 0
y <- 0
s <- 0
t <- 1
max <- c('max' = 2.5, 'tau' = 0.75)
besb_paths <- lapply(1:3, function(i) {
  return(max_Bessel_bridge_path_sampler(x = x,
                                        y = y,
                                        s = s,
                                        t = t,
                                        m = max['max'],
                                        tau = max['tau'],
                                        times = seq(s, t, 0.01))$full_path)
})
# plotting simulated paths
plot(x = besb_paths[[1]]['time',], y = besb_paths[[1]]['X',], type = 'l', ylim = c(-3, 3),
     xlab = '', ylab = '', lwd = 3, xaxt = 'n', yaxt = 'n')
lines(x = besb_paths[[2]]['time',], y = besb_paths[[2]]['X',], lwd = 3, lty = 2)
lines(x = besb_paths[[3]]['time',], y = besb_paths[[3]]['X',], lwd = 3, lty = 3)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)
points(x = c(s,t), y = c(x,y), pch = 19, cex = 1)
points(x = max['tau'], y = max['max'], pch = 19, cex = 1.5)

