library(layeredBB)

##### BESSEL BRIDGE MINIMUM SIMULATION #####

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

##### BESSEL BRIDGE MAXIMUM SIMULATION #####

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
