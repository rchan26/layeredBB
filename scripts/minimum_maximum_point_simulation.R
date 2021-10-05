library(layeredBB)

##### MINIMUM AND MAXIMUM POINT SIMULATION WITHOUT RESTRICTION #####

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
plot(x = sim_min[2,], y = sim_min[1,], ylim = c(-2.5, 2.5), pch = 20, lwd = 0.5, cex = 0.25,
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
points(x = sim_max[2,], y = sim_max[1,], pch = 20, lwd = 0.5, cex = 0.25)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)

##### MINIMUM AND MAXIMUM POINT SIMULATION WITH RESTRICTION #####

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
plot(x = sim_min_rest[2,], y = sim_min_rest[1,], ylim = c(-2.5, 2.5), pch = 20, lwd = 0.5, cex = 0.25,
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
points(x = sim_max_rest[2,], y = sim_max_rest[1,], pch = 20, lwd = 0.5, cex = 0.25)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 1), labels=seq(-3, 3, 1), font = 2, cex = 1.5)
axis(2, at=seq(-3, 3, 0.5), labels=rep("", 13), lwd.ticks = 0.5)
abline(h=c(min_low_bound, min_up_bound, max_low_bound, max_up_bound), lwd = 2)
