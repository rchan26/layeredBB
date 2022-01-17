library(layeredBB)

##### CROSSING #####

set.seed(11)
x1 <- -2
y1 <- 4
x2 <- 1.5
y2 <- -0.2
s <- 0
t <- 1
times <- seq(s, t, 0.001)
bridge1 <- Brownian_bridge_path_sampler(x1, y1, s, t, times)
bridge2 <- Brownian_bridge_path_sampler(x2, y2, s, t, times)
plot(x = bridge1$full_path['time',], y = bridge1$full_path['X',], type = 'l',
     ylim = c(-3, 5), lwd = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
lines(x = bridge2$full_path['time',], y = bridge2$full_path['X',], type = 'l',
      lwd = 3)
abline(v=c(s,t))
points(x=c(s,s,t,t), y=c(x1,x2,y1,y2), pch = 20, lwd = 3)
axis(1, at=c(s,t), labels=c('0', 'T'), font = 2, cex = 1.5)
axis(2, at=c(x1,x2), labels=c(expression(bold(x^"(1)")), expression(bold(x^"(2)"))), las = 2)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('X', 2, 2.75, font = 2, cex = 1.5)

##### COALESCING #####

set.seed(10)
C <- 5
x <- runif(C, -5, 5)
y <- rnorm(1, mean(x), sqrt(t/C))
s <- 0
t <- 1
times <- seq(s, t, 0.001)
bridge <- lapply(1:C, function(c) Brownian_bridge_path_sampler(x[c], y, s, t, times))
plot(x = bridge[[1]]$full_path['time',], y = bridge[[1]]$full_path['X',], type = 'l',
     ylim = c(min(x)-1, max(x)+1), lwd = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
for (c in 2:C) {
  lines(x = bridge[[c]]$full_path['time',], y = bridge[[c]]$full_path['X',], type = 'l',
        lwd = 3)
}
abline(v=c(s,t))
points(x=c(rep(s,5), t), y=c(x,y), pch = 20, lwd = 3)
axis(1, at=c(s,t), labels=c('0', 'T'), font = 2, cex = 1.5)
axis(2, at=c(x), labels=c(expression(bold(x^"(1)")),
                          expression(bold(x^"(2)")),
                          expression(bold(x^"(3)")),
                          expression(bold(x^"(4)")),
                          expression(bold(x^"(5)"))), las = 2)
axis(4, at=y, labels=expression(bold(y)), las = 2)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('X', 2, 2.75, font = 2, cex = 1.5)

##### COALESCING (with temporal partition) #####

set.seed(10)
C <- 5
x <- runif(C, -5, 5)
y <- rnorm(1, mean(x), sqrt(t/C))
s <- 0
t <- 1
times <- seq(s, t, 0.001)
temporal_partition <- seq(s, t, (t-s)/10)
bridge <- lapply(1:C, function(c) Brownian_bridge_path_sampler(x[c], y, s, t, times))
plot(x = bridge[[1]]$full_path['time',], y = bridge[[1]]$full_path['X',], type = 'l',
     ylim = c(min(x)-1, max(x)+1), lwd = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
for (c in 2:C) {
  lines(x = bridge[[c]]$full_path['time',], y = bridge[[c]]$full_path['X',], type = 'l',
        lwd = 3)
}
abline(v=c(s,t))
abline(v=temporal_partition, lty = 2)
points(x=c(rep(s,5), t), y=c(x,y), pch = 20, lwd = 3)
axis(1, at=c(s,t), labels=c(expression(t[0]), expression(t[n])), font = 2, cex = 2)
axis(1, at=c(temporal_partition[2:4], temporal_partition[length(temporal_partition)-1]),
     labels=c(expression(t[1]), expression(t[2]), '...', expression(t[n-1])), font = 2, cex = 2, lwd.ticks = 0)
axis(1, at=temporal_partition, labels=rep("", length(temporal_partition)), font = 2, cex = 1, lwd.ticks = 0.5)
axis(2, at=c(x), labels=c(expression(bold(x[0]^"(1)")),
                          expression(bold(x[0]^"(2)")),
                          expression(bold(x[0]^"(3)")),
                          expression(bold(x[0]^"(4)")),
                          expression(bold(x[0]^"(5)"))),
     las = 2)
axis(4, at=y, labels=expression(bold(y)), las = 2)
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('X', 2, 2.75, font = 2, cex = 1.5)
