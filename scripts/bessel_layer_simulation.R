library(layeredBB)

##### BESSEL LAYER SIMULATION #####

set.seed(2021)
x <- -0.4
y <- 0.2
s <- 0
t <- 1
bes_layer <- bessel_layer_simulation(x = x, y = y, s = s, t = t, mult = 0.2)

plot(x = c(s,t), y = c(x,y), ylim = c(-1, 1), pch = 20, lwd = 4,
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
mtext('Time', 1, 2.75, font = 2, cex = 1.5)
mtext('W', 2, 2.75, font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.2), labels=c("0.0", seq(0.2, 0.8, 0.2), "1.0"), font = 2, cex = 1.5)
axis(1, at=seq(0, 1, 0.1), labels=rep("", 11), lwd.ticks = 0.5)
axis(2, at=seq(-1, 1, 0.5), labels=seq(-1, 1, 0.5), font = 2, cex = 1.5)
axis(2, at=seq(-1, 1, 0.1), labels=rep("", 21), lwd.ticks = 0.5)
abline(h=c(bes_layer$L, bes_layer$U), lwd = 3)
abline(h=c(bes_layer$l, bes_layer$u), lwd = 3, lty = 2)
