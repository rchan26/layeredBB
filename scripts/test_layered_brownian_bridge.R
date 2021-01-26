# tests to check that:
# 1) the number of points that are returned from layered_brownian_bridge $simulated_path equals the number of points requested
# 2) the simulated Bridge actually lies within the simulated Bessel layer [L,U] 
# test over several number of points to return
seed <- 2021
set.seed(seed)
number_of_points <- c(10, 50, 500, 1000, 2000, 5000, 10000)
for (n in number_of_points) {
  returned_correct_n_of_points <- rep(NA, 1000)
  returned_BB_in_layer <- rep(NA, 1000)
  for (i in 1:1000) {
    x <- runif(1, -1, 1)
    y <- runif(1, -1, 1)
    times <- sort(runif(2, 0, 2))
    mult <- runif(1, 0, 2)
    requested_times <- runif(n, times[1], times[2])
    bes_layer <- bessel_layer_simulation(x = x,
                                         y = y,
                                         s = times[1],
                                         t = times[2],
                                         mult = mult)
    bridge <- layeredBB::layered_brownian_bridge(x = x,
                                                 y = y,
                                                 s = times[1],
                                                 t = times[2],
                                                 bessel_layer = bes_layer,
                                                 times = requested_times)
    # check $simulated_path returns the same number of points as requested (i.e. n = length(requested_times))
    returned_correct_n_of_points[i] = (ncol(bridge$simulated_path) == n)
    # check that the simulated layered BB actually lies within the simulated Bessel layer
    returned_BB_in_layer[i] <- ((min(bridge$simulated_path['X',]) > bes_layer$L) & (max(bridge$simulated_path['X',]) < bes_layer$U))
  }
  print('##########')
  print(paste('n =', n, '|| all returned correct number of points:', all(returned_correct_n_of_points)))
  print(paste('n =', n, '|| all returned a BB within the bessel layer:', all(returned_BB_in_layer)))
}
