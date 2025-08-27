# Loads the library
library(epihawkes)

double_ray_kernel_fun <- function(t, parameters){
  t1 <- t - parameters$delay1
  t2 <- t - parameters$delay2
  k1 <- ifelse(t1 > 0, parameters$alpha1 * t1 * exp(-0.5 * parameters$delta1 * t1^2), 0)
  k2 <- ifelse(t2 > 0, parameters$alpha2 * t2 * exp(-0.5 * parameters$delta2 * t2^2), 0)
  k1 + k2
}

# 单个 Rayleigh 分量在未来的最大值
ray_component_future_max <- function(x, alpha, delta, delay){
  a <- max(0, x - delay)
  peak_at <- 1 / sqrt(delta)
  if (a <= peak_at){
    alpha / sqrt(delta) * exp(-0.5)
  } else {
    alpha * a * exp(-0.5 * delta * a^2)
  }
}

# 在当前 time，下界 previous_event_time 之前的历史事件对未来强度的“核上界”之和
double_ray_future_kernel_max_sum <- function(time, previous_event_time, events, parameters){
  idx <- which(events <= previous_event_time)
  if (length(idx) == 0) return(0)
  
  xvec <- time - events[idx]  # 已经过的时间（相对每个历史事件）
  m1 <- vapply(xvec, ray_component_future_max,
               numeric(1),
               alpha = parameters$alpha1, delta = parameters$delta1, delay = parameters$delay1)
  m2 <- vapply(xvec, ray_component_future_max,
               numeric(1),
               alpha = parameters$alpha2, delta = parameters$delta2, delay = parameters$delay2)
  sum(m1 + m2)
}

# 给定时刻的实际强度 λ(t) = μ(t) + sum_g g(t - t_i)
lambda_at_double_ray <- function(time, previous_event_time, events, parameters, mu_fn){
  mu_t <- mu_fn(time, parameters = parameters)
  if (is.null(mu_t)) mu_t <- 0
  idx <- which(events <= previous_event_time)
  if (length(idx) == 0) return(mu_t)
  contrib <- vapply(time - events[idx], double_ray_kernel_fun, numeric(1), parameters = parameters)
  mu_t + sum(contrib)
}

# 采样下一事件Ogata thinning
simulate_next_event_double_ray <- function(time, events, T_max, num_children = 1,
                                           parameters,
                                           mu_fn = mu_none,
                                           mu_t_max = NULL,
                                           imported = FALSE,
                                           print_level = 1){
  i <- 1
  previous_event_time <- time
  ts <- c()
  
  while (i < num_children + 1){
    if (!imported){
      # 全强度：μ 的未来上界 + 双峰核的未来上界
      mu_upper <- if (is.null(mu_t_max)) mu_fn(time, parameters = parameters) else mu_t_max(time = time, T_max = T_max, parameters = parameters)
      if (is.null(mu_upper)) mu_upper <- 0
      lambda_max <- mu_upper +
        double_ray_future_kernel_max_sum(time = time,
                                         previous_event_time = previous_event_time,
                                         events = events,
                                         parameters = parameters)
    } else {
      # 只模拟 importations：只需 μ 的未来上界
      if (is.null(mu_t_max))
        stop("For imported = TRUE, please provide mu_t_max(time, T_max, parameters).")
      lambda_max <- mu_t_max(time = time, T_max = T_max, parameters = parameters)
    }
    
    if (lambda_max <= 0 || !is.finite(lambda_max)){
      stop("lambda_max is non-positive or non-finite. Check parameters and mu_t_max.")
    }
    
    # 拟议时间
    u <- stats::runif(1)
    tau <- -log(u) / lambda_max
    time <- time + tau
    if (time >= T_max) return(ts)
    
    # 接受率
    if (!imported){
      lambda <- lambda_at_double_ray(time = time,
                                     previous_event_time = previous_event_time,
                                     events = events,
                                     parameters = parameters, mu_fn = mu_fn)
    } else {
      lambda <- mu_fn(time, parameters = parameters); if (is.null(lambda)) lambda <- 0
    }
    s <- stats::runif(1)
    if (s <= lambda / lambda_max){
      ts <- c(ts, time)
      i <- i + 1
    }
  }
  ts
}

# 主循环
hawkes_simulation_double_ray <- function(events, T_max = Inf, N_max = Inf,
                                         parameters,
                                         mu_fn = mu_none, mu_t_max = NULL,
                                         imported = FALSE, print_level = 1){
  # 起始时间 = 第一事件
  current_time <- events[1]
  while (TRUE){
    new_times <- simulate_next_event_double_ray(time = current_time, events = events,
                                                T_max = T_max, num_children = 1,
                                                parameters = parameters,
                                                mu_fn = mu_fn, mu_t_max = mu_t_max,
                                                imported = imported, print_level = print_level)
    if (length(new_times) == 0) break

    current_time <- new_times[length(new_times)]
    events <- sort(c(events, current_time))
    if (current_time > T_max || length(events) > N_max || is.infinite(current_time)) break
  }
  events
}


for (seed in 1:10){
  set.seed(seed)
  
  mu_t_max_sinusoidal_linear <- function(time, T_max, parameters){
    roots <- all_roots(mu_diff_sinusoidal_linear, interval = c(time, T_max), parameters = parameters)
    mu_values <- mu_sinusoidal_linear(roots, parameters = parameters)
    
    if(length(roots) == 0){
      max_value <- mu_sinusoidal_linear(time, parameters = parameters)
    } else {
      max_value <- max(mu_values)
    }
    return(max_value)
  }
  
  #-------------------------------------------------------------------------------------------------------
  pv_params  = readRDS("output/pv_double_fit_2016.rds")
  parameters = (pv_params)
  parameters$delay1  = 15
  parameters$delay2  = 45
  
  T_max <- 1820
  print_level <- 1
  
  N_runs  =  300
  N_time <- 5001
  
  mu_term <- "sinusoidal_linear"
  mu_fn <- mu_sinusoidal_linear
  mu_t_max <- mu_t_max_sinusoidal_linear
  
  intensities <- matrix(nrow = N_runs, ncol = N_time)
  counts <- matrix(nrow = N_runs, ncol = N_time)
  event_times <- vector(length = N_runs, "list")
  
  intensities_imp <- matrix(nrow = N_runs, ncol = N_time)
  counts_imp <- matrix(nrow = N_runs, ncol = N_time)
  event_times_imp <- vector(length = N_runs, "list")
  
  for (i in 1:N_runs){
    events <- c(0)
    events <- hawkes_simulation_double_ray(
      events = c(0), T_max = T_max,
      parameters = parameters, mu_fn = mu_fn,
      mu_t_max = mu_t_max,
      imported = FALSE, print_level = print_level
    )
    
    print(sprintf("%d: Number of events: %d", i, length(events)))
    event_times[[i]] <- events
    
    # Compute intensity function
    intensities[i, ] <- compute_intensity_function(
      events = events, kernel = double_ray_kernel_fun,
      parameters = parameters, mu_fn = mu_fn, T_max = T_max, N = N_time
    )$intensity
    
    # Compute count function
    counts[i,] <- (compute_count_function(events = events, T_max = T_max, N = N_time))$counts
    
    
    events_imp <- c(0)
    events_imp <- hawkes_simulation_double_ray(
      events = c(0), T_max = T_max,
      parameters = parameters, mu_fn = mu_fn,
      mu_t_max = mu_t_max,
      imported = TRUE, print_level = print_level
    )
    
    print(sprintf("Number of importations: %d", length(events_imp)))
    event_times_imp[[i]] <- events_imp
    
    # Compute intensity function
    intensities[i, ] <- compute_intensity_function(
      events = events, kernel = double_ray_kernel_fun,
      parameters = parameters, mu_fn = mu_fn, T_max = T_max, N = N_time
    )$intensity
    
    # Compute count function
    counts_imp[i,] <- (compute_count_function(events = events_imp, T_max = T_max, N = N_time))$counts
  }
  
  time_vec <-  seq(0, T_max, length.out = N_time) 
  
  save(time_vec, intensities, counts, event_times, 
       intensities_imp, counts_imp, event_times_imp, 
       file=paste0("output/pv_double_sims_", seed, "_2016.Rdata"))
}