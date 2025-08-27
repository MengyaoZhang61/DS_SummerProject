# Loads the library
library(epihawkes)

# Set seed for reproducibility

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
  pv_params  = readRDS("output/pv_fit_2016.rds")
  parameters = (pv_params)
  parameters$delay  = 15
  
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
    events <- hawkes_simulation(events = c(0), kernel = ray_kernel, 
                                T_max = T_max,
                                parameters = parameters, mu_fn = mu_fn,
                                print_level = print_level)
    
    print(sprintf("%d: Number of events: %d", i, length(events)))
    event_times[[i]] <- events
    
    # Compute intensity function
    intensities[i,] <- (compute_intensity_function(events = events, kernel = ray_kernel, 
                                                   parameters = parameters, mu_fn = mu_fn, 
                                                   T_max = T_max, N = N_time))$intensity
    # Compute count function
    counts[i,] <- (compute_count_function(events = events, T_max = T_max, N = N_time))$counts
    
    
    events_imp <- c(0)
    events_imp <- hawkes_simulation(events = c(0), kernel = ray_kernel, 
                                    T_max = T_max,
                                    parameters = parameters, mu_fn = mu_fn,
                                    mu_t_max  = mu_t_max,
                                    print_level = print_level, imported  =  T)
    
    print(sprintf("Number of importations: %d", length(events_imp)))
    event_times_imp[[i]] <- events_imp
    
    # Compute intensity function
    intensities_imp[i,] <- (compute_intensity_function(events = events_imp, kernel = ray_kernel, 
                                                       parameters = parameters, mu_fn = mu_fn, 
                                                       T_max = T_max, N = N_time))$intensity
    # Compute count function
    counts_imp[i,] <- (compute_count_function(events = events_imp, T_max = T_max, N = N_time))$counts
  }
  
  time_vec <-  seq(0, T_max, length.out = N_time) 
  
  save(time_vec, intensities, counts, event_times, 
       intensities_imp, counts_imp, event_times_imp, 
       file=paste0("output/pv_sims_", seed, "_2016.Rdata"))
}