library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

# 外源项
mu_term <- "sinusoidal_linear"
mu_fn   <- mu_sinusoidal_linear

num_procs <- 10

plot_n <- 60

fig_out <- "figures/fig_5_cumulative_counts_2016_doublePV(300).pdf"

double_ray_kernel_fun <- function(t, parameters){
  t1 <- t - parameters$delay1
  t2 <- t - parameters$delay2
  k1 <- ifelse(t1 > 0, parameters$alpha1 * t1 * exp(-0.5 * parameters$delta1 * t1^2), 0)
  k2 <- ifelse(t2 > 0, parameters$alpha2 * t2 * exp(-0.5 * parameters$delta2 * t2^2), 0)
  k1 + k2
}

real_events       <- readRDS("data/pv_onset_2016.rds")
real_event_times  <- real_events$time
imported_events   <- readRDS("output/pv_imported_2016.rds")
T_max             <- max(real_event_times)

# 双峰拟合参数
parameters <- readRDS("output/pv_double_fit_2016.rds")
parameters$delay1 <- 15
parameters$delay2 <- 45

real_intensity <- compute_intensity_function(
  events     = real_event_times,
  kernel     = double_ray_kernel_fun,
  parameters = parameters,
  N = 5000, T_max = T_max, mu_fn = mu_fn
)

parameters_mu <- parameters
parameters_mu$alpha1 <- 0
parameters_mu$alpha2 <- 0

real_intensity_imported <- compute_intensity_function(
  events     = imported_events,
  kernel     = double_ray_kernel_fun,
  parameters = parameters_mu,
  N = 5000, T_max = T_max, mu_fn = mu_fn
)

real_counts           <- compute_count_function(events = real_event_times, N = 5000, T_max = T_max)
real_counts_imported  <- compute_count_function(events = imported_events,   N = 5000, T_max = T_max)

intensities_list     <- vector("list", length = num_procs)
counts_list          <- vector("list", length = num_procs)
intensities_list_imp <- vector("list", length = num_procs)
counts_list_imp      <- vector("list", length = num_procs)
time_vec_master      <- NULL

for (seed in 1:num_procs){
  f <- sprintf("output/pv_double_sims_%d_2016.Rdata", seed)
  if (!file.exists(f)) stop("Missing simulation file: ", f)
  load(f)
  
  #[N_time x N_runs]
  intensities_list    [[seed]] <- t(intensities)
  counts_list         [[seed]] <- t(counts)
  intensities_list_imp[[seed]] <- t(intensities_imp)
  counts_list_imp     [[seed]] <- t(counts_imp)
  

  if (is.null(time_vec_master)) time_vec_master <- time_vec
}

# 合并为 [N_time x (N_runs * num_procs)]
intensities     <- do.call(cbind, intensities_list)
counts          <- do.call(cbind, counts_list)
intensities_imp <- do.call(cbind, intensities_list_imp)
counts_imp      <- do.call(cbind, counts_list_imp)


to_long <- function(mat, times){
  df <- as.data.frame(mat)
  df$times <- times
  long <- tidyr::gather(df, -times, key = "key", value = "values")
  long$rowid <- seq_len(nrow(long))
  long
}

counts_long     <- to_long(counts,     time_vec_master)
counts_long_imp <- to_long(counts_imp, time_vec_master)

if (!is.null(plot_n) && plot_n > 0 && plot_n < ncol(counts)) {
  keep_cols <- sample(colnames(as.data.frame(counts)), plot_n)
  counts_long     <- subset(counts_long,     key %in% keep_cols)
  counts_long_imp <- subset(counts_long_imp, key %in% keep_cols)
}

p_double <- ggplot(counts_long) +
  geom_line(aes(times, values, group = key), alpha = 0.12) +
  geom_line(data = counts_long_imp, aes(times, values, group = key), alpha = 0.12, colour = "blue") +
  geom_line(data = real_counts,          aes(x = t, y = counts), colour = "red") +
  geom_line(data = real_counts_imported, aes(x = t, y = counts), colour = "green") +
  xlab("Time (days)") + ylab("Case count") + theme_bw() +
  scale_x_continuous(expand = c(0, 0))

ggsave(fig_out, p_double, width = 12, height = 5.9)

#counts / counts_imp 维度都是 [N_time x total_runs]
last_idx         <- nrow(counts)
counts_end       <- counts     [last_idx, ]
importations_end <- counts_imp [last_idx, ]
ratio            <- importations_end / counts_end

cat(sprintf("Prop importations (vivax, double-peak): %.2f%% [95 CI %.2f%% - %.2f%%]\n",
            mean(ratio, na.rm = TRUE) * 100,
            quantile(ratio, probs = 0.025, na.rm = TRUE) * 100,
            quantile(ratio, probs = 0.975, na.rm = TRUE) * 100))

rm(list = ls(all.names = TRUE))
print(gc())

# Pv plots
library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)

mu_term <- "sinusoidal_linear"
mu_fn   <- mu_sinusoidal_linear
num_procs <- 10                 
plot_n   <- 60              
fig_out  <- "figures/fig_5_cumulative_counts_2016_singlePV(300).pdf"

real_events      <- readRDS("data/pv_onset_2016.rds")
real_event_times <- real_events$time
imported_events  <- readRDS("output/pv_imported_2016.rds")
T_max            <- max(real_event_times)

parameters <- readRDS("output/pv_fit_2016.rds")
parameters$delay <- 15

real_intensity <- compute_intensity_function(
  events     = real_event_times,
  kernel     = ray_kernel,
  parameters = parameters,
  N = 5000, T_max = T_max, mu_fn = mu_fn
)

parameters_mu <- parameters
parameters_mu$alpha <- 0

real_intensity_imported <- compute_intensity_function(
  events     = imported_events,
  kernel     = ray_kernel,
  parameters = parameters_mu,
  N = 5000, T_max = T_max, mu_fn = mu_fn
)

real_counts          <- compute_count_function(events = real_event_times, N = 5000, T_max = T_max)
real_counts_imported <- compute_count_function(events = imported_events,   N = 5000, T_max = T_max)

intensities_list     <- vector("list", length = num_procs)
counts_list          <- vector("list", length = num_procs)
intensities_list_imp <- vector("list", length = num_procs)
counts_list_imp      <- vector("list", length = num_procs)
time_vec_master      <- NULL

for (seed in 1:num_procs){
  f <- sprintf("output/pv_sims_%d_2016.Rdata", seed)
  if (!file.exists(f)) stop("Missing simulation file: ", f)
  load(f) 
  
  intensities_list    [[seed]] <- t(intensities)
  counts_list         [[seed]] <- t(counts)
  intensities_list_imp[[seed]] <- t(intensities_imp)
  counts_list_imp     [[seed]] <- t(counts_imp)
  
  if (is.null(time_vec_master)) time_vec_master <- time_vec
}

# 合并为 [N_time x (N_runs * num_procs)]
intensities     <- do.call(cbind, intensities_list)
counts          <- do.call(cbind, counts_list)
intensities_imp <- do.call(cbind, intensities_list_imp)
counts_imp      <- do.call(cbind, counts_list_imp)

to_long <- function(mat, times){
  df <- as.data.frame(mat)
  df$times <- times
  long <- tidyr::gather(df, -times, key = "key", value = "values")
  long$rowid <- seq_len(nrow(long))
  long
}

counts_long     <- to_long(counts,     time_vec_master)
counts_long_imp <- to_long(counts_imp, time_vec_master)

if (!is.null(plot_n) && plot_n > 0 && plot_n < ncol(counts)) {
  keep_cols <- sample(colnames(as.data.frame(counts)), plot_n)
  counts_long     <- subset(counts_long,     key %in% keep_cols)
  counts_long_imp <- subset(counts_long_imp, key %in% keep_cols)
}


p2 <- ggplot(counts_long) +
  geom_line(aes(times, values, group = key), alpha = 0.12) +
  geom_line(data = counts_long_imp, aes(times, values, group = key), alpha = 0.12, colour = "blue") +
  geom_line(data = real_counts,          aes(x = t, y = counts), colour = "red") +
  geom_line(data = real_counts_imported, aes(x = t, y = counts), colour = "green") +
  xlab("Time (days)") + ylab("Case count") + theme_bw() +
  scale_x_continuous(expand = c(0, 0))

ggsave(fig_out, p2, width = 12, height = 5.9)

last_idx         <- nrow(counts)
counts_end       <- counts     [last_idx, ]
importations_end <- counts_imp [last_idx, ]
ratio <- importations_end / counts_end

cat(sprintf("Prop importations vivax (single): %.2f%% [95 CI %.2f%% - %.2f%%]\n",
            mean(ratio, na.rm = TRUE) * 100,
            quantile(ratio, probs = 0.025, na.rm = TRUE) * 100,
            quantile(ratio, probs = 0.975, na.rm = TRUE) * 100))

