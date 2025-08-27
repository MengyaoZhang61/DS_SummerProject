# Goodness of fits
library(epihawkes)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(goftest)

# Choose distribution
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear
mu_diff_fn <- mu_diff_sinusoidal_linear
mu_int_fn  <- mu_int_sinusoidal_linear

# Vivax (bimodal)
real_events <- readRDS("data/pv_onset_2016.rds")
real_event_times <- as.numeric(real_events$time)

pv1_par <- as.list(readRDS("output/pv_double_fit_2016.rds"))

delay1 <- 15
delay2 <- 45

double_ray_kernel <- function(t, alpha1, alpha2, delta1, delta2, delay1, delay2) {
  t1 <- t - delay1; t2 <- t - delay2
  k1 <- ifelse(t1 > 0, pv1_par$alpha1 * t1 * exp(-0.5 * pv1_par$delta1 * t1^2), 0)
  k2 <- ifelse(t2 > 0, pv1_par$alpha2 * t2 * exp(-0.5 * pv1_par$delta2 * t2^2), 0)
  k1 + k2
}

# 对每个历史事件j, 到时间t的核积分
int_double_ray_upto_t <- function(t, t_j, alpha1, alpha2, delta1, delta2, delay1, delay2) {
  # 对峰1
  x1 <- t - t_j - delay1
  I1 <- if (x1 > 0) (alpha1/delta1) * (1 - exp(-0.5*delta1*x1^2)) else 0
  # 对峰2
  x2 <- t - t_j - delay2
  I2 <- if (x2 > 0) (alpha2/delta2) * (1 - exp(-0.5*delta2*x2^2)) else 0
  I1 + I2
}

if (is.null(pv1_par$alpha1) || is.null(pv1_par$alpha2) || is.null(pv1_par$delta1) || is.null(pv1_par$delta2)) {

  stopifnot(!is.null(pv1_par$alpha), !is.null(pv1_par$delta))
  pv1_par$alpha1 <- 0.6 * pv1_par$alpha
  pv1_par$alpha2 <- 0.4 * pv1_par$alpha
  pv1_par$delta1 <- 0.8 * pv1_par$delta
  pv1_par$delta2 <- 0.6 * pv1_par$delta
}

alpha1 <- pv1_par$alpha1; alpha2 <- pv1_par$alpha2
delta1 <- pv1_par$delta1; delta2 <- pv1_par$delta2

# 计算补偿过程Λ(t_i)并做 KS/QQ
n_events <- length(real_event_times)
Lambda <- numeric(n_events)
for (i in seq_len(n_events)) {
  ti <- real_event_times[i]
  # μ 积分到 ti
  mu_int <- mu_int_fn(ti, parameters = list(A = pv1_par$A, B = pv1_par$B, M = pv1_par$M, N = pv1_par$N))
  # 核积分到 ti
  if (i > 1) {
    past <- real_event_times[1:(i-1)]
    kern_sum <- 0
    for (tj in past) {
      kern_sum <- kern_sum + int_double_ray_upto_t(ti, tj, alpha1, alpha2, delta1, delta2, delay1, delay2)
    }
  } else {
    kern_sum <- 0
  }
  Lambda[i] <- mu_int + kern_sum
}

taus <- diff(Lambda)
zks  <- 1 - exp(-taus)
sorted_zks <- sort(zks)
sorted_zks <- sorted_zks[sorted_zks > 0]
n <- length(sorted_zks)
k <- 1:n
bk <- (k - 0.5)/n

ks_data_pv1_double <- data.frame(
  zk = sorted_zks,
  bk = bk,
  kernel = rep("double-ray"), 
  n = length(sorted_zks)
)

# KS 图
p_ks_pv1 <- ggplot(ks_data_pv1_double) + 
  geom_point(aes(zk, bk)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_abline(intercept = 0 + 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  geom_abline(intercept = 0 - 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Cumulative Distribution Function")

quantiles <- qbeta(sorted_zks, shape1 = k, shape2 = n - k + 1)
qq_data_pv1_double <- data.frame(
  zk = sorted_zks,
  qs = quantiles, 
  li = sorted_zks + 1.96 * sqrt(sorted_zks * (1 - sorted_zks) / n),
  ui = sorted_zks - 1.96 * sqrt(sorted_zks * (1 - sorted_zks) / n),
  kernel = rep("double-ray"), 
  n = length(sorted_zks)
)

p_qq_pv1 <- ggplot(qq_data_pv1_double) + 
  geom_point(aes(zk, qs)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_line(aes(zk, li), col = 'black', linetype = "dashed") +
  geom_line(aes(zk, ui), col = 'black', linetype = "dashed") +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Empirical Quantiles")

# Vivax (Unimodal)
#-------------------------------------------------------------------------------------------------------------
# Loads in real data
real_events <- readRDS(paste0("data/pv_onset_2016.rds"))
parameters <- as.list(readRDS("output/pv_fit_2016.rds"))
parameters$delay <- 15
real_event_times <- as.vector(real_events$time)

# Calculates integral
integral <- vector(length = length(real_event_times))
for (i in 2:length(real_event_times)){
  events_sub <- real_event_times[1:i]
  integral[i] = integral_intensity(events = events_sub, int_kernel = int_ray,
                                   parameters = parameters, mu_fn = mu_sinusoidal_linear,
                                   mu_diff_fn = mu_diff_sinusoidal_linear, 
                                   mu_int_fn = mu_int_sinusoidal_linear,
                                   print_level = 1)
}

taus <- diff(integral)
zks <-  1-exp(-taus)
sorted_zks <- sort(zks)
sorted_zks = sorted_zks[sorted_zks > 0]
n = length(sorted_zks)
k = 1:n
bk = (k-0.5)/n

ks_data_pv_ray <- data.frame(zk = sorted_zks,
                             bk = bk,
                             kernel = rep("ray"), length(sorted_zks))

# KS Plot
p_ks_pv <- ggplot(ks_data_pv_ray) + 
  geom_point(aes(zk, bk)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_abline(intercept = 0 + 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  geom_abline(intercept = 0 - 1.36/n^0.5, slope = 1, col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Cumulative Distribution Function")

quantiles = qbeta(sorted_zks, shape1 = k, shape2 = n-k+1)
qq_data_pv_ray <- data.frame(zk = sorted_zks,
                             qs = quantiles, 
                             li = sorted_zks + 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                             ui = sorted_zks - 1.96*(sorted_zks*(1-sorted_zks)/n)^0.5,
                             kernel = rep("ray"), length(sorted_zks))
p_qq_pv <- ggplot(qq_data_pv_ray ) + 
  geom_point(aes(zk, qs)) + 
  geom_abline(intercept = 0, slope = 1, col = 'black') +
  geom_line(aes(zk, li), col = 'black', linetype = "dashed") +
  geom_line(aes(zk, ui), col = 'black', linetype = "dashed") +
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  xlab("Quantiles") + ylab("Empiracle Quantiles")


#----------------------------------------------------------------------------------------------------
p_ks_pv1 <- p_ks_pv1 + ggtitle("Bimodal - KS Test")
p_ks_pv  <- p_ks_pv  + ggtitle("Unimodal - KS Test")
p_qq_pv1 <- p_qq_pv1 + ggtitle("Bimodal - QQ Plot")
p_qq_pv  <- p_qq_pv  + ggtitle("Unimodal - QQ Plot")


p2 <- ggarrange(p_ks_pv1, p_ks_pv,
                p_qq_pv1, p_qq_pv,                
                labels = "AUTO",
                common.legend = TRUE, 
                legend = "bottom")
ggsave("figures/fig_goodness_2016_pv.pdf", p2, width = 8, height = 8)



# 1) KS 检验
ks_double <- ks.test(ks_data_pv1_double$zk, "punif", 0, 1)
ks_single <- ks.test(ks_data_pv_ray$zk,         "punif", 0, 1)

# 2) Cramér–von Mises检验
cvm_double <- cvm.test(ks_data_pv1_double$zk, null = "punif")
cvm_single <- cvm.test(ks_data_pv_ray$zk,     null = "punif")

# 3) KS 图上的误差
mse_double <- mean((ks_data_pv1_double$zk - ks_data_pv1_double$bk)^2)
mse_single <- mean((ks_data_pv_ray$zk         - ks_data_pv_ray$bk        )^2)

# 4) QQ 图上的误差
rmse_double <- sqrt(mean((qq_data_pv1_double$zk - qq_data_pv1_double$qs)^2))
rmse_single <- sqrt(mean((qq_data_pv_ray$zk      - qq_data_pv_ray$qs     )^2))

comp_tbl <- data.frame(
  Kernel        = c("Double Rayleigh", "Single Rayleigh"),
  KS_D          = c(unname(ks_double$statistic), unname(ks_single$statistic)),
  KS_pvalue     = c(ks_double$p.value,            ks_single$p.value),
  CvM_W2        = c(unname(cvm_double$statistic), unname(cvm_single$statistic)),
  CvM_pvalue    = c(cvm_double$p.value,           cvm_single$p.value),
  KS_MSE        = c(mse_double,                   mse_single),
  QQ_RMSE       = c(rmse_double,                  rmse_single),
  n             = c(ks_data_pv1_double$n[1],      ks_data_pv_ray$kernel %>% length())
)

print(comp_tbl)
write.csv(comp_tbl, "output/gof_numeric_comparison_2016_pv.csv", row.names = FALSE)
