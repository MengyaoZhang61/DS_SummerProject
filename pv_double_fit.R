library(epihawkes)
library(ggplot2)
library(optimx)

data <- readRDS("data/pv_onset_2016.rds")
events <- as.numeric(data$time)
T_max <- max(events)
set.seed(1)

# 外源项函数
mu_fn <- mu_sinusoidal_linear
mu_diff_fn <- mu_diff_sinusoidal_linear
mu_int_fn <- mu_int_sinusoidal_linear

# 双峰Rayleigh核函数
double_ray_kernel <- function(t, alpha1, alpha2, delta1, delta2, delay1, delay2) {
  t1 <- t - delay1
  t2 <- t - delay2
  k1 <- ifelse(t1 > 0, alpha1 * t1 * exp(-0.5 * delta1 * t1^2), 0)
  k2 <- ifelse(t2 > 0, alpha2 * t2 * exp(-0.5 * delta2 * t2^2), 0)
  return(k1 + k2)
}

# 计算条件强度
compute_conditional_intensity_direct <- function(events, events_history, alpha1, alpha2, delta1, delta2, delay1, delay2) {
  n <- length(events)
  conditional_intensities <- numeric(n)
  
  for (i in 1:n) {
    current_time <- events[i]
    past_events <- events_history[events_history < current_time]
    
    if (length(past_events) > 0) {
      kernel_values <- double_ray_kernel(current_time - past_events, alpha1, alpha2, delta1, delta2, delay1, delay2)
      conditional_intensities[i] <- sum(kernel_values)
    } else {
      conditional_intensities[i] <- 0
    }
  }
  
  return(conditional_intensities)
}

# 计算核函数积分
int_double_ray_direct <- function(T_max, events, alpha1, alpha2, delta1, delta2, delay1, delay2) {
  integrals <- numeric(length(events))
  
  for (i in 1:length(events)) {
    t_event <- events[i]
    max_time <- T_max - t_event
    
    if (max_time > 0) {
      # 第一个峰的积分
      if (max_time > delay1) {
        t1_max <- max_time - delay1
        integral1 <- (alpha1 / delta1) * (1 - exp(-0.5 * delta1 * t1_max^2))
      } else {
        integral1 <- 0
      }
      
      # 第二个峰的积分
      if (max_time > delay2) {
        t2_max <- max_time - delay2
        integral2 <- (alpha2 / delta2) * (1 - exp(-0.5 * delta2 * t2_max^2))
      } else {
        integral2 <- 0
      }
      
      integrals[i] <- integral1 + integral2
    } else {
      integrals[i] <- 0
    }
  }
  
  return(integrals)
}

# 负对数似然函数
neg_log_likelihood_clean <- function(par) {
  if (length(par) != 8) {
    return(1e10)
  }
  
  alpha1 <- par[1]
  alpha2 <- par[2]
  delta1 <- par[3]
  delta2 <- par[4]
  A <- par[5]
  B <- par[6]
  M <- par[7]
  N <- par[8]
  
  # 参数约束
  if (alpha1 <= 0 || alpha2 <= 0 || delta1 <= 0 || delta2 <= 0) {
    return(1e10)
  }
  
  # 创建参数列表
  parameters <- list(A = A, B = B, M = M, N = N)
  
  # 计算外源项
  mu_ts <- mu_fn(events, parameters = parameters)
  if (any(mu_ts <= 0)) {
    return(1e10)
  }
  
  # 计算条件强度
  events_history <- c(0, events)
  conditional_intensities <- compute_conditional_intensity_direct(
    events, events_history, alpha1, alpha2, delta1, delta2, delay1, delay2
  )
  
  # 总强度
  lambdas <- mu_ts + conditional_intensities
  if (any(lambdas <= 0)) {
    return(1e10)
  }
  
  # 对数似然的第一部分
  sum_log_lambdas <- sum(log(lambdas))
  
  # 外源项积分
  mu_integral <- mu_int_fn(T_max, parameters = parameters)
  
  # 核函数积分
  kernel_integral <- sum(int_double_ray_direct(T_max, events, alpha1, alpha2, delta1, delta2, delay1, delay2))
  
  # 总积分
  total_integral <- mu_integral + kernel_integral
  
  # 负对数似然
  nll <- -(sum_log_lambdas - total_integral)
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}

# 参数设置
delay1 <- 15 # 固定急性期延迟
delay2 <- 45 # 固定复发期延迟


# 读取单峰拟合结果作为初始参考
single_fit <- readRDS("output/pv_fit_2016.rds")
cat(sprintf("Unimodal Example: α=%.4f, δ=%.4f\n", single_fit$alpha, single_fit$delta))

# 多组初始值尝试拟合
initial_params_list <- list(
  c(0.01, 0.008, 0.05, 0.01, 0.1, 0, 0.1, 0), 
  c(0.008, 0.006, 0.03, 0.008, 0.08, 0.01, 0.08, 0.01), 
  c(0.015, 0.01, 0.08, 0.015, 0.12, -0.01, 0.12, -0.01), 
  c(single_fit$alpha * 0.6, single_fit$alpha * 0.4, single_fit$delta * 0.8, single_fit$delta * 0.6, single_fit$A, single_fit$B, single_fit$M, single_fit$N) # 基于单峰
)

best_fit <- NULL
best_nll <- Inf

for (i in 1:length(initial_params_list)) {
  
  current_params <- initial_params_list[[i]]
  
  tryCatch(
    {
      fit_result <- optimx(
        par = current_params,
        fn = neg_log_likelihood_clean,
        method = "BFGS",
        control = list(trace = 0, maxit = 500)
      )
      
      current_nll <- as.numeric(fit_result$value)
      
      if (is.finite(current_nll) && current_nll < best_nll) {
        best_nll <- current_nll
        best_fit <- fit_result
        cat(sprintf("NLL = %.2f\n", current_nll))
      } else {
        cat(sprintf("NLL = %.2f\n", current_nll))
      }
    },
    error = function(e) {
      cat(sprintf("  拟合失败: %s\n", e$message))
    }
  )
}

if (!is.null(best_fit)) {
  # 提取最佳参数
  best_params <- as.numeric(best_fit[1, 1:8])
  names(best_params) <- c("alpha1", "alpha2", "delta1", "delta2", "A", "B", "M", "N")
  
  cat("\n=== 最佳拟合结果 ===\n")
  cat(sprintf("Neg_Log_Negative: %.2f\n", best_nll))
  cat(sprintf("α₁: %.6f\n", best_params["alpha1"]))
  cat(sprintf("α₂: %.6f\n", best_params["alpha2"]))
  cat(sprintf("δ₁: %.6f\n", best_params["delta1"]))
  cat(sprintf("δ₂: %.6f\n", best_params["delta2"]))
  cat(sprintf("A: %.6f\n", best_params["A"]))
  cat(sprintf("B: %.6f\n", best_params["B"]))
  cat(sprintf("M: %.6f\n", best_params["M"]))
  cat(sprintf("N: %.6f\n", best_params["N"]))
  
  # 计算分支因子
  b1_fitted <- best_params["alpha1"] / best_params["delta1"]
  b2_fitted <- best_params["alpha2"] / best_params["delta2"]
  total_b_fitted <- b1_fitted + b2_fitted
  
  cat(sprintf("(n₁*): %.4f\n", b1_fitted))
  cat(sprintf("(n₂*): %.4f\n", b2_fitted))
  cat(sprintf("(n*): %.4f\n", total_b_fitted))
  
  # 与单峰比较
  single_b <- single_fit$alpha / single_fit$delta
  cat(sprintf("\nUnimodal Branching Factor: %.4f\n", single_b))
  cat(sprintf("\nBimodal Branching Factor: %.4f\n", total_b_fitted))
  cat(sprintf("\nImproved: %.1f%%\n", (total_b_fitted - single_b) / single_b * 100))
  
  # 保存拟合结果
  fitted_params <- as.list(best_params)
  fitted_params$delay1 <- delay1
  fitted_params$delay2 <- delay2
  fitted_params$nll <- best_nll
  fitted_params$branching_factor <- total_b_fitted
  
  saveRDS(fitted_params, "output/pv_double_fit_2016.rds")
  
  # 生成核函数可视化
  t_seq <- seq(0, 200, length.out = 1000)
  kernel_total <- double_ray_kernel(
    t_seq, best_params["alpha1"], best_params["alpha2"],
    best_params["delta1"], best_params["delta2"], delay1, delay2
  )
  peak1 <- ifelse(
    t_seq > delay1,
    best_params["alpha1"] * (t_seq - delay1) * exp(-0.5 * best_params["delta1"] * (t_seq - delay1)^2),
    0
  )
  peak2 <- ifelse(
    t_seq > delay2,
    best_params["alpha2"] * (t_seq - delay2) * exp(-0.5 * best_params["delta2"] * (t_seq - delay2)^2),
    0
  )

  df_relapse <- data.frame(time = t_seq, value = peak2)
  df_acute <- data.frame(time = t_seq, value = peak1)
  ymax <- max(c(df_acute$value, df_relapse$value))
  
  # Acute Peak
  AXIS_TITLE <- 16
  AXIS_TICK  <- 12
  
  theme_update(
    axis.title.x = element_text(size = AXIS_TITLE, face = "bold"),
    axis.title.y = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text.x  = element_text(size = AXIS_TICK),
    axis.text.y  = element_text(size = AXIS_TICK)
  )
  
  p_acute <- ggplot(df_acute, aes(x = time, y = value)) +
    geom_line(linewidth = 1.2, color = "red") +
    labs(x = "Time (days)", y = "Kernel Intensity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    scale_x_continuous(breaks = seq(0, 200, 40)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                       limits = c(0, ymax)) 
  ggsave("figures/pv_double_peak_acute.png", p_acute, width = 10, height = 6, dpi = 300)
  
  # Relapse Peak
  AXIS_TITLE <- 16
  AXIS_TICK  <- 12
  
  theme_update(
    axis.title.x = element_text(size = AXIS_TITLE, face = "bold"),
    axis.title.y = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text.x  = element_text(size = AXIS_TICK),
    axis.text.y  = element_text(size = AXIS_TICK)
  )
  
  p_relapse <- ggplot(df_relapse, aes(x = time, y = value)) +
    geom_line(linewidth = 1.2, color = "blue") +
    labs(x = "Time (days)", y = "Kernel Intensity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    scale_x_continuous(breaks = seq(0, 200, 40)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                       limits = c(0, ymax))   
  ggsave("figures/pv_double_peak_relapse.png", p_relapse, width = 10, height = 6, dpi = 300)
  
  # Acute vs Relapse
  AXIS_TITLE <- 16
  AXIS_TICK  <- 12
  
  theme_update(
    axis.title.x = element_text(size = AXIS_TITLE, face = "bold"),
    axis.title.y = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text.x  = element_text(size = AXIS_TICK),
    axis.text.y  = element_text(size = AXIS_TICK)
  )
  
  df_two <- rbind(
    data.frame(time = t_seq, value = peak1, type = "Acute Peak"),
    data.frame(time = t_seq, value = peak2, type = "Relapse Peak")
  )
  df_two$type <- factor(df_two$type, levels = c("Acute Peak", "Relapse Peak"))
  p_two <- ggplot(df_two, aes(x = time, y = value, color = type)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    labs(x = "Time (days)", y = "Kernel Intensity", color = "Component") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    scale_x_continuous(breaks = seq(0, 200, 40)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    scale_color_manual(values = c("Acute Peak" = "red", "Relapse Peak" = "blue"))
  ggsave("figures/pv_double_two_peaks.png", p_two, width = 10, height = 6, dpi = 300)
  
  # Combined
  AXIS_TITLE <- 16
  AXIS_TICK  <- 12
  
  theme_update(
    axis.title.x = element_text(size = AXIS_TITLE, face = "bold"),
    axis.title.y = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text.x  = element_text(size = AXIS_TICK),
    axis.text.y  = element_text(size = AXIS_TICK)
  )
  
  df_combined <- data.frame(time = t_seq, value = kernel_total)
  p_combined <- ggplot(df_combined, aes(x = time, y = value)) +
    geom_line(linewidth = 1.5, linetype = "solid", color = "black") +
    labs(x = "Time (days)", y = "Kernel Intensity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    scale_x_continuous(breaks = seq(0, 200, 40)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
  ggsave("figures/pv_double_combined.png", p_combined, width = 10, height = 6, dpi = 300)
  
} else {
  cat("Failed\n")
}
