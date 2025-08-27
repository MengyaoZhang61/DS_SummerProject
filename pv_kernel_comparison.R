library(epihawkes)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

pv1_parameters <- as.list(readRDS("output/pv_fit_2016.rds"))
pv1_parameters$delay <- 15

# 双峰：alpha1、alpha2、delta1、delta2、A、B、M、N、delay1、delay2
pv2_parameters <- as.list(readRDS("output/pv_double_fit_2016.rds"))
if (is.null(pv2_parameters$delay1)) pv2_parameters$delay1 <- 15
if (is.null(pv2_parameters$delay2)) pv2_parameters$delay2 <- 45

#定义双峰Rayleigh核
double_ray_kernel_fun <- function(t, parameters){
  t1 <- t - parameters$delay1
  t2 <- t - parameters$delay2
  k1 <- ifelse(t1 > 0, parameters$alpha1 * t1 * exp(-0.5 * parameters$delta1 * t1^2), 0)
  k2 <- ifelse(t2 > 0, parameters$alpha2 * t2 * exp(-0.5 * parameters$delta2 * t2^2), 0)
  k1 + k2
}

#核函数对比：单峰 vs 双峰
t_end <- max(80, pv2_parameters$delay2 + 30) 
t <- seq(0, t_end, length.out = 601)

y_pv_single <- vapply(t, function(x) ray_kernel(x, parameters = pv1_parameters), numeric(1))
y_pv_double <- vapply(t, function(x) double_ray_kernel_fun(x, parameters = pv2_parameters), numeric(1))

df_k <- data.frame(
  t = t,
  `Vivax (Unimodal Kernel)` = y_pv_single,
  `Vivax (Bimodal Kernel)` = y_pv_double
)
df_k_long <- gather(df_k, key = "Model", value = "value", -t)

p1 <- ggplot(df_k_long, aes(t, value, group = Model)) +
  geom_line(aes(col = Model, linetype = Model)) +
  theme_bw() +
  ylab("Kernel intensity") +
  xlab("Time since symptom onset (days)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom", plot.margin = unit(c(.2,.5,.2,.2), "cm"))

#importation基线强度μ(t)对比
ts <- 1:1000
mus_pv_single <- mu_sinusoidal_linear(ts, parameters = pv1_parameters)
mus_pv_double <- mu_sinusoidal_linear(ts, parameters = pv2_parameters)

df_mu <- data.frame(
  t = ts,
  `Vivax (Single-peak)` = mus_pv_single,
  `Vivax (Double-peak)` = mus_pv_double
)
df_mu_long <- gather(df_mu, key = "Model", value = "value", -t)

p2 <- ggplot(df_mu_long, aes(t, value, group = Model)) +
  geom_line(aes(col = Model, linetype = Model)) +
  theme_bw() +
  ylab("Importation intensity") +
  xlab("Time (days)") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom", plot.margin = unit(c(.2,.5,.2,.2), "cm"))

library(dplyr)

peaks <- df_mu_long %>%
  group_by(Model) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup()

print(peaks)

cat(sprintf("Vivax (Single-peak) max μ = %.6f at t = %d\n",
            peaks$value[peaks$Model=="Vivax (Single-peak)"],
            peaks$t[peaks$Model=="Vivax (Single-peak)"]))
cat(sprintf("Vivax (Double-peak) max μ = %.6f at t = %d\n",
            peaks$value[peaks$Model=="Vivax (Double-peak)"],
            peaks$t[peaks$Model=="Vivax (Double-peak)"]))

p <- ggarrange(p1, p2, ncol = 2,
               common.legend = TRUE,
               legend = "bottom",
               labels = "AUTO")

ggsave("figures/fig_4_pv_single_vs_double.pdf", p, width = 12, height = 5)