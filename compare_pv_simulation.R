library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(readr)


num_procs   <- 10         
alpha_total <- 0.15      
alpha_imp   <- 0.18  
out_dir     <- "figures"
dir.create(out_dir, showWarnings = FALSE)

export_raw  <- TRUE

# 输出文件
f_side_pdf  <- "figures/fig_side_by_side_single_vs_double.pdf"
f_side_png  <- "figures/fig_side_by_side_single_vs_double.png"
f_band_pdf  <- "figures/fig_import_band_single_vs_double.pdf"
f_band_png  <- "figures/fig_import_band_single_vs_double.png"

real_events      <- readRDS("data/pv_onset_2016.rds")
real_event_times <- real_events$time
T_max            <- max(real_event_times)

real_counts <- epihawkes::compute_count_function(
  events = real_event_times, N = 5000, T_max = T_max
)

real_import_times <- readRDS("output/pv_imported_2016.rds")
real_counts_imported <- epihawkes::compute_count_function(
  events = real_import_times, N = 5000, T_max = T_max
)

interp_to <- function(df_real, times){
  approx(df_real$t, df_real$counts, xout = times, rule = 2)$y
}

# 读取并堆叠“原始模拟矩阵”
load_sims <- function(prefix_fmt, num_procs){
  counts_list     <- vector("list", num_procs)
  counts_imp_list <- vector("list", num_procs)
  time_master     <- NULL
  for(seed in 1:num_procs){
    f <- sprintf(prefix_fmt, seed)
    if(!file.exists(f)) stop("Missing file: ", f)
    load(f) # expects: counts, counts_imp, time_vec
    #[N_time x N_runs]
    counts_list[[seed]]     <- t(counts)
    counts_imp_list[[seed]] <- t(counts_imp)
    if(is.null(time_master)) time_master <- time_vec
  }
  total     <- do.call(cbind, counts_list)      # [N_time x N_traj_total]
  imported  <- do.call(cbind, counts_imp_list)  # [N_time x N_traj_import]
  list(times = time_master, total = total, imported = imported)
}

mse <- function(a,b) mean((a-b)^2)
metrics_block <- function(sim){
  real_total_vec   <- interp_to(real_counts, sim$times)
  real_import_vec  <- interp_to(real_counts_imported, sim$times)
  mse_total   <- apply(sim$total,    2, mse, b = real_total_vec)
  mse_import  <- apply(sim$imported, 2, mse, b = real_import_vec)
  end_real_total  <- tail(real_total_vec, 1)
  end_real_import <- tail(real_import_vec, 1)
  end_total_vec   <- sim$total   [nrow(sim$total), ]
  end_import_vec  <- sim$imported[nrow(sim$imported), ]
  ratio_sim   <- end_import_vec / end_total_vec
  ratio_real  <- end_real_import / end_real_total
  list(
    mse_total  = mse_total,
    mse_import = mse_import,
    end_abs_err_total = abs(end_total_vec - end_real_total),
    end_abs_err_ratio = abs(ratio_sim - ratio_real)
  )
}

to_long_all <- function(mat, times, type_label){
  df <- as.data.frame(mat); df$times <- times
  long <- tidyr::pivot_longer(df, -times, names_to="id", values_to="value")
  long$type <- type_label
  long
}

summ_from_mat <- function(mat, times){
  data.frame(
    times = times,
    mean  = rowMeans(mat),
    lower = apply(mat, 1, quantile, 0.025),
    upper = apply(mat, 1, quantile, 0.975)
  )
}

single <- load_sims("output/pv_sims_%d_2016.Rdata",        num_procs)
double <- load_sims("output/pv_double_sims_%d_2016.Rdata", num_procs)

if (export_raw) {
  write_rds(single, file.path(out_dir, "RAW_single_counts_and_imports.rds"))
  write_rds(double, file.path(out_dir, "RAW_double_counts_and_imports.rds"))
}

mb_single <- metrics_block(single)
mb_double <- metrics_block(double)
impr_mse_total  <- (mean(mb_single$mse_total)  - mean(mb_double$mse_total))  / mean(mb_single$mse_total)  * 100
impr_mse_import <- (mean(mb_single$mse_import) - mean(mb_double$mse_import)) / mean(mb_single$mse_import) * 100
impr_end_total  <- (mean(mb_single$end_abs_err_total) - mean(mb_double$end_abs_err_total)) / mean(mb_single$end_abs_err_total) * 100
impr_end_ratio  <- (mean(mb_single$end_abs_err_ratio) - mean(mb_double$end_abs_err_ratio)) / mean(mb_single$end_abs_err_ratio) * 100

long_single_total <- to_long_all(single$total,    single$times, "total")
long_single_imp   <- to_long_all(single$imported, single$times, "import")
long_double_total <- to_long_all(double$total,    double$times, "total")
long_double_imp   <- to_long_all(double$imported, double$times, "import")

ymax_total <- max(
  real_counts$counts,
  single$total,
  double$total
)

make_panel_raw <- function(long_total, long_imp, title_txt, ylim_top){
  ggplot() +
    geom_line(data = long_total, aes(times, value, group = id),
              colour = "black", alpha = alpha_total, linewidth = 0.25) +
    geom_line(data = long_imp,   aes(times, value, group = id),
              colour = "blue", alpha = alpha_imp, linewidth = 0.25) +
    geom_line(data = real_counts, aes(x = t, y = counts), colour = "red",   linewidth = 0.6) +
    geom_line(data = real_counts_imported, aes(x = t, y = counts), colour = "green4", linewidth = 0.6) +
    labs(x = "Time (days)", y = "Case count") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0, face = "bold")) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(ylim = c(0, ylim_top)) +
    ggtitle(title_txt)
}

p_single <- make_panel_raw(long_single_total, long_single_imp, "A  Unimodal–Rayleigh", ymax_total)
p_double <- make_panel_raw(long_double_total, long_double_imp, "B  Bimodal–Rayleigh", ymax_total)


grid_side  <- plot_grid(p_single, p_double, nrow = 1)
final_side <- ggdraw(grid_side)
ggsave(f_side_pdf, final_side, width = 12, height = 5.8)
ggsave(f_side_png, final_side, width = 12, height = 5.8, dpi = 300)


# 导入“均值+95%CI”
ymax_import <- max(
  real_counts_imported$counts,
  single$imported,
  double$imported
)

make_import_band_panel <- function(sim, tag_letter, ylim_top){
  real_imp_vec <- interp_to(real_counts_imported, sim$times)
  imp_summ <- summ_from_mat(sim$imported, sim$times)
  ggplot() +
    geom_ribbon(data = imp_summ, aes(x = times, ymin = lower, ymax = upper), fill = "royalblue", alpha = 0.25) +
    geom_line(data = imp_summ, aes(x = times, y = mean), colour = "royalblue4", size = 0.7) +
    geom_line(data = data.frame(times = sim$times, real = real_imp_vec),
              aes(x = times, y = real), colour = "green4", size = 0.7) +
    labs(x = "Time (days)", y = "Imported cases (cumulative)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0, face = "bold")) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(ylim = c(0, ylim_top)) +
    ggtitle(tag_letter)
}
p_band_single <- make_import_band_panel(single, "A Import mean ±95%CI (Unimodal)", ymax_import)
p_band_double <- make_import_band_panel(double, "B Import mean ±95%CI (Bimodal)", ymax_import)

grid_band <- plot_grid(p_band_single, p_band_double, nrow = 1)
ggsave(f_band_pdf, grid_band, width = 12, height = 5.2)
ggsave(f_band_png, grid_band, width = 12, height = 5.2, dpi = 300)

fmt_pct <- function(x) sprintf("%.1f%%", x)

improvement_summary <- tibble::tibble(
  metric = c("MSE (total cases)",
             "MSE (importations)",
             "End-point abs. error (total cases)",
             "End-point error (import ratio)"),
  improvement_bimodal_vs_unimodal = c(
    fmt_pct(impr_mse_total),
    fmt_pct(impr_mse_import),
    fmt_pct(impr_end_total),
    fmt_pct(impr_end_ratio)
  )
)

print(improvement_summary)

