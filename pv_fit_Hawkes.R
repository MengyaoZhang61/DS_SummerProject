library(epihawkes)
library(ggplot2)
library(optimx)

data = readRDS("data/pv_onset_2016.rds")

# Sets the seed
set.seed(1)

# Sets parameters for the exogenous term
mu_term <- "sinusoidal_linear"
mu_fn <- mu_sinusoidal_linear
mu_diff_fn <- mu_diff_sinusoidal_linear
mu_int_fn <- mu_int_sinusoidal_linear

# Prints log level
print_level <- 1

n_fits = 10
outputs <- vector("list", length = n_fits)
nlls <- vector(length = n_fits)
alpha <- vector(length = n_fits)
delta <- vector(length = n_fits)

for(i in 1:n_fits){
  parameters <- list(alpha = runif(1, 0, 1),
                     delta = runif(1, 0, 1),
                     A = runif(1, 0, 1),
                     B = runif(1, 0, 1),
                     M = runif(1, 0, 1), 
                     N = runif(1, 0, 1))
  
  outputs[[i]] <- optimx(par = unlist(parameters), fn = neg_log_likelihood, gr = ray_derivatives,
                         method = "BFGS",
                         events = as.numeric(data$time), 
                         delay = 15,
                         kernel = ray_kernel, 
                         mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                         mu_int_fn = mu_int_fn,
                         print_level = print_level)
  
  nlls[i] = outputs[[i]]$value
  alpha[i] = outputs[[i]]$alpha
  delta[i] = outputs[[i]]$delta
  print(sprintf("neg LL: %f", nlls[i]))
  print(paste(c("Optimal parameters:", outputs[[i]][1:6])))
  
}

idx = which(nlls == min(nlls[which(delta >0 & alpha >0)]))
parameters <- outputs[[idx]][1:6] 

saveRDS(parameters, "output/pv_fit_2016.rds")

print("Vivax parameters:")
print(outputs[[idx]][1:6])

parameters <- as.list(readRDS("output/pv_fit_2016.rds"))
parameters$delay = 15
pv_decay_kernel <- plot_decay_kernel(ray_kernel, parameters = parameters, T_max = 50)
print(pv_decay_kernel)
ggsave("figures/pv_decay_kernel.png", pv_decay_kernel, width = 10, height = 6, dpi = 300)

ts = 1:max(as.numeric(data$time))
ys = mu_fn(ts, parameters = parameters)
df <- data.frame(ts, ys)
p <- ggplot(df) + geom_line(aes(ts, ys)) + 
  ylab(expression(mu)) + xlab('t') + theme_bw()
print(p)
ggsave("figures/pv_events.png", p, width = 10, height = 6, dpi = 300)

print("Branching factors")
print(parameters$alpha/parameters$delta)

best_row <- outputs[[idx]]


final_nll <- as.numeric(best_row$value)
final_ll  <- -final_nll

cat(sprintf("Final Negative Log-Likelihood (NLL): %.6f\n", final_nll))
cat(sprintf("Final Log-Likelihood (LL): %.6f\n", final_ll))
