# ---- Packages ----
# install.packages(c("survival", "survminer", "dplyr"))
# install.packages("pracma")  # for trapz (RMST calc)
library(survival)
library(survminer)
library(dplyr)
library(pracma)

# ---- Your data ----
# months to death/last-contact; censor: 1 = death, 0 = censored (alive)
months <- c(63.5, 42.5, 49.5, 50.5, 6.0, 6.7, 42.5, 31.4, 27.7, 14.7)
censor <- c(0,     0,    0,    0,    0,   1,   0,    0,    0,    1)   # 1=death, 0=censored

time_years <- months/12
status     <- as.integer(censor == 1)  # 1=event (death), 0=censored

df <- data.frame(id = paste0(c("01P","02H","03H","04H","05P","07P","08H","10H","11H","13P")),
                 time = time_years,
                 status = status)

# ---- KM fit ----
fit <- survfit(Surv(time, status) ~ 1, data = df)

# ---- Estimated (BODE-weighted) curve to overlay ----
est_t <- c(0, 1,    2,     3,     4)
est_s <- c(1.0, 0.86, 0.762, 0.675, 0.604)
est_df <- data.frame(time = est_t, surv = est_s)

# ---- Plot with ggsurvplot and overlay the estimate ----
p <- ggsurvplot(
  fit,
  conf.int = FALSE,
  risk.table = TRUE,
  break.time.by = 0.5,
  xlim = c(0, 6),
  ylim = c(0, 1.05),
  xlab = "Years since procedure",
  ylab = "Survival probability",
  ggtheme = theme_minimal(base_size = 12),
  surv.scale = "percent"
)

# add the estimated curve as a dashed step line
p$plot <- p$plot +
  geom_step(data = est_df, aes(x = time, y = surv),
            inherit.aes = FALSE, direction = "hv",
            linewidth = 1.05, linetype = "dashed") +
  annotate("text", x = 3.7, y = 0.62, label = "Estimated (BODE-weighted)", hjust = 0) +
  annotate("text", x = 0.2, y = 0.98, label = "Actual COVE KM", hjust = 0)

print(p)

# ---- Numbers to report ----

# Survival at 4 years
S4 <- tryCatch(summary(fit, times = 4)$surv, error = function(e) NA_real_)

# Median survival (will be NA if not reached)
km_median <- as.numeric(summary(fit)$table["median"])

# RMST to 4 years (numerical integration of KM step function)
# Get KM step function as a dense grid
sf <- summary(fit)
# Build a grid from 0 to tau, and piecewise-constant survival between event times
tau <- 4
grid <- seq(0, tau, by = 0.01)
# Evaluate survival at grid via right-continuous step function
# use last survival value prior to (or at) each grid time
step_times <- c(0, sf$time)
step_surv  <- c(1, sf$surv)
# for grid points earlier than first event time, survival is 1
surv_on_grid <- sapply(grid, function(g) step_surv[max(which(step_times <= g))])
RMST_0_4 <- trapz(grid, surv_on_grid)

cat(sprintf("\n--- Summary ---\nS(4y): %s\nMedian: %s\nRMST(0–4y): %.2f years\n",
            ifelse(is.na(S4), "not estimable", sprintf("%.1f%%", 100*S4)),
            ifelse(is.na(km_median), "Not reached", sprintf("%.2f years", km_median)),
            RMST_0_4))

# ---- (Optional) Export figure ----
# ggsave("COVE_actual_vs_estimated.png", plot = p$plot, width = 8, height = 5, dpi = 300)
