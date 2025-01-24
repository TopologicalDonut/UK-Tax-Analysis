# ---- Setup ----
library(tidyverse)
library(lubridate)
library(ggplot2)
library(tsibble)
library(fable)
library(feasts)
library(tseries)
library(kableExtra)
library(patchwork)
library(here)

TAX_CHANGE_DATE <- as.Date('2021-01-01')
CACHE_DIR <- here('intermediate_outputs')
FIGURE_DIR <- here('figures')

dir.create(CACHE_DIR)
dir.create(FIGURE_DIR)

# ---- Table for Rebase Example ----

example_data_path <- here('data', 'example_ONS_rebasing.csv')
example_data <- read_csv(example_data_path, show_col_types = FALSE)

kbl(example_data,
    col.names = (c("Date", "Item ID", "Item",
                   "Price Index")),
    align = c('l', 'c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 2,
    caption = "Example of ONS Rebasing") %>%
  kable_styling(latex_options = c("hold_position"))

# ---- Analysis ----

analysis_data_path <- here('data', 'analysis_data.csv')
analysis_data <- read_csv(analysis_data_path)

create_item_data <- function(id) {
  analysis_data %>%
    filter(item_id == id) %>%
    mutate(month = tsibble::yearmonth(date),
           tax_dummy = as.integer(date >= as.Date("2021-01-01"))) %>%
    select(-date) %>%
    as_tsibble(index = month) %>%
    tsibble::fill_gaps()
}

analyze_item <- function(item_id) {
  item_data <- create_item_data(item_id)
  item_name <- item_data$item_desc[1]
  
  arima_model <- item_data %>%
    model(ARIMA(ln_rebased_index ~ tax_dummy, stepwise = FALSE))
  
  resid_plot <- gg_tsresiduals(arima_model)[[1]] + 
    ggtitle(paste("Residuals for", item_name))
  
  acf_plot <- gg_tsresiduals(arima_model)[[2]] + 
    ggtitle(paste("ACF Plot for", item_name))

  estimates_info <- tidy(arima_model) %>%
    filter(term == "tax_dummy") %>%
    select(estimate, std.error, p.value)
  
  # Return results
  list(
    item_id = item_id,
    item_name = item_name,
    data = item_data,
    model = arima_model,
    resid_plot = resid_plot,
    acf_plot = acf_plot,
    estimates_info = estimates_info
  )
}

tampon_analysis <- analyze_item(520206)

# ---- ACF Plot ----

tampon_data <- create_item_data(520206)

tampon_data %>% 
  ACF(ln_rebased_index, lag_max = 48) %>%
  autoplot()

# ---- Tampon Analysis Table----

kbl(tampon_analysis$estimates_info,
    col.names = c("Estimate", "Std. Error", "p-value"),
    align = c('c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 4,
    caption = "Estimate of Effect of Tampon Tax Abolition") %>%
  kable_styling(latex_options = c("hold_position"))

# ---- Tampon Resid Graphs ----

tampon_resid_and_ACF <- tampon_analysis$resid_plot / tampon_analysis$acf_plot +
  plot_layout(heights = c(1,1))
print(tampon_resid_and_ACF)

# ---- Tampon Price Graph ---- 

autoplot(tampon_analysis$data, ln_rebased_index) +
  geom_vline(xintercept = as.Date("2021-01-01"), color = "red", linetype = "dashed") +
  labs(title = "Plot of ln Price Index for Tampons",
       y = "ln of Rebased Index")

# ---- Robustness Check ----

item_ids <- c(520213, 430536, 520249, 520241)

output_filename <- "intermediate_outputs/other_items_analysis.rds"

if (!file.exists(output_filename)){
  other_items_analysis <- lapply(item_ids, analyze_item)
  saveRDS(other_items_analysis, file = output_filename)
} else {
  other_items_analysis <- readRDS(output_filename)
}

# ---- Robustness Check Table ----

summary_table <- do.call(rbind, lapply(other_items_analysis, function(x) {
  tibble(item_name = x$item_name,
         estimates_info = x$estimates_info)
  })) %>%
  unnest(cols = estimates_info)

kbl(summary_table,
    col.names = c("Item", "Estimate", "Std. Error", "p-value"),
    align = c('l', 'c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 4,
    caption = "Estimates of Placebo Tax Abolition on Other Items") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  add_header_above(c(" " = 1, "Statistics" = 3))

# ---- Resid Plots for Others ----

residual_plots <- lapply(other_items_analysis, function(x) {
  x$resid_plot})

acf_plots <- lapply(other_items_analysis, function(x) {
  x$acf_plot})

combined_residual_plots <- wrap_plots(residual_plots, nrow = 4)
combined_acf_plots <- wrap_plots(acf_plots, nrow = 4)

ggsave("Figures/residual_plots.png", plot = combined_residual_plots, width = 10, height = 12)
ggsave("Figures/acf_plots.png", plot = combined_acf_plots, width = 10, height = 12)