library(tidyverse)
library(lubridate)
library(ggplot2)
library(tsibble)
library(fable)
library(feasts)
library(tseries)
library(patchwork)
library(kableExtra)

# ---- Merge and Clean Raw Data ----

data_list <- list.files("Data/Raw", pattern = "*.csv", full.names = TRUE) %>%
  setdiff("Data/Raw/CPI.csv")

product_data <- data_list %>% 
  map_df(~read_csv(.x) %>% 
           rename_all(toupper))

product_data_clean <- product_data %>%
  select(INDEX_DATE, ITEM_ID, ITEM_DESC, ALL_GM_INDEX) %>%
  mutate(INDEX_DATE = ym(INDEX_DATE))

write_csv(product_data_clean, "Data/Processed/merged_product_data_clean.csv")

# ---- Table for Rebase Example ----

product_data <- read_csv("Data/Processed/merged_product_data_clean.csv", show_col_types = FALSE)

example_data <- product_data %>%
  filter(ITEM_ID == 610310) %>%
  filter(month(INDEX_DATE) %in% c(12, 1, 2, 3) & year(INDEX_DATE) < 2010) %>%
  slice(3:10)

kbl(example_data,
    col.names = (c("Index Date", "Item ID", "Item",
                   "Price Index")),
    align = c('l', 'c', 'c', 'c'),
    booktabs = T,
    linesep = "",
    digits = 2,
    caption = "Example of ONS Rebasing") %>%
  kable_styling(latex_options = c("hold_position"))

# ---- Analysis ----

create_item_data <- function(item_id) {
  product_data %>%
    filter(ITEM_ID == item_id) %>%
    mutate(Month = tsibble::yearmonth(INDEX_DATE),
           tax_dummy = as.integer(INDEX_DATE >= as.Date("2021-01-01"))) %>%
    select(-INDEX_DATE) %>%
    as_tsibble(index = Month)
}

rebase_cpi <- function(data){
  january_indices <- data %>%
    filter(month(Month) == 1) %>%
    mutate(year = year(Month),
           jan_index = ALL_GM_INDEX) %>%
    as_tibble() %>%
    select(year, jan_index) %>%
    mutate(
      cumulative_factor = cumprod(jan_index / 100)
    )
  
  data_rebased_cpi <- data %>%
    mutate(year = year(Month)) %>%
    left_join(january_indices, by = "year") %>%
    mutate(
      cumulative_factor = if_else(month(Month) == 1, lag(cumulative_factor), cumulative_factor) %>%
        replace_na(1),
      rebased_index = ALL_GM_INDEX * cumulative_factor,
      log_rebased_index = log(rebased_index)
    ) %>%
    select(-year, -ALL_GM_INDEX, -jan_index, -cumulative_factor)
  return(data_rebased_cpi)
}

analyze_item <- function(item_id) {
  item_data <- create_item_data(item_id)
  item_name <- item_data$ITEM_DESC[1]
  item_data_rebased <- rebase_cpi(item_data)
  
  arima_model <- item_data_rebased %>%
    model(ARIMA(log_rebased_index ~ tax_dummy, stepwise = FALSE))
  
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
    data = item_data_rebased,
    model = arima_model,
    resid_plot = resid_plot,
    acf_plot = acf_plot,
    estimates_info = estimates_info
  )
}

tampon_analysis <- analyze_item(520206)

# ---- ACF Plot ----

tampon_data <- create_item_data(520206) %>%
  rebase_cpi()

tampon_data %>% 
  ACF(log_rebased_index, lag_max = 48) %>%
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

autoplot(tampon_analysis$data, log_rebased_index) +
  geom_vline(xintercept = as.Date("2021-01-01"), color = "red", linetype = "dashed") +
  labs(title = "Plot of ln Price Index for Tampons",
       y = "ln of Rebased Index")

# ---- Robustness Check ----

item_ids <- c(520213, 430536, 520249, 520241)

output_filename <- "Intermediate_outputs/other_items_analysis.rds"

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