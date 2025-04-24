# Libraries commonly used -------------------------------------------------

library(lubridate)
library(reshape2)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(zoo)
library(purrr)


# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")  # Adjust the path to your pandoc installation


# Load obs data -----------------------------------------------------------


dayflow.csv <- "//nasbdo/Modeling_Data/dayflow/dayflow_1983_2023.csv"
dayflow.df <- read.csv(dayflow.csv, skip=1)
dayflow.df$datetime <- lubridate::ymd(dayflow.df$datetime)

head(dayflow.df)

# Get X2 and Outflow ------------------------------------------------------

x2.out.df <- dayflow.df[,c('datetime','X2','OUT')]
x2.out.df <- x2.out.df[complete.cases(x2.out.df),]
x2.out.df$OUT <- x2.out.df$OUT/1000 # change outflow units to thousand-cfs


# Analysis ----------------------------------------------------------------

# Create monthly.out.df with monthly average of OUT
monthly.out.df <- x2.out.df %>%
  mutate(year_month = format(datetime, "%Y-%m")) %>%  # Extract year and month
  group_by(year_month) %>%  # Group by year and month
  summarize(out_monthly = mean(OUT, na.rm = TRUE)) %>%  # Calculate monthly average of OUT
  ungroup() %>% 
  mutate(datetime = as.Date(paste0(year_month, "-01"))) %>%  # Convert year_month to datetime
  select(datetime, out_monthly)  # Select relevant columns

# Function to compute NDAYS for a given threshold, standardized to a 30-day month
compute_ndays_monthly <- function(df, threshold) {
  df %>%
    mutate(
      below_thresh = X2 < threshold,
      year_month = format(datetime, "%Y-%m")  # Extract year and month
    ) %>%
    group_by(year_month) %>%
    summarize(
      NDAYS = sum(below_thresh, na.rm = TRUE),  # Count days below the threshold
      actual_days = n()  # Count the actual number of days in the month
    ) %>%
    ungroup() %>%
    mutate(
      standardized_NDAYS = NDAYS * (30 / actual_days)  # Standardize to a 30-day month
    ) %>%
    select(year_month, standardized_NDAYS) %>%  # Keep only relevant columns
    rename(!!paste0(threshold, "_km") := standardized_NDAYS)  # Rename column for the threshold
}

# Apply the function for different thresholds and combine results into a single dataframe
thresholds <- c(80, 74)
monthly_results <- thresholds %>%
  lapply(function(thresh) compute_ndays_monthly(x2.out.df, thresh)) %>%
  reduce(full_join, by = "year_month")  # Combine all results by year_month

# Convert year_month to datetime for better handling
monthly_results <- monthly_results %>%
  mutate(datetime = as.Date(paste0(year_month, "-01"))) %>%
  select(datetime, everything(), -year_month)  # Reorder columns

monthly_results  <- merge(monthly_results, monthly.out.df, by='datetime')

# Add out-1month and out-2month columns
monthly_results <- monthly_results %>%
  arrange(datetime) %>%  # Ensure data is sorted by datetime
  mutate(
    `out-1month` = lag(out_monthly, 1),  # Previous month's OUT_monthly_avg
    `out-2month` = lag(out_monthly, 2)   # Two months ago's OUT_monthly_avg
  )

# Reshape the data to long format for plotting
plot.df <- melt(
  monthly_results,
  id.vars = c("datetime", "out_monthly", "out-1month", "out-2month"),  # Include out-1month and out-2month as identifiers
  measure.vars = c("80_km", "74_km"),  # *_km columns
  variable.name = "threshold",  # Name for the threshold column
  value.name = "days_under_threshold"  # Name for the number of days column
)

# # Log-transform the out_monthly column for ggplotly compatibility
# plot.df <- plot.df %>%
#   mutate(log_out_monthly = log10(out_monthly))  # Add a log-transformed column

# Plot Results ---------------------------------------------------------

# Create the plot
plt <- ggplot(plot.df, aes(x = out_monthly, y = days_under_threshold, color = `out-1month`)) +
  geom_point(size = 3, alpha = 0.8) +  # Scatter plot with transparency
  labs(
    x = "Monthly Outflow (log10 scale, 1000-cfs)",
    y = "Number of Days Under Threshold",
    color = "Previous Month (1000-cfs)"
  ) +  
  scale_color_gradientn(
    colors = c("beige", "khaki", "darkgreen", "darkblue"),  # Define a gradient from beige to dark blue
    trans = "log",  # Apply a logarithmic transformation to the color scale
    limits = c(1, 75),  # Set the range of the color scale (log scale cannot include 0)
    oob = scales::squish,  # Squish values outside the range into the limits
    breaks = c(1,5,10,25,50,70)
  ) +
  scale_x_log10() +
  scale_y_continuous(
    limits = c(0, 31),  # Set y-axis limits
    breaks = seq(0, 30, 5)  # Define y-axis breaks
  ) +
  facet_wrap(~threshold, ncol = 1, scales = "fixed") +  # Ensure matching scales across facets
  theme(
    text = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = 'grey80', linetype = 'dashed'),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = 'grey80', linetype = 'dashed'),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(color = 'black'),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NULL)
  )
# plt

gplt <- ggplotly(plt, dynamicTicks=TRUE, tooltip=c('days_under_threshold','out_monthly',"out-1month", "out-2month"))
gplt

# Save the interactive plot as an HTML file -------------------------
pltl_name <- "x2_outflow.html"
saveWidget(gplt, pltl_name, selfcontained = TRUE)
file.rename(pltl_name, paste0("./plots/", pltl_name))

