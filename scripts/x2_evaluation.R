
# Libraries commonly used -------------------------------------------------

library(lubridate)
library(reshape2)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(zoo)


# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory


# Load obs data -----------------------------------------------------------


dayflow.csv <- "//nasbdo/Modeling_Data/dayflow/dayflow_1983_2023.csv"
dayflow.df <- read.csv(dayflow.csv, skip=1)
dayflow.df$datetime <- lubridate::ymd(dayflow.df$datetime)

head(dayflow.df)

# Get X2 and Outflow ------------------------------------------------------

x2.out.df <- dayflow.df[,c('datetime','X2','OUT')]
x2.out.df <- x2.out.df[complete.cases(x2.out.df),]


# Analysis ----------------------------------------------------------------


# Function to compute NDAYS for a given threshold
compute_ndays <- function(df, threshold) {
  df %>%
    mutate(below_thresh = X2 < threshold) %>%
    group_by(group = cumsum(c(1, diff(below_thresh)) * below_thresh)) %>%
    mutate(!!paste0("NDAYS", threshold) := ifelse(below_thresh, n(), 0)) %>%
    ungroup() %>%
    select(-group, -below_thresh)
}

# Apply the function for different thresholds
thresholds <- c(75, 70, 65, 60, 55, 50, 45)
for (thresh in thresholds) {
  x2.out.df <- compute_ndays(x2.out.df, thresh)
}


# 70 km threshold ---------------------------------------------------------


df70 <- x2.out.df[,c("datetime", "X2", "OUT", paste0("NDAYS", 70))]

# Compute 30-day rolling average
df70 <- df70 %>%
  arrange(datetime) %>%  # Ensure data is ordered by date
  mutate(OUT_30 = rollmean(OUT, k = 30, fill = NA, align = "right"))
df70 <- df70[complete.cases(df70),]


plt <- ggplot(df70, aes(x = NDAYS70, y = OUT, color = OUT_30)) +
  geom_point(size = 3, alpha = 0.8) +  # Scatter plot with transparency
  labs(
    x = "Number of days under 70km",
    y = "Outflow (cfs)",
    color = "30 Antecedent Average (cfs)"
  ) +
  scale_x_continuous(breaks=seq(0,1500,100), limits=c(0,1500)) +
  scale_y_continuous(breaks=seq(0,575000,75000), limits=c(0,575000)) +
  scale_color_viridis_c(
    breaks = c(3000, 100000, 200000, 270000),  # Custom legend breaks
    option = "C"  # Viridis color palette option
  ) +
  theme(text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        legend.key=element_blank(),
        panel.grid.major.x = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(color='black'),
        legend.position.inside=c(0.975,0.975),
        legend.justification=c(0.975,0.975),
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
  )  
plt

ggplotly(plt)

gplt <- ggplotly(plt, dynamicTicks=TRUE)

pltl_name <- "x2_outflow.html"
saveWidget(gplt, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("./plots/",pltl_name))

