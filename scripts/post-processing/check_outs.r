
# Libraries commonly used -------------------------------------------------

library(lubridate)
library(reshape2)
library(ggplot2)
library(plotly)
# library(zoo)
# library(hash)
# library(readxl)
library(htmlwidgets)
# library(stringr)
# library(scales)


# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory

# Load data ---------------------------------------------------------------

casanntra.dir <- '../../../casanntra/data'
# data.dir <- '../ann_training/data_out'

cases <- seq(1,100)

for (case in cases) {
  df <- read.csv(paste0(casanntra.dir,'/dsm2_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df <- df[,c('datetime','x2','ndo','case')]
  # df <- df[-(1:30),c('datetime','x2','ndo','case')] # remove first 30 days (spinup) and subset variables
  
  if (case==1) {
    x2_ndo_df <- df
  } else {
    x2_ndo_df <- rbind(x2_ndo_df, df)
  }
}
x2_ndo_df$case <- as.factor(x2_ndo_df$case)

plt <- ggplot(x2_ndo_df, aes(x=ndo, y=x2, color=case)) +
  geom_point()

plt

ggplotly(plt)


# check specific cases ----------------------------------------------------

# cases in the 90s look like there's an anomoly for high ndo and high x2 (expected inverse relationship)

plt_cases <- x2_ndo_df[x2_ndo_df$case %in% seq(91,100),]
# plt_cases <- x2_ndo_df

plt <- ggplot(plt_cases, aes(x=datetime, y=x2, color=case)) +
  geom_line()

plt

ggplotly(plt)

# the anomoly was shown in the first 30 days. Looking at the model spinup, this is a remnant of high salinity


# Plot a few EC out locs --------------------------------------------------


cases <- c(1,11,21,31,41,51,61,71,81,91)
# cases <- c(1,11,21,31,41,51)
# cases <- seq(1,10)
cases <- c(23,45)
cases <- c(23,45,seq(101,107))
cases <- c(23,seq(101,104))
cases <- c(45,seq(105,107))

for (case in cases) {
  df <- read.csv(paste0(casanntra.dir,'/dsm2_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df <- df[,c('datetime','case',"mrz_tidal_energy", "mrz_tidal_filter", 
              'sf_tidal_energy','sf_tidal_filter')] # subset variables
  # df <- df[,c('datetime','case','ndo', 'cu_flow', 'sjr_flow','sf_tidal_energy','sf_tidal_filter',"mrz_tidal_energy", "mrz_tidal_filter",'emm2')] # subset variables
  
  if (case==cases[1]) {
    var_df <- df
  } else {
    var_df <- rbind(var_df, df)
  }
}
var_df$case <- as.factor(var_df$case)
var_df <- melt(var_df, id.vars=c('datetime','case'))

plt <- ggplot(var_df, aes(x=datetime, y=value, color=variable)) +
  facet_wrap(~case, ncol=1, scales='free_y') +
  geom_line()

# plt

ggplotly(plt, dynamicTicks=TRUE)

# plt <- ggplot(var_df[var_df$variable %in% c('sf_tidal_energy','sf_tidal_filter',"mrz_tidal_energy", "mrz_tidal_filter"),], aes(x=datetime, y=value, color=variable)) +
#   geom_line()
# 
# ggplotly(plt)
