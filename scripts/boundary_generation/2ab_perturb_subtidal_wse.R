# Perturbing historical tides

# WARNING: THIS WILL OVERWRITE TIMESERIES RESULTS AND THUS CHANGE THE BOUNDARY INPUTS MOVING FORWARD

# Libraries commonly used -------------------------------------------------

library(lubridate)
library(reshape2)
library(ggplot2)
library(plotly)
library(zoo)
library(hash)
# library(readxl)
library(htmlwidgets)
library(stringr)
# library(scales)

# Changing this requires an install from GitHub
# https://github.com/CADWRDeltaModeling/stochCycle.git
# install.packages('devtools')
# devtools::install_github("CADWRDeltaModeling/stochCycle")
suppressMessages(library(stochCycle, quietly=TRUE))
doingStochCycle=TRUE

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
source("./functions/pulse_events.R")
# source("./functions/filter.R")



# Load data ---------------------------------------------------------------

version <- 'v1'

# load mrz stage
# mrz_stage <- read.table("./input/dsm2_mrz_stage_ft_15min.csv", skip=4, sep=',', header=FALSE)
# 
# mrz_stage$datetime <- lubridate::dmy_hm(paste0(mrz_stage$V2,' ',mrz_stage$V3))
# mrz_stage <- mrz_stage[,c('datetime','V4')]
# names(mrz_stage) <- c('datetime','stage (ft)')
# mrz_stage$datetime <- format(mrz_stage$datetime, "%Y-%m-%d %H:%M:%S")
# write.table(mrz_stage, "./input/dsm2_mrz_stage_ft_15min_clean.csv", sep=',',
#             row.names=FALSE, col.names=TRUE)

mrz_stage <- read.table("./input/dsm2_mrz_stage_ft_15min_clean.csv", sep=',', header=TRUE)
names(mrz_stage) <- c('datetime','stage (ft)')
mrz_stage$datetime <- lubridate::ymd_hms(mrz_stage$datetime)

mrz_filter <- read.table("./data_out/dsm2_mrz_stage_ft_15min_clean_filtered.csv", sep=',', header=TRUE)
mrz_filter$datetime <- lubridate::ymd_hm(mrz_filter$datetime)

mrz_tidal <- merge(mrz_stage, mrz_filter, by='datetime')
names(mrz_tidal) <-  c('datetime','inst','filtered')
mrz_tidal$tidal <- mrz_tidal$inst - mrz_tidal$filtered
mrz_tidal <- mrz_tidal[,c('datetime','tidal')]

# Plot filtered ts --------------------------------------------------------

plt.filt <- merge(mrz_stage, mrz_filter, by='datetime')
names(plt.filt) <-  c('datetime','inst','filtered')
plt.filt <- melt(plt.filt, id.vars='datetime')


plt.filt <- plt.filt[plt.filt$datetime>=lubridate::dmy_hm('15-01-2000 00:00') &
                       plt.filt$datetime<=lubridate::dmy_hm('30-06-2000 00:00'),] 

plt <- ggplot() +
  geom_line(data=plt.filt, mapping=aes(x=datetime, y=value, color=variable)) +
  scale_x_datetime(date_breaks = "1 week", date_labels = "%b-%d-%Y", expand=c(0,0)) +
  labs(x='',y='Stage (ft)') +
  theme(text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        panel.grid.major.x = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        legend.position='bottom',
        legend.box='vertical',
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        axis.title.y = element_text(color='black'),
        axis.text.x = element_text(angle=90)) +
  ggtitle("Altered Stage at mrz")
plt


# Set parameters ----------------------------------------------------------

nstep <- nrow(mrz_filter)
ints <- mrz_filter$`stage..ft.`
nstep_buf <- nstep+2000
rel_mag <- 0.5
center <- mean(ints)
sigma2_kappa <- c(1e-10,1e-11,1e-12,1e-13)

lambda = 2*pi/(c(2,4,7,14) * (24 * 60 * 15))

stoch <- stoch_tide(nstep, ints, rel_mag, center, lambda,
                    nstep_buf=nstep_buf, sigma2_kappa=sigma2_kappa)

plt.df <- mrz_filter
plt.df$stoch <- stoch

# Plots
plt.df <- melt(plt.df, id.vars='datetime')
plt.df <- plt.df[plt.df$datetime>=lubridate::dmy_hm('15-01-2000 00:00') &
                   plt.df$datetime<=lubridate::dmy_hm('30-12-2002 00:00'),] 

plt.df$value <- as.numeric(plt.df$value)

plt <- ggplot() +
  geom_line(data=plt.df, mapping=aes(x=datetime, y=value, color=variable)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%b-%d-%Y", expand=c(0,0)) +
  scale_color_manual(values=c('black','green','red'),
                     breaks=c('stage (ft)','stoch','Edit')) +
  labs(x='',y='Stage (ft)') +
  theme(text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        panel.grid.major.x = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        legend.position='bottom',
        legend.box='vertical',
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        axis.title.y = element_text(color='black'),
        axis.text.x = element_text(angle=90)) +
  ggtitle("Altered Stage at mrz")
# plt

mean(ints)
mean(stoch)

out.df <- mrz_tidal
out.df$pert_subtide <- stoch
out.df$Edit <- out.df$pert_subtide + out.df$tidal
out.df <- out.df[,c('datetime','Edit')]

# Write out results -------------------------------------------------------

csv_dir <- "./data_out/"

csv_df <- out.df[,c('datetime','Edit')]
csv_df$datetime <- format(csv_df$datetime, "%Y-%m-%d %H:%M:%S")

write.table(csv_df, file=paste0(csv_dir,'/perturb_historical_subtide_',version,'.csv'),
            sep=',', row.names=FALSE, col.names=FALSE)
