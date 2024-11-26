
# Set up --------------------------------------------------

# load libraries
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(reshape2)
library(lubridate)
library(stringr)
library(scales)

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory


# Define variables ---------------------------------------------------------

lhc_5_sim_ccfb_fn <- 'C:/Users/tomkovic/AppData/Local/Temp/scp56588/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/ccfb_gate.th'
lhc_5_ccfb_fn <- 'C:/Users/tomkovic/AppData/Local/Temp/scp23588/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/lhc_5/ccfb_lhc_5.th'

start_date <- lubridate::ymd_hm('2006-11-14 00:00')
time_fail <- as.POSIXct(start_date) + 43295670. 

cc.head <- c('datetime','install','ndup','op_down','op_up','elev','width','height')
# Read in data ------------------------------------------------------------

cc_sim.df <- read.table(lhc_5_sim_ccfb_fn)
cc.df <- read.table(lhc_5_ccfb_fn)

names(cc_sim.df) <- cc.head
names(cc.df) <- cc.head

cc_sim.df$datetime <- as.POSIXct(start_date) + as.numeric(cc_sim.df$datetime)
cc.df$datetime <- as.POSIXct(start_date) + as.numeric(cc.df$datetime)


# Refine & combine data ---------------------------------------------------

cc_sim.df$eff_op <- cc_sim.df$ndup * cc_sim.df$height
cc.df$eff_op <- cc.df$ndup * cc.df$height

cc_sim.df <- melt(cc_sim.df[,c('datetime','eff_op')], id.vars='datetime')
cc_sim.df$case <- 'in sim folder'
cc.df <- melt(cc.df[,c('datetime','eff_op')], id.vars='datetime')
cc.df$case <- 'in case folder'

plt.cc <- rbind(cc_sim.df,cc.df)
plt.cc <- plt.cc[plt.cc$datetime<lubridate::ymd_hm('2009-01-01 00:00'),]

plt <- ggplot(data=plt.cc) + 
  geom_line(aes(x=datetime, y=value, color=case)) +
  geom_vline(xintercept = time_fail, color = "red", linetype = "dashed") +
  scale_x_datetime(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle=90),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        legend.key=element_blank(),
        panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        legend.position.inside=c(0.975,0.975),
        legend.justification=c(0.975,0.975),
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        legend.title=element_blank(),
        axis.title.y = element_text(color='black')
  )

ggplotly(plt, dynamicTicks=TRUE)
