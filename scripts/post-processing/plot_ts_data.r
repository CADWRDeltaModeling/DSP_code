# Checking


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

out.dir <- './data_out/'
experiment <- 'schism_lhc_v3'

# SUISUN ------------------------------------------------------------

suisun_basebc <- read.csv(paste0(out.dir,experiment,'/suisun_basebc_lhc_4.csv'))
suisun_basebc$Time <- lubridate::ymd_hms(suisun_basebc$X0)

suisun_suisbc <- read.csv(paste0(out.dir,experiment,'/suisun_lhc_4.csv'))
suisun_suisbc$Time <- lubridate::ymd_hms(suisun_suisbc$X0)

names(suisun_basebc)
extraction_vars <- c('Time','Net.Delta.Outflow.Out','Rio.Vista.Flow.Out',
                     'emm2_upper.EC','mal_upper.EC','carqb_upper.EC',
                     'gzm_default.EC','msl_default.EC','oh4_default.EC')

melt_sbb <- melt(suisun_basebc[,extraction_vars], id.vars='Time')
melt_sbb$Barotropic <- 'Baseline'
melt_ssb <- melt(suisun_suisbc[,extraction_vars], id.vars='Time')
melt_ssb$Barotropic <- 'Suisun'

plt.df <- rbind(melt_sbb, melt_ssb)

plt <- ggplot(data=plt.df, aes(x=Time, y=value, color=Barotropic)) +
  geom_line() +
  facet_wrap(~variable, ncol=1, scales='free_y')


pltl <- ggplotly(plt, dynamicTicks=TRUE)

pltl_name <- paste0("Suisun_Barotropic_Check_",experiment,".html")

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))


# Calculate differences ---------------------------------------------------

diff.df <- data.frame(Time=suisun_basebc$Time)

for (var in extraction_vars[2:length(extraction_vars)]) {
  print(var)
  varname <- paste0(var,' (Baseline-Suisun)')
  interm.df <- merge(suisun_basebc[,c('Time',var)],suisun_suisbc[,c('Time',var)], by='Time')
  interm.df[,varname] <- interm.df[2] - interm.df[3]
  diff.df <- merge(diff.df, interm.df[,c('Time',varname)], by='Time')
  print(paste0("Max= ",round(max(diff.df[,varname]),2),
               " at ",diff.df$Time[which.max(diff.df[,varname])], 
               ", Min= ",round(min(diff.df[,varname]),2),
               " at ",diff.df$Time[which.min(diff.df[,varname])]))
}
