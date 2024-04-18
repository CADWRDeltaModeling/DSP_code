# Perturbing historical dataset to create variability in TS

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
source("./functions/load_delta_vars.R")
source("./functions/pulse_events.R")


# Specify overwrite -------------------------------------------------------

overwrite <- c('SJR Flow','Exports','Sacramento')

# Load data ---------------------------------------------------------------

# load delta state vars
dsm2_filename <- '../../data/dsm2_ann_inputs_base.xlsx'

sheetlist <- c("northern_flow","sjr_flow","exports","dxc_gate_fraction",
               "net_delta_cu","mtz_daily_max-min_stage",
               "sjr_vernalis_ec","sac_ec","base_ec_output")
sheetnames <- c('Northern Flow', 'SJR Flow','Exports','DCC Gate',
                'Consump Use','Tidal Amp', 
                'San Joaquin EC','Sacramento EC')

for (sheet in sheetlist){
  assign(sheet,load_delta_vars(dsm2_filename, sheet))
}

# Combine all but base_ec_output
for (var in sheetlist[-length(sheetlist)]) {
  vardf <- get(var)
  if (var==sheetlist[1]){
    delta_state <- vardf
  } else {
    delta_state <- merge(delta_state, vardf, by='Time')
  } # end if var==sheetlist[1]
} # end for var in sheetlist
names(delta_state) <- append(c('Time'), sheetnames)

experiment <- 'perturbhist'

# load sacramento river flow alone:
# sac_only <- read.csv('../../data/dsm2_rsac155_flow_histdss.csv', skip=4, header=FALSE)
sac_only <- read.csv('../../data/dsm2_rsac155_flow_histdss.csv', header=FALSE)
sac_only$Time <- lubridate::ymd(sac_only$V1) + lubridate::days(1)
sac_only <- sac_only[,c('Time','V2')]
names(sac_only) <- c('Time','Sacramento')

delta_state <- merge(delta_state, sac_only, by='Time')
delta_state$NF_nonSac <- delta_state$`Northern Flow`-delta_state$Sacramento

# Set parameters ----------------------------------------------------------

scale_denom <- 2476156.199 # through some testing this seems to be the scaling ratio where the first number is the approximate increase in flow on the peaks
pulse_params <- hash() # flow peak scaling, random shift mean, random shift standard deviation, 
#                       min criteria, max criteria, number of pulses
#                       apply stochify

# pulse_params[['Northern Flow']] <- c(15000/scale_denom, 200, 1500, 5000, NA)
pulse_params[['Sacramento']] <- c(15000/scale_denom, 200, 1500, 
                                  4000, NA, 2,
                                  TRUE, 1)
pulse_params[['SJR Flow']] <- c(2000/scale_denom, 50, 500, 
                                200, NA, 2,
                                TRUE, .9)
pulse_params[['Exports']] <- c(100/scale_denom, 1000, 2000, 
                               50, max(delta_state$Exports), 2,
                               FALSE, 0)
# pulse_params[['Consump Use']] <- c(200/scale_denom, 1000, 2000, 
#                                     -15000, 5000, 6,
#                                    FALSE)

# Create data variability -------------------------------------------------

runyrs <- seq(2006,2016,1)
runyr <- 2006

for (runyr in runyrs) {
  print(paste0("  - ",runyr,"  ====="))
  delta_df <- delta_state[lubridate::year(delta_state$Time)==runyr,]
  delta_df$`Net Delta Outflow` <- delta_df$NF_nonSac + delta_df$Sacramento + delta_df$`SJR Flow` - 
    delta_df$Exports - delta_df$`Consump Use`
  
  edit.df <- perturb_all(delta_df, pulse_params)
  
  if (runyr==runyrs[1]){
    full.df <- edit.df
  } else {
    full.df <- dplyr::bind_rows(full.df, edit.df)
  }
}


# Net Delta Outflow
full.df$`Net Delta Outflow` <- full.df$NF_nonSac + full.df$Sacramento + full.df$`SJR Flow` - full.df$Exports - 
  delta_state$`Consump Use`[delta_state$Time>=min(full.df$Time) & delta_state$Time<=max(full.df$Time)]
delta_state$`Net Delta Outflow` <- delta_state$NF_nonSac + delta_state$Sacramento + delta_state$`SJR Flow` - delta_state$Exports - 
  delta_state$`Consump Use`


# Plots
plt.vars <- append(keys(pulse_params),c('Time', 'Net Delta Outflow'))
plt.df <- melt(full.df, id.vars='Time')
plt.df$Version <- 'Edited'
plt.df$value <- as.numeric(plt.df$value)
plt.df <- plt.df[plt.df$variable %in% plt.vars,]
og.plt.df <- melt(delta_state[,names(delta_state) %in% plt.vars], id.vars='Time')
og.plt.df$Version <- 'Original'
plt.df <- rbind(plt.df, og.plt.df)

plt <- ggplot(data=plt.df) + 
  geom_line(aes(x=Time, y=value, color=Version)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_x_date(limits=lubridate::ymd(c(paste0(runyrs[1],'-1-1'),
                                       paste0(max(runyrs),'12-31'))))
# plt

pltl <- ggplotly(plt, dynamicTicks=TRUE)
pltl

pltl_name <- "perturb_historical.html"

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))


# Write out results -------------------------------------------------------

csv_dir <- "./data_out/"

for (name in overwrite) {

    out_df <- full.df[,c('Time',name)]
    out_csv <- str_replace_all(name," ","_")
    write.table(out_df, file=paste0(csv_dir,'/',out_csv,'_',min(runyrs),'-',max(runyrs),'_perturb_historical.csv'),
                sep=',', row.names=FALSE, col.names=FALSE)

}
