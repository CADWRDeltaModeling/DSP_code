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
# library(stochCycle)
# doingStochCycle=TRUE

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
source("./functions/load_delta_vars.R")
source("./functions/pulse_events.R")

# Load data ---------------------------------------------------------------

# load delta state vars
dsm2_filename <- '../../data/dsm2_ann_inputs_base.xlsx'

sheetlist <- c("northern_flow","sjr_flow","exports","dxc_gate_fraction",
               "net_delta_cu","mtz_daily_max-min_stage",
               "sjr_vernalis_ec","sac_ec","base_ec_output")
sheetnames <- c('Northern Flow', 'SJR Flow','Exports','DCC Gate',
                'Consump. Use','Tidal Amp', 
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

# Set parameters ----------------------------------------------------------

scale_denom <- 2476156.199 # through some testing this seems to be the scaling ratio where the first number is the approximate increase in flow on the peaks
pulse_params <- hash() # flow peak scaling, random shift mean, random shift standard deviation, min criteria, max criteria

pulse_params[['Northern Flow']] <- c(15000/scale_denom, 200, 1500, 5000, NA)
pulse_params[['SJR Flow']] <- c(2000/scale_denom, 50, 500, 200, NA)
pulse_params[['Exports']] <- c(100/scale_denom, 1000, 2000, 50, max(delta_state$Exports))
pulse_params[['DCC Gate']] <- c(1/scale_denom, 0, 1.5, NA, NA)

# Create data variability -------------------------------------------------

years <- seq(2006,2016,1)

for (year in years) {
  delta_df <- delta_state[lubridate::year(delta_state$Time)==year,]
  delta_df$`Net Delta Outflow` <- delta_df$`Northern Flow` + delta_df$`SJR Flow` - 
    delta_df$Exports - delta_df$`Consump. Use`
  
  edit.df <- perturb_all(delta_df, pulse_params)
  
  if (year==years[1]){
    full.df <- edit.df
  } else {
    full.df <- dplyr::bind_rows(full.df, edit.df)
  }
}


# Fix DCC
full.df$`DCC Gate`[full.df$`DCC Gate`>=0] <- 1
full.df$`DCC Gate`[full.df$`DCC Gate`<0] <- 0

# Net Delta Outflow
full.df$`Net Delta Outflow` <- full.df$`Northern Flow` + full.df$`SJR Flow` - full.df$Exports - 
  delta_state$`Consump. Use`[delta_state$Time>=min(full.df$Time) & delta_state$Time<=max(full.df$Time)]
delta_state$`Net Delta Outflow` <- delta_state$`Northern Flow` + delta_state$`SJR Flow` - delta_state$Exports - 
  delta_state$`Consump. Use`

# Shift tide data
shift.tide <- delta_state[,c('Time','Tidal Amp')]
shift.tide$Time <- lubridate::days(100) + shift.tide$Time # tidal cycle ~27.3 days
full.df <- merge(full.df, shift.tide, by='Time') 


# Plots
plt.df <- melt(full.df, id.vars='Time')
plt.df$Version <- 'Edited'
plt.df$value <- as.numeric(plt.df$value)
og.plt.df <- melt(delta_state[,names(delta_state) %in% append(keys(pulse_params),c('Time', 'Net Delta Outflow', 'Tidal Amp'))], id.vars='Time')
og.plt.df$Version <- 'Original'
plt.df <- rbind(plt.df, og.plt.df)

plt <- ggplot(data=plt.df) + 
  geom_line(aes(x=Time, y=value, color=Version)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_x_date(limits=lubridate::ymd(c(paste0(years[1],'-1-1'),
                                       paste0(max(years),'12-31'))))

# plt
pltl <- ggplotly(plt, dynamicTicks=TRUE)

pltl_name <- "perturb_historical.html"

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))


# Write out results -------------------------------------------------------

csv_dir <- "./data_out/"

for (name in names(edit.df)[!(names(edit.df) %in% c('Time','Net Delta Outflow'))]) {
  out_df <- edit.df[,names(edit.df) %in% c('Time',name)]
  out_csv <- str_replace_all(name," ","_")
  write.table(out_df, file=paste0(csv_dir,'/',out_csv,'_',min(years),'-',max(years),'_perturb_historical.csv'), 
              sep=',', row.names=FALSE, col.names=FALSE)
}
