## This script is to check the Net Delta Outflow of a simulation that has 
## already run using the boundaries generated by the previous scripts


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
library(scales)
library(yaml)
library(RColorBrewer)
library(tidyr)

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
source("./functions/load_delta_vars.R")

# Define variables ---------------------------------------------------------

dsm2.dir <- '../../model/dsm2/DSP_DSM2_202307/'

sheetlist <- c("northern_flow","sjr_flow","exports","dxc_gate_fraction",
               "suisun_gate_fraction","net_delta_cu","mtz_tidal_nrg",
               "base_ec_output")
sheetnames <- c('Northern Flow', 'SJR Flow','Exports','DCC Gate',
                'Suisun Gate','Consump Use','Tidal Nrg')

# ``` load historical data ----------------------------------------------------

dsm2_filename <- paste0(dsm2.dir,'historical/anninputs/dsm2_ann_inputs_historical.xlsx')

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

ds.df <- melt(delta_state, id.vars='Time')
ds.df$case <- 'hist'

# ``` load markov perturb data ----------------------------------------------------

datout.dir <- './data_out/suisun_gates_lhc_v4/'

mark1_filename <- paste0(datout.dir,'/MTZSL_markov_pert_v1.csv')
mark2_filename <- paste0(datout.dir,'/MTZSL_markov_pert_v2.csv')

colnames <- c('Time',"radial_op", "radial_up", "radial_down", "flashboard", "boat_lock")

mark1 <- read.csv(mark1_filename, col.names = colnames)
mark2 <- read.csv(mark2_filename, col.names = colnames)

mark1$Time <- lubridate::ymd_hms(mark1$Time)
mark1$radial_op[mark1$radial_op==1] <- 0
mark1$radial_op[mark1$radial_op==-10] <- 1
mark2$Time <- lubridate::ymd_hms(mark2$Time)
mark2$radial_op[mark2$radial_op==1] <- 0
mark2$radial_op[mark2$radial_op==-10] <- 1


# plot -------------------

suis.df <- merge(delta_state[,c('Time','Suisun Gate')], mark1[,c('Time','radial_op')], by='Time', all=TRUE)
suis.df <- merge(suis.df, mark2[,c('Time','radial_op')], by='Time', all=TRUE)
colnames(suis.df) <- c('Time','Historical','Markov1','Markov2')
suis.df <- suis.df %>% fill(Markov1)
suis.df <- suis.df %>% fill(Markov2)
suis.df <- suis.df[complete.cases(suis.df),]

plt.df <- melt(suis.df, id.vars='Time')

plt <- ggplot(plt.df, aes(x=Time,y=value,color=variable)) +
  facet_wrap(~variable, ncol=1) +
  xlab("") +
  ylab("Operation") +
  geom_line()

plt
ggplotly(plt)



# check ala test_smscg.py -------------------------------------------------

mark1$boatcheck <- (mark1$radial_up == 0) & (mark1$boat_lock == 1) |
  (mark1$radial_up == 1)


