# Cleaning up cluster_by_wy, this is the code used as of 5/19 for year selection criteria for ANN


# Set up --------------------------------------------------

# load libraries
library(lubridate)
library(readxl)
library(ggplot2)
library(plotly)
library(reshape2)
library(stringr)
library(htmlwidgets)
library(scales)

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory

# read in data ------------------------------------------------------------

cases <- c('lhc_3','lhc_4','lhc_5')
case_inputs <- paste0('./input/vsource',c('_case3','_case4','.lhc_5'),'.th')
ref_dates <- c('2013-10-24','2006-11-14','2006-11-14')

ss_names <- read.table('./input/source_sink.in',sep="!",header=FALSE,nrows=1)
ss_names <- read.table('./input/source_sink.in',sep="!",header=FALSE,nrows=ss_names$V1[1],skip=1)
ss_names$V2 <- trimws(ss_names$V2)
ss_names$V2 <- sub(",.*$", "", ss_names$V2)

rsl_ss <- read.table('./input/rsl_area_source_sinks.csv',header=TRUE)
rsl_src <- rsl_ss[rsl_ss$stype=='source',]
rsl_src$ss_index <- match(rsl_src$site, ss_names$V2)

# Load th files -----------------------------------------------------------

for (i in seq_along(cases)) {
  th_in <- read.table(case_inputs[i], header=FALSE, sep="")
  th_in$datetime <- lubridate::ymd(ref_dates[i]) + seconds(th_in$V1)
  th_in <- th_in[, c(which(names(th_in) == "datetime"), rsl_src$ss_index+1)]
  names(th_in) <- c("datetime",rsl_src$site)
  th_in <- melt(th_in, id.vars='datetime')
  th_in$case <- cases[i]
  
  if (i==1){
    src_sink_plt_df <- th_in
  } else {
    src_sink_plt_df <- rbind(src_sink_plt_df, th_in)
  }
}

rm(th_in)

# Plot --------------------------------------------------------

plt <- ggplot() +
  geom_line(src_sink_plt_df, mapping=aes(x=datetime, y=value, color=case)) +
  facet_wrap(~variable, ncol=2)

gplt <- ggplotly(plt, tooltip=c('datetime','value'), dynamicTicks=TRUE)
gplt <- gplt %>%
  layout(hovermode = "x unified")

gplt
