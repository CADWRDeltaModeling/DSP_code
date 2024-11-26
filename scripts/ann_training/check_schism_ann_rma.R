
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

schism_ann_data_dir <- "D:/projects/delta_salinity/scripts/DSP_code/scripts/ann_training/data_out/schism_lhc_v3/"
rma_ann_data_dir <- "D:/projects/delta_salinity/scripts/casanntra/data/"

# schism_cases <- c('lhc_1','lhc_2','lhc_3','lhc_4','lhc_5','lhc_6','lhc_7','mss')
cases <- seq(1,7)

compare_std <- c("datetime", "model", "scene", "case")
compare_flux <- c("sac_flow", "sjr_flow", "exports", "ndo", "cu_flow")
compare_unitless <- c("sf_tidal_energy", "x2", "dcc", "smscg")
compare_ec <- c('emm2','rsl','jer','mal','vcu','wci','bdl','cse','anh') #,'tms','anh','cll')
#c("trp", "wci", "vcu", "uni", "rsl", "old", "pct", "mal", "cll", "emm2", "srv", "anc", "jer", "sal", "ppt", "rri2", "bdt", "lps", "snc", "dsj", "bdl", "nsl2", "vol", "tss", "sss", "tms", "anh", "oh4", "rsl.1", "vcu.1", "mtz")
all_vars <- c(compare_std, compare_flux, compare_unitless, compare_ec)
# Read in data ------------------------------------------------------------

for (case in cases) {
  schism.dat.in <- read.csv(paste0(schism_ann_data_dir,'calsim_ann_baseline_lhc_',case,'.csv'))
  schism.dat.in <- schism.dat.in[,all_vars]
  schism.dat.in$datetime <- lubridate::ymd(schism.dat.in$datetime)
  
  rma.dat.in <- read.csv(paste0(rma_ann_data_dir,'rma_base_',case,'.csv'))
  rma.dat.in <- rma.dat.in[,all_vars]
  rma.dat.in$datetime <- lubridate::ymd(rma.dat.in$datetime)
  
  case.df <- rbind(schism.dat.in, rma.dat.in)
  case.df$case <- case
  
  if (case==cases[1]){
    comb.df <- case.df
  } else {
    comb.df <- rbind(comb.df, case.df)
  }
}
comb.df <- melt(comb.df, id.vars=compare_std)

# plot outputs  ---------------------------------------------------

for (case in cases) {
  
  # plot
  plt.df <- comb.df[comb.df$case == case, ]
  
  plt <- ggplot(data=plt.df) + 
    geom_line(aes(x=datetime, y=value, color=model)) +
    facet_wrap(~variable, ncol=2, scales='free_y') +
    scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
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
  
  pltl <- ggplotly(plt, dynamicTicks=TRUE)
  
  pltl_name <- paste0("case_",case,".html")
  
  saveWidget(pltl, pltl_name, selfcontained=TRUE)
  file.rename(pltl_name, paste0("plots/check_ann_inputs/",pltl_name))
  
} 

