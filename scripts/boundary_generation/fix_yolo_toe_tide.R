

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

lhc_4_flux_fn <- './input/flux.th'
lhc_5_flux_fn <- './input/flux.lhc_5.th'
original_flux_fn <- "D:/schism/repositories/BayDeltaSCHISM/data/time_history/flux.th"
flux_head_fn <- './input/flux.lhc_5.dated.th'
lhc_4_monterey_fn <- './input/monterey.csv' 
lhc_5_monterey_fn <- './input/monterey_shift_forward.csv' 

start_date <- lubridate::ymd_hm('2006-11-14 00:00')

# Read in data ------------------------------------------------------------

flux_4.df <- read.table(lhc_4_flux_fn)
flux_5.df <- read.table(lhc_5_flux_fn)
original.df <- read.table(original_flux_fn, header=TRUE)
names(flux_4.df) <- names(read.table(flux_head_fn, nrows=1, header=TRUE))
names(flux_5.df) <- names(read.table(flux_head_fn, nrows=1, header=TRUE))
original.df$datetime <- lubridate::ymd_hm(original.df$datetime) 
flux_4.df$datetime <- as.POSIXct(start_date) + as.numeric(flux_4.df$datetime)
flux_5.df$datetime <- as.POSIXct(start_date) + as.numeric(flux_5.df$datetime)

mont_4.df <- read.delim(lhc_4_monterey_fn, sep=',', header=TRUE, comment.char='#')
mont_5.df <- read.delim(lhc_5_monterey_fn, sep=',', header=TRUE, comment.char='#')
mont_4.df$Date.Time <- lubridate::ymd_hm(mont_4.df$Date.Time)
mont_5.df$Date.Time <- lubridate::ymd_hms(mont_5.df$Date.Time)


# Refine & combine data ---------------------------------------------------

toedrain_4 <- flux_4.df[,c('datetime','yolo_toedrain')]
toedrain_5 <- flux_5.df[,c('datetime','yolo_toedrain')]
toedrain_orig <- original.df[,c('datetime','yolo_toedrain')]
tide_4 <- mont_4.df[,c('Date.Time', 'Water.Level')]
tide_5 <- mont_5.df[,c('Date.Time','Water.Level')]
names(tide_4) <- c('datetime','monterey_wse_m')
names(tide_5) <- c('datetime','monterey_wse_m')

toedrain_4 <- melt(toedrain_4,id.vars='datetime')
toedrain_4$case <- 'LHC 4'
toedrain_5 <- melt(toedrain_5,id.vars='datetime')
toedrain_5$case <- 'LHC 5'
toedrain_orig <- melt(toedrain_orig,id.vars='datetime')
toedrain_orig$case <- 'Full'

tide_4 <- melt(tide_4, id.vars='datetime')
tide_4$case <- 'LHC 5'
tide_5 <- melt(tide_5, id.vars='datetime')
tide_5$case <- 'LHC 5'

toe_df <- rbind(toedrain_4, toedrain_5)
toe_df <- rbind(toe_df, toedrain_orig)
tide_df <- rbind(tide_4, tide_5)


# facet_wrap mixed log scale ----------------------------------------------

plt <- ggplot(data=toe_df) + 
  geom_line(aes(x=datetime, y=value, color=case)) +
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

ggplotly(plt)


# Filter out flood pulses -------------------------------------------------

# Shift everything "forward" by 100 days

toedrain_edit <- original.df[,c('datetime','yolo_toedrain')]
toedrain_edit$datetime <- toedrain_edit$datetime + lubridate::as.difftime(2.5, units='hours')
toedrain_edit <- merge(original.df[,c('datetime','yolo_toedrain')], toedrain_edit, 
                       by='datetime',
                       all=TRUE)
names(toedrain_edit) <- c('datetime','Original','2.5h Shift')
toedrain_100d <- original.df[,c('datetime','yolo_toedrain')]
toedrain_100d$datetime <- toedrain_100d$datetime + lubridate::as.difftime(100, units='days')
names(toedrain_100d) <- c('datetime','100d Shift')

toedrain_edit <- merge(toedrain_edit, toedrain_100d, by='datetime', all=TRUE)

# Filter flood pulses
toedrain_edit$Edit <- toedrain_edit$`2.5h Shift`
toedrain_edit$Edit[(toedrain_edit$Original < -19 | toedrain_edit$`2.5h Shift` < -19) & 
                     !is.na(toedrain_edit$Original) & !is.na(toedrain_edit$`2.5h Shift`)] <- 
  toedrain_edit$Original[(toedrain_edit$Original < -19 | toedrain_edit$`2.5h Shift` < -19) & 
                           !is.na(toedrain_edit$Original) & !is.na(toedrain_edit$`2.5h Shift`)]

toe_plt <- melt(toedrain_edit, id.vars='datetime')
toe_plt <- toe_plt[toe_plt$datetime<lubridate::ymd_hm('2009-01-01 00:00'),]

plt <- ggplot(data=toe_plt) + 
  geom_line(aes(x=datetime, y=value, color=variable)) +
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

ggplotly(plt)

flux_5.edit.df <- merge(flux_5.df, toedrain_edit[,c('datetime','Edit')], by='datetime')
flux_5.edit.df$yolo_toedrain <- flux_5.edit.df$Edit
flux_5.edit.df$datetime <-  as.numeric(flux_5.edit.df$datetime - as.POSIXct(start_date))
flux_5.edit.df <- flux_5.edit.df[,names(flux_5.df)]

flux_5.edit.df[] <- lapply(flux_5.edit.df, function(x) if(is.numeric(x)) sprintf("%.2f", x) else x)

write.table(flux_5.edit.df, file='./data_out/flux_lhc_5.edit_toedrain.th', 
            col.names=FALSE, row.names=FALSE, sep=' ', quote=FALSE)
