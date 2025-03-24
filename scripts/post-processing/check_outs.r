
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


casanntra.dir <- '../../../casanntra/data'

# Load obs data -----------------------------------------------------------

obs_stas<- c("des_bdl", "des_cse", "ncro_emm2","usbr_jer", "dwr_rsl")
obs_label <- c("bdl", "cse", "emm2", "jer", "rsl")

for (s in seq(1,length(obs_stas))) {
  df <- read.csv(paste0("./data_out/",obs_stas[s],"_obs_ec_1990_2025.csv"))
  df$datetime <- lubridate::ymd_hms(df$datetime)
  df <- df %>%
    mutate(datetime = as.Date(datetime)) %>%  # Extract the date
    group_by(datetime) %>%                    # Group by date
    summarise(!!sym(obs_label[s]) := mean(value, na.rm = TRUE))
  df <- melt(df, id.vars=c("datetime"))
  
  if (s==1) {
    obs_df <- df
  } else {
    obs_df <- rbind(obs_df, df)
  }
}

# Load data ---------------------------------------------------------------

# data.dir <- '../ann_training/data_out'

cases <- c(seq(1,107), seq(1001,1007))

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


# Plot NDO ----------------------------------------------------------------


plt <- ggplot(x2_ndo_df, aes(x=ndo, y=x2, color=case)) +
  geom_point()

plt

ggplotly(plt)


# check specific cases ----------------------------------------------------

# cases in the 90s look like there's an anomoly for high ndo and high x2 (expected inverse relationship)
for (case in seq(91,100)) {
  df <- read.csv(paste0(casanntra.dir,'/dsm2_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  
  if (case==91) {
    dsm2_check_df <- df
  } else {
    dsm2_check_df <- rbind(dsm2_check_df, df)
  }
}
plt_cases <- dsm2_check_df[dsm2_check_df$case %in% seq(91,100),c('datetime','case',"mrz_tidal_energy", "mrz_tidal_filter", 
                                                                 'sf_tidal_energy','sf_tidal_filter')]
plt_cases <- melt(plt_cases, id.vars=c('datetime','case'))
plt_cases$case <- as.factor(plt_cases$case)

plt <- ggplot(plt_cases, aes(x=datetime, y=value, color=case)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  geom_line()

# plt

pltl <- ggplotly(plt, dynamicTicks=TRUE)

pltl_name <- paste0("check_dsm2_91-100.html")

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))
# the anomoly was shown in the first 30 days. This is not model spinup....


# Plot a few EC out locs --------------------------------------------------


# cases <- c(1,11,21,31,41,51,61,71,81,91)
# cases <- c(23,45,seq(101,107))
# cases <- c(23,seq(101,104))
# cases <- c(45,seq(105,107))
# good test for smcg
cases <- c(5,6,9, # 6 & 9 same except for smscg ops and 5 is diff from both
           11,12,13,
           21,22,23,
           37,38,39,
           45,46,47,
           53,54,55,
           64,65,66,
           76,78,79,
           82,83,84,
           96,97,98,
           seq(1001,1007))  
color_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A",
  "#FF7F00", "#FFFF33", "#A65628",
  "#984EA3", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F",
  "#B3B3B3", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E",
  "#A6761D", "#666666", "#17BECF",
  "#BCBD22", "#D62728", "#9467BD",
  "#2CA02C", "#FF7F0E", "#1F77B4",
  "#E41A1C", "#377EB8", "#4DAF4A",
  "#FF7F00", "#FFFF33", "#A65628",
  "#984EA3", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB"
)
for (case in cases) {
  df <- read.csv(paste0(casanntra.dir,'/dsm2_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df <- df[,c('datetime','case',"mrz_tidal_filter",'dcc','northern_flow','cu_flow','smscg','emm2')] # subset variables
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

plt <- ggplot(var_df, aes(x=datetime, y=value, color=case)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_color_manual(values = color_palette) +
  geom_line()
# plt

ggplotly(plt, dynamicTicks=TRUE)

# plt <- ggplot(var_df[var_df$variable %in% c('sf_tidal_energy','sf_tidal_filter',"mrz_tidal_energy", "mrz_tidal_filter"),], aes(x=datetime, y=value, color=variable)) +
#   geom_line()
# 
# ggplotly(plt)



# Compare SCHISM to DSM2 --------------------------------------------------

casanntra.dir <- '../../../casanntra/data'

# SCHISM
sch_cases <- seq(1,7)
# sch_cases <- c(1)

for (case in sch_cases) {
  df <- read.csv(paste0(casanntra.dir,'/schism_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df$case <- case + 1000
  
  if (case==1) {
    sch_df <- df
  } else {
    sch_df <- rbind(sch_df, df)
  }
}

# SCHISM SLR

for (case in sch_cases) {
  df <- read.csv(paste0(casanntra.dir,'/schism_slr_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df$case <- case + 1000
  
  if (case==1) {
    sch_slr_df <- df
  } else {
    sch_slr_df <- rbind(sch_slr_df, df)
  }
}

# RMA
rma_cases <- seq(1,7)
# sch_cases <- c(1)

for (case in rma_cases) {
  df <- read.csv(paste0(casanntra.dir,'/rma_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df$case <- case + 1000
  
  if (case==1) {
    rma_df <- df
  } else {
    rma_df <- rbind(rma_df, df)
  }
}

# DSM2
dsm2_cases <- seq(1001,1007)
# dsm2_cases <- c(1001)

for (case in dsm2_cases) {
  df <- read.csv(paste0(casanntra.dir,'/dsm2_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  
  if (case==1001) {
    dsm2_df <- df
  } else {
    dsm2_df <- rbind(dsm2_df, df)
  }
}

compare_std <- c("datetime", "model", "scene", "case")
compare_flux <- c("sac_flow", "sjr_flow", "exports", "cu_flow", "ndo")
# compare_flux <- c("northern_flow", "sjr_flow", "exports", "cu_flow", "ndo")
compare_tide <- c("sf_tidal_energy","sf_tidal_filter")
compare_x2 <- c("x2")
compare_gates <- c("smscg","dcc")
compare_ec <- c('bdl','emm2','rsl','jer','cse') #,'anh','mal','vcu','wci','tms','anh','cll')
#c("trp", "wci", "vcu", "uni", "rsl", "old", "pct", "mal", "cll", "emm2", "srv", "anc", "jer", "sal", "ppt", "rri2", "bdt", "lps", "snc", "dsj", "bdl", "nsl2", "vol", "tss", "sss", "tms", "anh", "oh4", "rsl.1", "vcu.1", "mtz")
all_vars <- c(compare_std, compare_flux, compare_tide, 
              compare_x2, compare_gates, compare_ec)
plt_vars <- c(compare_flux, compare_tide, 
              compare_x2, compare_gates, compare_ec)
plt_ec_vars <- c('ndo',compare_ec)

dsm2_df <- dsm2_df[,names(dsm2_df) %in% all_vars]
sch_df <- sch_df[,names(sch_df) %in% all_vars]
rma_df <- rma_df[,names(rma_df) %in% all_vars]
sch_slr_df <- sch_slr_df[,names(sch_slr_df) %in% all_vars]
sch_slr_df$model <- "SLR-SCHISM"
# setdiff(names(rma_df), names(sch_df))
comb.df <- rbind(dsm2_df, sch_df)
comb.df <- rbind(comb.df, rma_df)
comb.df <- rbind(comb.df, sch_slr_df)
comb.df <- melt(comb.df, id.vars=compare_std)
comb.df$variable <- factor(comb.df$variable, levels=plt_vars)

plt.obs.df <- obs_df
plt.obs.df$model <- 'Observed'
plt.obs.df$scene <- 'base'
plt.obs.df$case <- NA

compare_cases <- seq(1001,1007)

for (case in compare_cases) {
  
  # plot
  plt.df <- comb.df[comb.df$case == case, ]
  plt.df <- rbind(plt.df, plt.obs.df[plt.obs.df$datetime >= 
                                       min(plt.df$datetime) & 
                                       plt.obs.df$datetime <= 
                                       max(plt.df$datetime), ])
  plt.df$model <- factor(plt.df$model, levels=c('Observed','dsm2','SCHISM','SLR-SCHISM','RMA'))
  
  plt <- ggplot(data=plt.df) + 
    geom_line(aes(x=datetime, y=value, color=model)) +
    facet_wrap(~variable, ncol=2, scales='free_y') +
    scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
    scale_color_manual(values=c("black", "#377EB8", "#E41A1C", "#E7298A", "#66A61E"),
                       breaks=c("Observed","dsm2","SCHISM",'SLR-SCHISM',"RMA")) +
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
  
  pltl_name <- paste0("check_dsm2_v_schism_",case,".html")
  
  saveWidget(pltl, pltl_name, selfcontained=TRUE)
  file.rename(pltl_name, paste0("plots/",pltl_name))
  
  # plot ec out only
  plt.ec.df <- plt.df[plt.df$variable %in% plt_ec_vars,]
  plt <- ggplot(data=plt.ec.df) + 
    geom_line(aes(x=datetime, y=value, color=model)) +
    facet_wrap(~variable, ncol=1, scales='free_y') +
    scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
    scale_color_manual(values=c("black", "#377EB8", "#E41A1C", "#E7298A", "#66A61E"),
                       breaks=c("Observed","dsm2","SCHISM",'SLR-SCHISM',"RMA")) +
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
  
  pltl_name <- paste0("ec_check_dsm2_v_schism_",case,".html")
  
  saveWidget(pltl, pltl_name, selfcontained=TRUE)
  file.rename(pltl_name, paste0("plots/",pltl_name))
} 



# Check National Steel ----------------------------------------------------

sch_out_4_df <- read.csv("./data_out/schism_lhc_v3/baseline_lhc_4.csv")
sch_out_5_df <- read.csv("./data_out/schism_lhc_v3/baseline_lhc_5.csv")

nsl_4_df <- sch_out_4_df %>%
  mutate(
    datetime = as.Date(X0),  # Convert datetime to date format
    model = "SCHISM",
    scene = "base",
    case = 1004,
    variable = "nsl_cfs"
  ) %>%
  group_by(datetime, model, scene, case, variable) %>%  # Group by datetime & constants
  summarise(value = mean(nsl.Out..cms., na.rm = TRUE) * 35.31, .groups = "drop")  # Compute daily mean & convert to CFS

nsl_5_df <- sch_out_5_df %>%
  mutate(
    datetime = as.Date(X0),  # Convert datetime to date format
    model = "SCHISM",
    scene = "base",
    case = 1005,
    variable = "nsl_cfs"
  ) %>%
  group_by(datetime, model, scene, case, variable) %>%  # Group by datetime & constants
  summarise(value = mean(nsl.Out..cms., na.rm = TRUE) * 35.31, .groups = "drop")  # Compute daily mean & convert to CFS
nsl_df <- rbind(nsl_4_df, nsl_5_df)

plt.obs.nsl <- plt.obs.df %>%
  filter(variable == "bdl" &
           datetime >= min(plt_nsl$datetime) & datetime <= max(plt_nsl$datetime)) %>%
  mutate(case = "Observed")

plt_nsl <- sch_df  %>%
  select(all_of(compare_std), smscg, bdl) %>%
  melt(id.vars = compare_std) %>%
  filter(case %in% c(1004, 1005))
plt_nsl <- rbind(plt_nsl, nsl_df)
plt_nsl$case <- factor(plt_nsl$case)
plt_nsl <- plt_nsl[complete.cases(plt_nsl),]
plt_nsl <- plt_nsl[plt_nsl$datetime < max(plt_nsl$datetime),]
plt_nsl <- rbind(plt_nsl, plt.obs.nsl)

plt <- ggplot(data=plt_nsl) + 
  geom_line(aes(x=datetime, y=value, color=case)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
  scale_color_manual(values=c("black", "#377EB8", "#E41A1C"),
                     breaks=c("Observed","1004","1005")) +
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
pltl

pltl_name <- paste0("nsl_check_schism_4-5.html")

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))


# Check SLR ---------------------------------------------------------------

# SCHISM
sch_cases <- seq(1,7)
# sch_cases <- c(1)

for (case in sch_cases) {
  df <- read.csv(paste0(casanntra.dir,'/schism_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df$case <- case + 1000
  
  if (case==1) {
    sch_df <- df
  } else {
    sch_df <- rbind(sch_df, df)
  }
}

# SCHISM SLR

for (case in sch_cases) {
  df <- read.csv(paste0(casanntra.dir,'/schism_slr_base_',case,'.csv'))
  df$datetime <- lubridate::ymd(df$datetime)
  df$case <- case + 1000
  
  if (case==1) {
    sch_slr_df <- df
  } else {
    sch_slr_df <- rbind(sch_slr_df, df)
  }
}


compare_std <- c("datetime", "model", "scene", "case")
compare_tide <- c("sf_tidal_energy","sf_tidal_filter")


tide_sch_df <- sch_df[,names(sch_df) %in% c(compare_std, compare_tide)]
tide_sch_slr_df <- sch_slr_df[,names(sch_slr_df) %in% c(compare_std, compare_tide)]
tide_sch_slr_df$model <- "SLR-SCHISM"

comb.df <- rbind(tide_sch_df, tide_sch_slr_df)
comb.df <- melt(comb.df, id.vars=compare_std)
comb.df$variable <- factor(comb.df$variable, levels=compare_tide)

case <- 1004

plt.df <- comb.df[comb.df$case == case, ]
plt.df$model <- factor(plt.df$model, levels=c('Observed','dsm2','SCHISM','SLR-SCHISM','RMA'))

plt <- ggplot(data=plt.df) + 
  geom_line(aes(x=datetime, y=value, color=model)) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
  scale_color_manual(values=c("black", "#377EB8", "#E41A1C", "#E7298A", "#66A61E"),
                     breaks=c("Observed","dsm2","SCHISM",'SLR-SCHISM',"RMA")) +
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
pltl
