# Libraries commonly used -------------------------------------------------

library(lubridate)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory


casanntra.dir <- '../../../scripts/casanntra/data'

compare_std <- c("datetime", "model", "scene", "case")
compare_flux <- c("sac_flow", "sjr_flow", "exports", "cu_flow", "ndo")
# compare_flux <- c("northern_flow", "sjr_flow", "exports", "cu_flow", "ndo")
compare_tide <- c("sf_tidal_energy","sf_tidal_filter")
compare_x2 <- c("x2")
compare_gates <- c("smscg","dcc")
compare_ec <- c('mrz','bdl','emm2','rsl','jer','cse') #,'anh','mal','vcu','wci','tms','anh','cll')
#c("trp", "wci", "vcu", "uni", "rsl", "old", "pct", "mal", "cll", "emm2", "srv", "anc", "jer", "sal", "ppt", "rri2", "bdt", "lps", "snc", "dsj", "bdl", "nsl2", "vol", "tss", "sss", "tms", "anh", "oh4", "rsl.1", "vcu.1", "mtz")
all_vars <- c(compare_std, compare_flux, compare_tide, 
              compare_x2, compare_gates, compare_ec)
plt_vars <- c(compare_flux, compare_tide, 
              compare_x2, compare_gates, compare_ec)
plt_ec_vars <- c('ndo',compare_ec)

model_colors <- c("black", "#377EB8", "#E41A1C", "#984EA3", "#E7298A", "#66A61E")
model_levels <- c('Observed','dsm2','SCHISM','SCHISM-SUISUN','SCHISM-SLR','RMA')

# Load obs data -----------------------------------------------------------

obs_stas<- c("des_bdl", "des_cse", "ncro_emm2","usbr_jer", "dwr_rsl", "des_mrz")
obs_label <- c("bdl", "cse", "emm2", "jer", "rsl", "mrz")

obs_df <- map2_dfr(obs_stas, obs_label, ~{
  df <- read_csv(paste0("./data_out/", .x, "_obs_ec_1990_2025.csv")) %>%
    mutate(datetime = ymd_hms(datetime),
           datetime = as.Date(datetime)) %>%
    group_by(datetime) %>%
    summarise(!! .y := mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(cols = -datetime, names_to = "variable", values_to = "value")
})

dayflow.csv <- "//nasbdo/Modeling_Data/dayflow/dayflow_1983_2023.csv"
dayflow.df <- read.csv(dayflow.csv, skip=1)
dayflow.df$datetime <- lubridate::ymd(dayflow.df$datetime)


plt.obs.df <- bind_rows(
  obs_df,
  dayflow.df %>%
    transmute(datetime, variable = "ndo", value = OUT)
) %>%
  mutate(model = factor("Observed"),
         scene = "base",
         case = NA)
plt.obs.df$model <- as.factor('Observed')
plt.obs.df$scene <- 'base'
plt.obs.df$case <- NA



# Helper function to load one case ----------------------------------------

read_case <- function(prefix, cases, model = NULL, suffix = "", case_offset = 0) {
  map_dfr(cases, function(case) {
    file <- paste0(casanntra.dir, "/", prefix, case, suffix, ".csv")
    df <- read_csv(file, show_col_types = FALSE) %>%
      mutate(datetime = ymd(datetime),
             case = case + case_offset)
    if (!is.null(model)) {
      df$model <- model
    }
    return(df)
  })
}


# Load all model outputs --------------------------------------------------

sch_cases <- 1:7
dsm2_cases <- 1001:1007

sch_df         <- read_case("schism_base_",      sch_cases, model = "SCHISM",        case_offset = 1000)
sch_slr_df     <- read_case("schism_slr_base_",  sch_cases, model = "SCHISM-SLR",    case_offset = 1000)
sch_suisun_df  <- read_case("schism_suisun_",    sch_cases, model = "SCHISM-SUISUN", case_offset = 1000)
rma_df         <- read_case("rma_base_",         sch_cases, model = "RMA",           case_offset = 1000)
dsm2_df        <- read_case("dsm2_base_",        dsm2_cases)  # already has case numbers

# Keep only desired columns
filter_cols <- function(df) df[, names(df) %in% all_vars]
dsm2_df       <- filter_cols(dsm2_df)
sch_df        <- filter_cols(sch_df)
rma_df        <- filter_cols(rma_df)
sch_slr_df    <- filter_cols(sch_slr_df)
sch_suisun_df <- filter_cols(sch_suisun_df)

# Combine and reshape
comb.df <- bind_rows(dsm2_df, sch_df, rma_df, sch_slr_df, sch_suisun_df) %>%
  pivot_longer(
    cols = -all_of(compare_std),      # melt everything except datetime, model, scene, case
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(variable, levels = plt_vars),
    model = factor(model, levels = model_levels)
  )


# Plot --------------------------------------------------------------------

compare_cases <- seq(1001,1007)

case <- 1001

for (case in compare_cases) {
  
  print(paste0("Plotting ",case))
  
  plt.df <- comb.df[comb.df$case == case, ]
  plt.df <- rbind(plt.df, plt.obs.df[plt.obs.df$datetime >=
                                       min(plt.df$datetime) &
                                       plt.obs.df$datetime <=
                                       max(plt.df$datetime), ])
  
  # plot ec out only
  plt.ec.df <- plt.df[plt.df$variable %in% plt_ec_vars,]
  plt <- ggplot(data=plt.ec.df) + 
    geom_line(aes(x=datetime, y=value, color=model)) +
    facet_wrap(~variable, ncol=1, scales='free_y') +
    scale_x_date(date_breaks = "6 month", date_labels = "%b-%Y", expand=c(0,0)) +
    xlab("") +
    ylab("EC (psu)") +
    scale_color_manual(values=model_colors,
                       breaks=model_levels,
                       drop=FALSE) +
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
  
  plt
  pltl <- ggplotly(plt, dynamicTicks=TRUE)
  
  pltl_name <- paste0("ec_check_dsm2_v_schism_",case,".html")
  
  saveWidget(pltl, pltl_name, selfcontained=TRUE)
  file.rename(pltl_name, paste0("plots/",pltl_name))
} 
