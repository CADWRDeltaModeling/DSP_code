# Cleaning up cluster_by_wy, this is the code used as of 5/19 for year selection criteria for ANN


# Set up --------------------------------------------------

# load libraries
library(zoo)
library(dplyr)
library(readxl)
library(ggplot2)
library(plotly)
library(corrplot)
library(cluster)
library(rgl)
library(fpc)
library(reshape2)
library(hash)
library(stringr)
library(ggh4x)
library(htmlwidgets)
library(scales)
library(colorspace)
library(grid)
library(gtable)

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
source('./functions/compare_clusters_fxn.R')
source("./functions/whipple_adapt_clustergram.R")
source("./functions/load_delta_vars.R")

# read in data ------------------------------------------------------------

# DSM2 Inputs/Outputs for ANN
dsm2_filename <- '../../data/dsm2_ann_inputs_base.xlsx'
# ec outputs are reported in micro-Siemens/cm. 1000 mS/cm = 1 mmhos/cm
# sheetlist <- excel_sheets(path=dsm2_filename)
northern_flow <- load_delta_vars(dsm2_filename, 'northern_flow')
exports <- load_delta_vars(dsm2_filename, 'exports')
ec_outputs <- load_delta_vars(dsm2_filename, 'base_ec_output')

cols <- c("navy","royalblue1","mediumorchid","orangered3","palegreen3","goldenrod","green4","paleturquoise1","yellow1","pink1","slateblue2")
col.clust <- data.frame(clust=c(1:length(cols)), color=cols)
# library(scales)
# show_col(cols)

may1 <- lubridate::ymd("2001/05/01")
jun1 <- lubridate::ymd("2001/06/01")
nov30 <- lubridate::ymd("2001/11/30")

flow.brks <- c(5000,10000, 20000, 30000, 40000, 50000,100000)
y.scale <- facetted_pos_scales(y=list(
  scale_y_continuous(limits=c(0,14000), expand=c(0,0)),
  scale_y_log10(breaks=flow.brks, labels=comma(flow.brks, accuracy=1))))

plt.out.dir <- 'plots'

ec_eval_locs <- c("RSAN018-JERSEYPOINT","RSAC081-COLLINSVILLE","RSAC092-EMMATON","RSAN007-ANTIOCH","ROLD024-OLD RIVER AT BACON ISLAND","CHSWP003-CCFB_INTAKE")

# Pre-prepare Data --------------------------------------------------------

refdate <- '2000/11/30'
# reshape exports
exp.df.rshp <- reshape.df(exports, monthstart=12)
exp.dec.dry.rshp <- exp.df.rshp[,names(exp.df.rshp) %in% c('wy',as.numeric(may1-lubridate::ymd(refdate)):as.numeric(nov30-lubridate::ymd(refdate)))]

# reshape northern flow from December
nfl.df.rshp <- reshape.df(northern_flow, monthstart=12) # reshapes dataframe starting julian day from December 1
nfl.dec.dry.rshp <- nfl.df.rshp[,names(nfl.df.rshp) %in% c('wy',as.numeric(may1-lubridate::ymd(refdate)):as.numeric(nov30-lubridate::ymd(refdate)))]

# combine exports and flows
nfex.full.rshp <- merge(nfl.df.rshp, exp.df.rshp, by='wy') # full dataset for plotting
nfex.decdry.rshp <- merge(nfl.dec.dry.rshp, exp.dec.dry.rshp, by='wy') # columns will be named "julian day".x or "julian_day".y, where x=flow, y=exports

# normalize flow and export data so variation is equally weighted across exports and flows
nf.dd.norm <- norm.biv.df(nfex.decdry.rshp, 2)
nf.dd.norm.stat <- nf.dd.norm[,!(names(nf.dd.norm) %in% c('wy'))]


wss <- (nrow(nf.dd.norm.stat)-1)*sum(apply(nf.dd.norm.stat,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(nf.dd.norm.stat,centers=i)$withinss)
x11()
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# Clustergram
png(filename="plots/flowsexports_decdry_normed_clustergram.png", type="cairo", units="in",width=8,height=6,pointsize=10, res=200)
set.seed(20000)
clustergram(as.matrix(nf.dd.norm.stat), k.range=2:10)
dev.off()

# NOTICE: Looks as though 4 clusters would be most appropriate when including May in analysis

# K-means clustering with 4 -----------------------------------------------


nclust <- 4 # Using 4 4 might be most appropriate
bivars <- c('Flow','Export')
nfex.decdry4 <- clust.plt(nclust, nf.dd.norm.stat, nfex.full.rshp, plt.out.dir, refdate=refdate, bivars=bivars, y.scale=y.scale, only.clust.legend=TRUE) # cluster with nfex.dd.norm.stat but plot with nfex.full.rshp (un-log-normalized full time series for plotting)
nfex.decdry4.km <- nfex.decdry4$km.cb

# save(nfex.decdry4.km, file=paste0('output/',plt.out.dir,'/nfexdecdry4km.RData'))

# Determine selected water years ------------------------------------------

slct.df <- nf.dd.norm
bivars <- c('Flow','Export')
obs <- ncol(slct.df) - 1 # number of days times number of bivariates
clusts <- nfex.decdry4.km$partition

wy.metrics <- data.frame(wy=slct.df$wy, clust=clusts, skill=NA, rmse=NA, mmskill=NA, mmrmse=NA)

slct.df$clust <- clusts
jds <- names(slct.df)[!(names(slct.df) %in% c('wy','clust'))]
for (c in sort(unique(clusts))) {

  # get average values for cluster
  mean.clust <- data.frame(jday=jds, mean=NA)
  for (jd in jds) {
    mean.clust$mean[mean.clust$jday==jd] <- mean(slct.df[slct.df$clust==c,names(slct.df) %in% c(jd)])
  } # end jd mean loop

  # compute metrics for each wy within cluster
  for (wy in slct.df$wy[slct.df$clust==c]) {
    wy.df <- data.frame(jday=jds,
                        pred=unlist(slct.df[slct.df$wy==wy, jds]))
    wy.df <- merge(wy.df, mean.clust, by='jday')
    wy.metrics$rmse[wy.metrics$wy==wy] <- sqrt(sum((wy.df$pred-wy.df$mean)^2)/nrow(wy.df))
    wy.metrics$skill[wy.metrics$wy==wy] <- 1-(sum((wy.df$pred-wy.df$mean)^2)/
                                                (sum((abs(wy.df$pred-mean(wy.df$mean))+abs(wy.df$mean-mean(wy.df$mean)))^2)))
  } # end wy loop

  # compute max and min skill/rmse
  wy.metrics$mmskill[wy.metrics$skill==
                       max(wy.metrics$skill[wy.metrics$clust==c & wy.metrics$wy > 2007], na.rm=TRUE)] <- paste0('Max Skill cluster ', c)
  wy.metrics$mmskill[wy.metrics$skill==
                       min(wy.metrics$skill[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Min Skill cluster ', c)
  wy.metrics$mmrmse[wy.metrics$rmse==
                      max(wy.metrics$rmse[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Max RMSE cluster ', c)
  wy.metrics$mmrmse[wy.metrics$rmse==
                      min(wy.metrics$rmse[wy.metrics$clust==c & wy.metrics$wy > 2007], na.rm=TRUE)] <- paste0('Min RMSE cluster ', c)

} # end cluster loop

wy.metrics <- wy.metrics[order(wy.metrics$clust, wy.metrics$rmse),]
wy.metrics[grepl("Min RMSE",wy.metrics$mmrmse),]
wy.slcts <- wy.metrics[grepl("Min RMSE",wy.metrics$mmrmse),'wy']
wy.slcts <- wy.slcts[!(wy.slcts %in% c(2019))] # chose to only include 2017 and not 2019 in cluster 6 to keep limited to 6 years
sort(wy.slcts)

# save(wy.metrics, file=paste0('output/',plt.out.dir,'/wymetrics_cluster.RData'))

df.plt.rshp <- nfex.decdry4$df.plt.rshp
plt.cst.rshp <- nfex.decdry4$plt.cst.rshp
plots <- c('Min Skill','Max Skill','Min RMSE','Max RMSE')

wy.to.clust <- slct.df[,c('wy','clust')]
# save(wy.to.clust, file=paste0('output/',plt.out.dir,'/wytoclust.RData'))

# Plot final selection ----------------------------------------------------

slct.wys <- wy.slcts
slct.wys <- c(2008, 2010, 2012, 2014)

wycol <- merge(wy.to.clust, col.clust, by='clust')
unslcted <- sort(unique(df.plt.rshp$variable)[!(unique(df.plt.rshp$variable) %in% slct.wys)])
c.colordf <- data.frame(breaks=c(unslcted,
                                 wycol$wy[wycol$wy %in% slct.wys]),
                        values=c(rep('grey',length(unslcted)),
                                 wycol$color[wycol$wy %in% slct.wys]))
flow.brks <- c(5000,10000, 20000, 30000, 40000, 50000,100000)
figname <-  paste0(plt.out.dir,'/Final_WY_Selection.png')

df.plt.rshp$variable <- as.character(df.plt.rshp$variable)
other_levels <- setdiff(unique(df.plt.rshp$variable), as.character(slct.wys))
ordered_levels <- c(as.character(slct.wys), sort(other_levels))
df.plt.rshp$variable <- factor(df.plt.rshp$variable, levels = ordered_levels)

unanalyzed.rect <- data.frame(xmin=rep(as.Date("2001-05-01"),2),
                              xmax=rep(as.Date("2001-05-01"),2),
                              ymin=c(min(df.plt.rshp$value[df.plt.rshp$bivar==bivars[1]]),
                                     min(df.plt.rshp$value[df.plt.rshp$bivar==bivars[2]])),
                              ymax=c(max(df.plt.rshp$value[df.plt.rshp$bivar==bivars[1]]),
                                     max(df.plt.rshp$value[df.plt.rshp$bivar==bivars[2]])),
                              bivar=bivars)
# create facet wrap plot for all years with the selected years plotted in respective cluster color
plt <- ggplot() +
  geom_line(data=df.plt.rshp, mapping=aes(x=dt,y=value,color=variable), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
  geom_line(data=df.plt.rshp[df.plt.rshp$variable %in% slct.wys,], mapping=aes(x=dt,y=value,color=variable), linewidth=1) +
  geom_rect(data=unanalyzed.rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = alpha("grey", 0.5), color = NA, inherit.aes = FALSE) +
  facet_wrap(~bivar, nrow=length(bivars), scales='free') +
  scale_color_manual(breaks=as.character(slct.wys),
                     values=setNames(c.colordf$values, as.character(c.colordf$breaks)),
                     # guide=guide_legend(override.aes=list(linetype=0)),
                     name="Cluster") +
  # scale_y_continuous(limits=c(0,500000), breaks=seq(0,500000,100000), expand=c(0,0)) +
  # scale_y_continuous(limits=c(0,25000), expand=c(0,0)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) +
  labs(x='',y='Discharge (cfs)') +
  guides(color=guide_legend(nrow=1)) +
  facetted_pos_scales(y=list(
    scale_y_continuous(limits=c(0,14000), expand=c(0,0)),
    scale_y_log10(breaks=flow.brks, labels=comma(flow.brks, accuracy=1)))) +
  theme(text=element_text(size=12),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth=0.5),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        # legend.key=element_blank(),
        panel.grid.major.y = element_line(size=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        # legend.position=c(0.975,0.975),
        legend.position='bottom',
        legend.box='vertical',
        # legend.justification=c(0.975,0.975),
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        # legend.title=element_blank(),
        axis.title.y = element_text(color='black'),
        axis.text.x = element_text(angle=90))
plt
gt = ggplot_gtable(ggplot_build(plt))
which.ylab = grep('ylab-l', gt$layout$name)
which.axes = grep('axis-l', gt$layout$name)
axis.rows  = gt$layout$t[which.axes]
label.col  = gt$layout$l[which.ylab]
gt = gtable::gtable_add_grob(gt, rep(gt$grobs[which.ylab], length(axis.rows)), axis.rows, label.col)
gt = gtable::gtable_filter  (gt, 'ylab-l', invert = TRUE) 
grid::grid.draw(gt)
# plt
png(filename = figname, width = 12.5, height = 7, units = "in", res = 300)
grid::grid.draw(gt)
dev.off()
# ggsave(figname, width=12.5, height=7)


# K-means clustering with N# -----------------------------------------------


nclust <- 10 # Using 10 for greater variability in larger DSM2 dataset
bivars <- c('Flow','Export')
# these have been saved (10 clusters) to the RData file. Just load
# nfex.decdry10 <- clust.plt(nclust, nf.dd.norm.stat, nfex.full.rshp, plt.out.dir, refdate=refdate, bivars=bivars, y.scale=y.scale, only.clust.legend=TRUE) # cluster with nfex.dd.norm.stat but plot with nfex.full.rshp (un-log-normalized full time series for plotting)
# nfex.decdry10.km <- nfex.decdry10$km.cb

# save(list=c('nfex.decdry10', 'nfex.decdry10.km'), file=paste0('data_out/nfex.decdry',nclust,'.RData'))

load(paste0('data_out/nfex.decdry',nclust,'.RData'))

# Determine selected water years for 10 clusters ------------------------------------------

slct.df <- nf.dd.norm
bivars <- c('Flow','Export')
obs <- ncol(slct.df) - 1 # number of days times number of bivariates
clusts <- nfex.decdry10.km$partition

wy.metrics <- data.frame(wy=slct.df$wy, clust=clusts, skill=NA, rmse=NA, mmskill=NA, mmrmse=NA)

slct.df$clust <- clusts
jds <- names(slct.df)[!(names(slct.df) %in% c('wy','clust'))]
for (c in sort(unique(clusts))) {
  
  # get average values for cluster
  mean.clust <- data.frame(jday=jds, mean=NA)
  for (jd in jds) {
    mean.clust$mean[mean.clust$jday==jd] <- mean(slct.df[slct.df$clust==c,names(slct.df) %in% c(jd)])
  } # end jd mean loop
  
  # compute metrics for each wy within cluster
  for (wy in slct.df$wy[slct.df$clust==c]) {
    wy.df <- data.frame(jday=jds,
                        pred=unlist(slct.df[slct.df$wy==wy, jds]))
    wy.df <- merge(wy.df, mean.clust, by='jday')
    wy.metrics$rmse[wy.metrics$wy==wy] <- sqrt(sum((wy.df$pred-wy.df$mean)^2)/nrow(wy.df))
    wy.metrics$skill[wy.metrics$wy==wy] <- 1-(sum((wy.df$pred-wy.df$mean)^2)/
                                                (sum((abs(wy.df$pred-mean(wy.df$mean))+abs(wy.df$mean-mean(wy.df$mean)))^2)))
  } # end wy loop
  
  # compute max and min skill/rmse
  wy.metrics$mmskill[wy.metrics$skill==
                       max(wy.metrics$skill[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Max Skill cluster ', c)
  wy.metrics$mmskill[wy.metrics$skill==
                       min(wy.metrics$skill[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Min Skill cluster ', c)
  wy.metrics$mmrmse[wy.metrics$rmse==
                      max(wy.metrics$rmse[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Max RMSE cluster ', c)
  wy.metrics$mmrmse[wy.metrics$rmse==
                      min(wy.metrics$rmse[wy.metrics$clust==c], na.rm=TRUE)] <- paste0('Min RMSE cluster ', c)
  
} # end cluster loop

wy.metrics <- wy.metrics[order(wy.metrics$clust, wy.metrics$rmse),]
wy.metrics[grepl("Min RMSE",wy.metrics$mmrmse),]
wy.slcts <- wy.metrics[grepl("Min RMSE",wy.metrics$mmrmse),'wy']
wy.slcts <- wy.slcts[!(wy.slcts %in% c(2019))] # chose to only include 2017 and not 2019 in cluster 6 to keep limited to 6 years
sort(wy.slcts) 

# save(wy.metrics, file=paste0('output/',plt.out.dir,'/wymetrics_cluster.RData'))

df.plt.rshp <- nfex.decdry10$df.plt.rshp
plt.cst.rshp <- nfex.decdry10$plt.cst.rshp
plots <- c('Min Skill','Max Skill','Min RMSE','Max RMSE')

wy.to.clust <- slct.df[,c('wy','clust')]
# save(wy.to.clust, file=paste0('output/',plt.out.dir,'/wytoclust.RData'))

# Plot final selection ----------------------------------------------------

slct.wys <- wy.slcts

wycol <- merge(wy.to.clust, col.clust, by='clust')
unslcted <- sort(unique(df.plt.rshp$variable)[!(unique(df.plt.rshp$variable) %in% slct.wys)])
c.colordf <- data.frame(breaks=c(as.numeric(as.character(unslcted)), 
                                 wycol$wy[wycol$wy %in% slct.wys]),
                        values=c(rep('grey',length(unslcted)),
                                 wycol$color[wycol$wy %in% slct.wys]))
flow.brks <- c(5000,10000, 20000, 30000, 40000, 50000,100000)
figname <-  paste0(plt.out.dir,'/Final_WY_Selection_nclusts',nclust,'.png')

# create facet wrap plot for all years with the selected years plotted in respective cluster color
plt <- ggplot() +
  geom_line(data=df.plt.rshp, mapping=aes(x=dt,y=value,color=variable), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
  geom_line(data=df.plt.rshp[df.plt.rshp$variable %in% slct.wys,], mapping=aes(x=dt,y=value,color=variable), linewidth=1) +
  facet_wrap(~bivar, nrow=length(bivars), scales='free') +
  scale_color_manual(breaks= c.colordf$breaks,
                     values=c.colordf$values,
                     # guide=guide_legend(override.aes=list(linetype=0)),
                     name="Cluster") +
  # scale_y_continuous(limits=c(0,500000), breaks=seq(0,500000,100000), expand=c(0,0)) +
  # scale_y_continuous(limits=c(0,25000), expand=c(0,0)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) +
  labs(x='',y='Discharge (cfs)') +
  guides(color=guide_legend(nrow=2)) + 
  facetted_pos_scales(y=list(
    scale_y_continuous(limits=c(0,14000), expand=c(0,0)),
    scale_y_log10(breaks=flow.brks, labels=comma(flow.brks, accuracy=1)))) +
  theme(text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        # legend.key=element_blank(),
        panel.grid.major.y = element_line(size=.25, colour='grey80', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        # legend.position=c(0.975,0.975),
        legend.position='bottom',
        legend.box='vertical',
        # legend.justification=c(0.975,0.975),
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        # legend.title=element_blank(),
        axis.title.y = element_text(color='black'),
        axis.text.x = element_text(angle=90)) + 
  ggtitle(paste0(nclust," Clusters"))
plt
ggsave(figname, width=12.5, height=7)


# Evaluate EC in these windows ---------------------------------------------------------------


# selection period used to determine the cluster 
slct_pers <- data.frame(
  year= sort(wy.slcts) ,
  start=lubridate::ymd(paste0(sort(wy.slcts), "-05-01")),
  end=lubridate::ymd(paste0(sort(wy.slcts), "-11-30"))
)
# simulation period windows determined by evaluating the ##_cluster_wy_slct.html plot
mod_pers <- data.frame(
  year= sort(wy.slcts), # 1991 1992 1994 1995 2000 2002 2005 2013 2015 2017
  start=lubridate::ymd(c('1991-4-1',
                         '1992-3-1',
                         '1993-8-1',
                         '1995-4-1',
                         '2000-2-1',
                         '2002-1-1',
                         '2005-4-1',
                         '2013-1-1',
                         '2015-1-1',
                         '2017-2-1')),
  end=lubridate::ymd(c('1992-5-1',
                        '1993-8-1',
                        '1995-2-1',
                        '1996-7-1',
                        '2001-7-1',
                        '2003-6-15',
                        '2006-6-1',
                        '2015-1-1',
                        '2016-8-1',
                        '2018-7-1'))
)

# create plot-able dataframes
plot_ec <- ec_outputs[,names(ec_outputs) %in% c('Time', ec_eval_locs)]
slct_ec <- plot_ec
mod_ec <- plot_ec
mod_ec$year <- NA
slct_dates <- rep(FALSE, nrow(plot_ec))
mod_dates <- rep(FALSE, nrow(plot_ec))
for (i in 1:nrow(mod_pers)) {
  slct_dates <- slct_dates | (plot_ec$Time >= slct_pers$start[i] & plot_ec$Time <= slct_pers$end[i])
  mod_dates <- mod_dates | (plot_ec$Time >= mod_pers$start[i] & plot_ec$Time <= mod_pers$end[i])
  mod_ec[(plot_ec$Time >= mod_pers$start[i] & plot_ec$Time <= mod_pers$end[i]), 'year'] <- mod_pers$year[i]
}
slct_ec[!slct_dates, -1] <- NA  # Exclude datetime column when replacing
mod_ec[!mod_dates, !(names(mod_ec) %in% c('Time','year'))] <- NA  # Exclude datetime column when replacing
mod_ec <- mod_ec[complete.cases(mod_ec),]
mod_ec$year <- as.factor(mod_ec$year)

plot_ec <- melt(plot_ec, id.vars='Time')
slct_ec <- melt(slct_ec, id.vars='Time')
mod_ec <- melt(mod_ec, id.vars=c('Time','year'))

ds.date.lims <- c(min(plot_ec$Time),
                  max(plot_ec$Time))

# melt_ec <- melt(ec_outputs, id.vars='Time')
# ggplot() +
#   geom_line(data=melt_ec, mapping=aes(x=Time,y=value), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
#   facet_wrap(~variable, ncol=2, scales='free_y')

plt <- ggplot() +
  geom_line(data=plot_ec, mapping=aes(x=Time,y=value), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
  geom_line(data=mod_ec, mapping=aes(x=Time,y=value,color=year), linewidth=1) +
  geom_line(data=slct_ec, mapping=aes(x=Time,y=value,color='selected cluster period'), linewidth=1, color='black', alpha=0.5) +
  facet_wrap(~variable, ncol=1, scales='free_y') +
  scale_x_date(breaks=date_breaks('6 months')) +
  xlab("") +
  ylab("EC (micro-Siemens/cm)") +
  theme(text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        legend.key=element_blank(),
        panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
        panel.grid.major.x = element_line(linewidth=.25, colour='grey80', linetype = 'solid'),
        panel.grid.minor.x = element_line(linewidth=.25, colour='grey90', linetype = 'dashed'),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=90),
        legend.position.inside=c(0.975,0.1),
        legend.justification=c(0.975,0.975),
        legend.spacing=unit(c(0,0,0,0),"null"),
        legend.background = element_rect(fill = "white", color = NULL),
        legend.title=element_blank(),
        axis.title.y = element_text(color='black')) 


plt

pltl <- ggplotly(plt, tooltip=c('value','Time','year'), dynamicTicks=TRUE)
pltl

pltl_name <- paste0(nclust,"_cluster_wy_slct_.html")

saveWidget(pltl, pltl_name, selfcontained=TRUE)
file.rename(pltl_name, paste0("plots/",pltl_name))
