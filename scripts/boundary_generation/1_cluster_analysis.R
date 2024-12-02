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
library(scales)

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
source('./functions/compare_clusters_fxn.R')
source("./functions/whipple_adapt_clustergram.R")
source("./functions/load_delta_vars.R")

# read in data ------------------------------------------------------------

# DSM2 Inputs/Outputs for ANN
dsm2_filename <- '../../data/dsm2_ann_inputs_base.xlsx'
# ec outputs are reported in micro-Siemens/cm. 1000 mS/cm = 1 mmhos/cm

northern_flow <- load_delta_vars(dsm2_filename, 'northern_flow')
exports <- load_delta_vars(dsm2_filename, 'exports')

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

remove(var, error=FALSE)
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

wycol <- merge(wy.to.clust, col.clust, by='clust')
unslcted <- sort(unique(df.plt.rshp$variable)[!(unique(df.plt.rshp$variable) %in% slct.wys)])
c.colordf <- data.frame(breaks=c(unslcted, 
                                 wycol$wy[wycol$wy %in% slct.wys]),
                        values=c(rep('grey',length(unslcted)),
                                 wycol$color[wycol$wy %in% slct.wys]))
flow.brks <- c(5000,10000, 20000, 30000, 40000, 50000,100000)
figname <-  paste0(plt.out.dir,'/Final_WY_Selection.png')

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
        axis.text.x = element_text(angle=90))
plt
ggsave(figname, width=12.5, height=7)


# K-means clustering with N# -----------------------------------------------


nclust <- 10 # Using 10 for greater variability in larger DSM2 dataset
bivars <- c('Flow','Export')
nfex.decdry10 <- clust.plt(nclust, nf.dd.norm.stat, nfex.full.rshp, plt.out.dir, refdate=refdate, bivars=bivars, y.scale=y.scale, only.clust.legend=TRUE) # cluster with nfex.dd.norm.stat but plot with nfex.full.rshp (un-log-normalized full time series for plotting)
nfex.decdry10.km <- nfex.decdry10$km.cb

# save(nfex.decdry10.km, file=paste0('data_out/nfex.decdry10.kmnfex.decdry10.km.RData'))

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

df.plt.rshp <- nfex.decdry4$df.plt.rshp
plt.cst.rshp <- nfex.decdry4$plt.cst.rshp
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
