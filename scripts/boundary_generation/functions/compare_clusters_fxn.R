
# reshape.df: function to reshape timeseries dataframe --------------------

# reshape dataframe into julian day columns and wy rows
# TODO: how to handle leap years uncertain, for now just eliminate 366th day and presume day shift is fine.
# which(is.na(nfl.rshp),arr.ind=TRUE)

### Debugging
# ts.df <- northern_flow
# monthstart <- 12
# interval <- 'monthly'

reshape.df <- function(ts.df, monthstart=10, interval='daily') {
  # ts.df is a dataframe with a "Time" column
  colnms <- names(ts.df)
  df.rshp <- ts.df
  if (interval=='daily') {
    df.rshp$jday <- lubridate::yday(df.rshp$Time)
    df.rshp$wyjday <- NA
  } else if (interval=='monthly') {
    df.rshp$mo <- lubridate::month(df.rshp$Time)
    df.rshp$wymoday <- NA
  }
  df.rshp$year <- lubridate::year(df.rshp$Time)
  df.rshp$wy <- NA
  for (year in unique(lubridate::year(df.rshp$Time))) {
    df.rshp[lubridate::year(df.rshp$Time)==(year) & lubridate::month(df.rshp$Time)>=monthstart,'wy'] <- year + 1
    df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart,'wy'] <- year + 1
    if (interval=='daily') {
      df.rshp[lubridate::year(df.rshp$Time)==year & lubridate::month(df.rshp$Time)>=monthstart,'wyjday'] <- 1 + 
        as.numeric(
          df.rshp[lubridate::year(df.rshp$Time)==year & lubridate::month(df.rshp$Time)>=monthstart,'Time'] - 
            lubridate::ymd(paste0(year,'/',monthstart,'/1')))
      df.rshp[lubridate::year(df.rshp$Time)==(year + 1) & lubridate::month(df.rshp$Time)<monthstart,'wyjday'] <- 1 + 
        as.numeric(
          df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart,'Time'] - 
            lubridate::ymd(paste0(year,'/',monthstart,'/1')))
    } else if (interval=='monthly') {
      # need to address leap year
      if (lubridate::leap_year(year+1)) {adj.day <- 1} else {adj.day <- 0}
      df.rshp[lubridate::year(df.rshp$Time)==year & lubridate::month(df.rshp$Time)>=monthstart,'wymoday'] <- 1 + 
        as.numeric(
          lubridate::ymd(df.rshp[lubridate::year(df.rshp$Time)==year & lubridate::month(df.rshp$Time)>=monthstart,'Time']) - 
            lubridate::ymd(paste0(year,'/',monthstart,'/1')))
      df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart & lubridate::month(df.rshp$Time)<3,'wymoday'] <- 1 + 
        as.numeric(
          lubridate::ymd(df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart& lubridate::month(df.rshp$Time)<3,'Time']) - 
            lubridate::ymd(paste0(year,'/',monthstart,'/1'))) # no leap year adjustment for jan/feb
      df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart & lubridate::month(df.rshp$Time)>=3,'wymoday'] <- 1 -
        adj.day + 
        as.numeric(
          lubridate::ymd(df.rshp[lubridate::year(df.rshp$Time)==(year+1) & lubridate::month(df.rshp$Time)<monthstart & lubridate::month(df.rshp$Time)>=3,'Time']) - 
            lubridate::ymd(paste0(year,'/',monthstart,'/1'))) # leap year adjustment for march onward
      
    }
  } # end year loop
  df.rshp <- df.rshp[!(is.na(df.rshp$wy)),]
  
  if (interval=='daily'){
    df.rshp <- melt(df.rshp[,c(colnms[2],'wyjday','wy')], id.vars=c('wyjday','wy'))
    df.rshp <- reshape2::dcast(df.rshp, wy~wyjday, value.var='value') # dataframe now with WYs as rows and julian days as columns, with flow values as the variable
    df.rshp <- df.rshp[,!(names(df.rshp) %in% c(366))] # remove leap year extra day
  } else if (interval=='monthly'){
    df.rshp <- melt(df.rshp[,c(colnms[2],'wymoday','wy')], id.vars=c('wymoday','wy'))
    df.rshp <- reshape2::dcast(df.rshp, wy~wymoday, value.var='value') # dataframe now with WYs as rows and julian days as columns, with flow values as the variable
  }
  df.rshp <- df.rshp[complete.cases(df.rshp),] # remove any years with missing data
  
  df.rshp
}


# Normalize bivariate data ------------------------------------------------

##### Debugging
# df.rshp <- nfex.decdry.rshp
# nbv <- 2

norm.biv.df <- function(df.rshp, nbv) {
  # df.rshp: reshaped df (columns are water year then julian days involved)
  # nbv: number of bivariate variables
  
  df.norm <- data.frame(matrix(nrow=nrow(df.rshp), ncol=ncol(df.rshp)))
  colnames(df.norm) <- colnames(df.rshp)
  df.norm$wy <- df.rshp$wy
  
  for (b in c(1:nbv)) {
    colinds <- c((ncol(df.rshp)-1)/2*(b-1)+2, (ncol(df.rshp)-1)/2*b+1)
    ts <- df.rshp[,colinds[1]:colinds[2]]
    divave.ts <- ts/mean(unlist(ts)) # divide timeseries by average
    log.ts <- log(divave.ts+10^-6) # log transform
    min.bv <- min(log.ts)
    max.bv <- max(log.ts)
    
    norm.ts <- (log.ts-min.bv)/(max.bv-min.bv) # normalize data between 0 and 1
    
    df.norm[,colinds[1]:colinds[2]] <- norm.ts
    
  } # end bivariate loop
  
  df.norm
  
} # end function


# clust.plt: cluster and plot TS clusters ---------------------------------

# create clusters, plot the time series clustered
# nclust is the number of kmeans clusters desired
# df.stat is the dataframe (columns=julian days, rows=wys)
# df.plt is df.stat with the wy column
# plt.out.dir is the directory that plots are saved to
# refdate is the character string which is the 0th julian day
# cols is the plot colors for the number of clusters
# bivars is whether you're performing a bivariate cluster analysis where the julian days are doubled per variable

##### debugging:
# df.stat <-nf.dd.norm.stat
# df.plt <- nfex.full.rshp
# only.clust.legend=TRUE


clust.plt <- function(nclust, df.stat, df.plt, plt.out.dir, 
                      refdate='2000/9/30',
                      cols=c("navy","royalblue1","mediumorchid","orangered3","palegreen3","goldenrod","green4"),
                      bivars=NA,
                      y.scale=NA,
                      only.clust.legend=FALSE) {
  # df.stat has just the julian day and values, df.rshp has the wy column in it
  
  # dir.create(file.path('output',plt.out.dir), showWarnings = FALSE)
  
  # K-means clustering 
  km.cb <- clusterboot(df.stat,B=1000,clustermethod=kmeansCBI,k=nclust,count=FALSE,seed=13541)
  
  # save(km.cb, file='wy_cluster.RData')
  
  df.plt$clust <- km.cb$partition # add cluster number
  wy.to.clust <- df.plt[,c('wy','clust')]
  
  if (all(is.na(bivars))){
    plt.clst <- data.frame(wyday=as.character(names(df.stat)))
    for (c in c(1:nclust)) {
      plt.clst[,paste0('cl.',c)] <- NA
      for (jd in names(df.stat)) {
        plt.clst[plt.clst$wyday==jd,paste0('cl.',c)] <- mean(df.plt[df.plt$clust==c,names(df.plt) %in% jd])
      } # end jd loop
    } # end cluster loop
    plt.cst.rshp <- melt(plt.clst, id.vars='wyday')
  } else {
    # this is a bi-variate cluster timeseries analysis
    colnms <- sapply(strsplit(names(df.stat), "[.]"),`[`,1)
    plt.clst <- data.frame(wyday=unique(colnms))
    for (b in c(1:length(bivars))) {
      bvar <- bivars[b]
      for (c in c(1:nclust)) {
        colname <- paste0('cl.',c, "_",bvar)
        plt.clst[,colname] <- NA
        for (jd in colnms) {
          jdcol <- names(df.plt)[grepl(paste0("^",jd,"."), names(df.plt))][b]
          plt.clst[plt.clst$wyday==jd,colname] <- mean(df.plt[df.plt$clust==c,jdcol])
        } # end jd loop
      } # end cluster loop
    } # end bivariate loop
    plt.cst.rshp <- melt(plt.clst, id.vars='wyday')
    plt.cst.rshp$bivar <- NA
    for (bvar in bivars) {
      plt.cst.rshp$bivar[grepl(paste0("_",bvar), plt.cst.rshp$variable)] <- bvar
    } # end bivariate naming loop
  } # end if bivariate
  # add dates for plotting
  plt.cst.rshp$dt <- lubridate::ymd('1990/1/1')
  for (wyd in unique(plt.cst.rshp$wyday)) {
    plt.cst.rshp$dt[plt.cst.rshp$wyday==wyd] <- lubridate::ymd(refdate) + 
      lubridate::days(plt.cst.rshp$wyday[plt.cst.rshp$wyday==wyd])
  }
  
  
  if (all(is.na(bivars))){
    df.plt.rshp <- melt(df.plt, id.vars=c('wy','clust'))
    df.plt.rshp <- reshape2::dcast(df.plt.rshp, variable~wy, value.var='value') # dataframe now with WYs as columns and julian days as columns, with flow values as the variable
    df.plt.rshp$variable <- as.numeric(as.character(df.plt.rshp$variable))
    df.plt.rshp <- melt(df.plt.rshp, id.vars=c('variable'))
    names(df.plt.rshp) <- c('wyday','variable','value')
  } else {
    # bivariate
    bvlist <- c('x','y')
    df.plt.rshp <- melt(df.plt, id.vars=c('wy','clust'))
    df.plt.rshp$wyday <- sapply(strsplit(as.character(df.plt.rshp$variable), "[.]"),`[`,1)
    # df.plt.rshp$bivar <- NA
    for (b in c(1:length(bivars))) {
      bvar <- bivars[b]
      # df.plt.rshp$bivar[grepl(paste0(".",bvlist[b]), df.plt.rshp$variable)] <- bvar
      tmp.plt.rshp <- reshape2::dcast(df.plt.rshp[grepl(paste0(".",bvlist[b]), df.plt.rshp$variable),], wyday~wy, value.var='value') # dataframe now with WYs as columns and julian days as columns, with flow values as the variable
      tmp.plt.rshp$wyday <- as.numeric(as.character(tmp.plt.rshp$wyday))
      tmp.plt.rshp <- melt(tmp.plt.rshp, id.vars=c('wyday'))
      tmp.plt.rshp$bivar <- bvar
      if (b==1){
        comb.plt.rshp <- tmp.plt.rshp
      } else {
        comb.plt.rshp <- rbind(comb.plt.rshp, tmp.plt.rshp)
      } # end if b==1
    } # end bivariate naming loop
    df.plt.rshp <- comb.plt.rshp
    names(df.plt.rshp) <- c('wyday','variable','value', 'bivar')
  } # end if bivars
  
  # add dates for plotting
  df.plt.rshp$dt <- lubridate::ymd('1990/1/1')
  for (wyd in unique(df.plt.rshp$wyday)) {
    df.plt.rshp$dt[df.plt.rshp$wyday==wyd] <- lubridate::ymd(refdate) + 
      lubridate::days(df.plt.rshp$wyday[df.plt.rshp$wyday==wyd])
  }
  
  # plot each cluster with just the WYs highlighted
  
  max.y.lim <- max(df.plt.rshp$value)*1.2
  
  for (c in 1:nclust) {
    wy.cl.subst <- df.plt.rshp[df.plt.rshp$variable %in% wy.to.clust$wy[wy.to.clust$clust==c],]
    # figname <-  paste0('output/',plt.out.dir,'/Cluster_',c,'.png')
    figname <-  paste0(plt.out.dir,'/Cluster_',c,'.png')
    plot.cluster(wy.cl.subst, df.plt.rshp, plt.cst.rshp, c, figname, bivars=bivars, cols=cols, y.scale=y.scale,,
                 only.clust.legend=only.clust.legend)
  } # end cluster loop
  
  return(list('km.cb'=km.cb, 'df.plt.rshp'=df.plt.rshp, 'plt.cst.rshp'=plt.cst.rshp)) 
}



# plot.cluster: function to plot cluster among all wys --------------------

##### Debugging
savefig=FALSE
plot.mean=FALSE

plot.cluster <- function(wy.cl.subst, df.plt.rshp, plt.cst.rshp, c,
                         figname,
                         savefig=TRUE,
                         plot.mean=TRUE,
                         bivars=NA,
                         plt.title=NA,
                         y.scale=NA,
                         cols=c("navy","royalblue1","mediumorchid","orangered3","palegreen3","goldenrod","green4"),
                         only.clust.legend=FALSE) {
  # wy.cl.subst <- df.plt.rshp[df.plt.rshp$variable %in% wy.to.clust$wy[wy.to.clust$clust==c],]
  wys.inclust <- unlist(as.list(as.vector(unique(wy.cl.subst$variable))))
  wys.notinclust <- levels(df.plt.rshp$variable)[!(levels(df.plt.rshp$variable) %in% wys.inclust)]
  
  # Plotting create initial plot
  plt <- ggplot() +
    geom_line(data=df.plt.rshp, mapping=aes(x=dt,y=value,color=variable), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
    geom_line(data=wy.cl.subst, mapping=aes(x=dt,y=value,color=variable), linewidth=1, alpha=0.3) 
  
  # create color dataframe for scale_color_manual
  if (all(is.na(bivars))) {
    c.colordf <- data.frame(breaks=c(wys.notinclust, 
                                     wys.inclust, 
                                     'Mean'),
                            values=c(rep('grey',length(wys.notinclust)),
                                     rep(cols[c], length(wys.inclust)),
                                     cols[c]))
    # clust_list <- c(paste0('cl.',c))
  } else {
    bv.clust.mean.df <- plt.cst.rshp[grepl(paste0('cl.',c), plt.cst.rshp$variable),]
    bv.clust.mean.df$variable <- 'Mean'
    # clust_list <- c()
    # for (bvar in bivars) {
    #   clust_list <- c(clust_list, paste0('cl.',c,"_",bvar))
    # }
    c.colordf <- data.frame(breaks=c(wys.notinclust, 
                                     wys.inclust, 
                                     rep('Mean', length(bivars))),
                            values=c(rep('grey',length(wys.notinclust)),
                                     rep(cols[c], length(wys.inclust)),
                                     rep(cols[c], length(bivars))))
    # create facet wrap for bivars
    plt <- plt +
      facet_wrap(~bivar, nrow=length(bivars), scales='free')
  
  } # end if bivars
  c.colordf <- c.colordf[order(c.colordf$breaks),]
  if (only.clust.legend) {
    c.colordf <- c.colordf[!(c.colordf$breaks %in% wys.notinclust),]
  }
  
  # Add thick mean line
  if (plot.mean) {
    if (all(is.na(bivars))) {
      plt <- plt +
        geom_line(data=plt.cst.rshp[plt.cst.rshp$variable==paste0('cl.',c),], mapping=aes(x=dt,y=value,color=variable), linewidth=2)
    } else {
      plt <- plt + 
        geom_line(data=bv.clust.mean.df, mapping=aes(x=dt,y=value,color=variable), linewidth=2)
    } # end if bivars
  } # end if plot.mean
  
  plt <- plt +
    scale_color_manual(breaks= c.colordf$breaks,
                       values=c.colordf$values,
                       limits=unique(c.colordf$breaks),
                       # guide=guide_legend(override.aes=list(linetype=0)),
                       name="Cluster") +
    # scale_y_continuous(limits=c(0,500000), breaks=seq(0,500000,100000), expand=c(0,0)) +
    # scale_y_continuous(limits=c(0,25000), expand=c(0,0)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) +
    labs(x='',y='Discharge (cfs)') +
    guides(color=guide_legend(nrow=2)) +
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
          axis.text.x = element_text(angle=90)
    ) # final plot object
  if (!is.na(plt.title)) {
    plt <- plt + ggtitle(plt.title)
  }
  if (all(!is.na(y.scale))) {
    plt <- plt + y.scale
  }
  if (savefig) {
    ggsave(figname, width=12.5, height=7)
  } else {
    return(plt)
  }
}



# compare.clusters: function to create and plot meta clusters -------------

# rshp.dfs is a list of reshaped dataframes (rows=WYs, cols=julian days) to be compared
# kms is a list of k-means cluster results from same order of reshaped dataframes, used for $partition value
# colnames is an array of names of the dataframes entered. must be same length and in order
# if plot=TRUE, then plots will be made of each meta cluster
# plot.dir is the directory where the plots will be saved
# facet is a list matching the colnames where each dataframe goes with respect to ggplot facets
# refdate is a character string (ymd) of the date with which to add the julian days to, default is 2000/9/30 for Oct 1 start date/analysis
# cols is an array of colors

###### Debugging
# rm(list=ls(all=TRUE)) #start with empty workspace
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory
# # rshp.dfs <- list(nfl.dec.dry.rshp, exp.dec.dry.rshp)
# # kms <- list(nfl.dec.dry.km, exp.dec.dry.km)
# load("compare_clusters_fxn_test.Rdata")
# colnames <- c('Flow','Exports')
# plot <- TRUE
# facet <- c('Flow','Exports')
# refdate <- '2000/11/30'
# cols=c("navy","royalblue1","mediumorchid","orangered3","palegreen3","goldenrod","green4","darkorchid4", "salmon")
# # library(scales)
# # show_col(cols)
# plot.dir <- 'output/compare_dec_dry_flow_exports'
# debuggging
plot=FALSE
plot.dir=NA

compare.clusters <- function(rshp.dfs, kms, colnames, 
                             plot=TRUE, plot.dir=NA, facet=NA, refdate='2000/9/30',
                             cols=c("navy","royalblue1","mediumorchid","orangered3","palegreen3","goldenrod","green4","darkorchid4", "salmon")) {
  library(hash)
  library(reshape2)
  library(lubridate)
  library(ggplot2)
  
  part.compare <- data.frame(wy=rshp.dfs[[1]]$wy)
  for (k in c(1:length(kms))) {
    # part.compare[,colnames[k]] <- kms[[k]]$partition
    tmp.df <- data.frame(wy=rshp.dfs[[k]]$wy, clust=kms[[k]]$partition)
    colnames(tmp.df) <- c('wy', colnames[k])
    part.compare <- merge(part.compare,tmp.df, by='wy')
  }
  
  # create meta cluster groupings
  meta.grp <- hash()
  mg <- 1
  meta.grp.md <- data.frame(wy=part.compare$wy, mg=NA, mg.no=NA) # metadata for metaclusters
  for (wy in part.compare$wy) {
    # get wys from each entry df/km
    grps <- hash()
    for (c in c(1:length(colnames))) {
      grp <- part.compare$wy[part.compare[,colnames[c]]==(part.compare[,colnames[c]][part.compare$wy==wy])]
      # assign(paste0(colnames[c],'.grp'), grp)
      grps[colnames[c]] <- grp
      meta.grp.md[meta.grp.md$wy==wy,colnames[c]] <- part.compare[,colnames[c]][part.compare$wy==wy]
      meta.grp.md[meta.grp.md$wy==wy,paste0(colnames[c],'.no')] <- length(grp)
    }
    
    commons <- Reduce(intersect, grps) # common wys accross different groupings
    if (length(commons)>1) {
      if (length(names(meta.grp))!=0) {
        # see if group already exists
        grp.exs <- FALSE
        for (key in names(meta.grp)) { 
          if (setequal(commons,meta.grp[[key]])) {
            grp.exs <- TRUE
            meta.grp.md[meta.grp.md$wy==wy,'mg'] <- key
          }
        }
        if (!grp.exs) {
          meta.grp[mg] <- commons
          meta.grp.md[meta.grp.md$wy==wy,'mg'] <- mg
          mg <- mg+1
        }
      } else {
        # assign first key to meta.grp hash
        meta.grp[mg] <- commons
        meta.grp.md[meta.grp.md$wy==wy,'mg'] <- mg
        mg <- mg+1
      } # end if meta.grp has names yet
    } # end if commons has more than 1
    meta.grp.md[meta.grp.md$wy==wy,'mg.no'] <- length(commons)
  } # end wy loop to create meta clusters
  
  # plot commonality clusters
  if (plot) {
    if (is.na(plot.dir)) {
      stop("There's no plot directory specified. Please add 'plot.dir' variable to funciton")
    } else {
      dir.create(file.path(plot.dir), showWarnings = FALSE)
    }
    # change shape and features of rshp df
    for (r in  c(1:length(rshp.dfs))) {
      rshp.df <- rshp.dfs[[r]]
      rshp.df <- melt(rshp.df, id.vars=c('wy'))
      
      rshp.df$dt <- lubridate::ymd('1990/1/1')
      for (wy in unique(rshp.df$wy)) {
        rshp.df$dt[rshp.df$wy==wy] <- lubridate::ymd(refdate) + 
          lubridate::days(as.character(rshp.df$variable[rshp.df$wy==wy]))
      }
      rshp.df$wy <- as.factor(rshp.df$wy)
      if (!is.na(facet)) {
        rshp.df$facet <- facet[r]
      } # end if facet
      if (r==1) {
        comb.plt <- rshp.df
      } else {
        comb.plt <- rbind(comb.plt, rshp.df)
      } # end comb.plt
    } # end rshp.df reshape loop
    
    # create plot object with all WYs
    plt <- ggplot() +
      geom_line(data=comb.plt, mapping=aes(x=dt,y=value,color=wy), linewidth=0.5, alpha=0.5, show.legend=FALSE) 
    if (facet>0) {
      plt <- plt +
        facet_wrap(~facet, nrow=length(unique(facet)), scales='free')
    } # end if facet
    
    # plot in meta groups
    meta.col.breaks <- c() # for color dataframe
    meta.col.values <- c() # for color dataframe
    if (length(meta.grp) > length(cols)) {
      stop("The number of meta clusters exceeds the number of colors specified, please specify more colors in the 'cols' variable")
    }
    for (key in names(meta.grp)) {
      plt.meta.clst <- comb.plt[comb.plt$wy %in% meta.grp[[key]],]
      plt <- plt +
        geom_line(data=plt.meta.clst, mapping=aes(x=dt,y=value,color=wy), linewidth=0.5)
      
      meta.col.breaks <- c(meta.col.breaks, meta.grp[[key]])
      meta.col.values <- c(meta.col.values, rep(cols[as.numeric(key)], length(meta.grp[[key]])))
    }
    
    uncat.wys <- unique(rshp.dfs[[1]]$wy)[!(unique(rshp.dfs[[1]]$wy) %in% meta.col.breaks)]
    
    colordf <- data.frame(breaks=c(uncat.wys, meta.col.breaks),
                          values=c(rep('grey',length(uncat.wys)),
                                   meta.col.values))
    plt <- plt + 
      scale_color_manual(breaks= colordf$breaks,
                         values=colordf$values,
                         name="Cluster") +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) +
      labs(x='',y='Discharge (cfs)') +
      guides(color=guide_legend(nrow=2)) +
      theme(text=element_text(size=12),
            panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
            panel.background = element_blank(),
            # legend.key=element_blank(),
            panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
            axis.line = element_line(colour = "black"),
            # legend.position=c(0.975,0.975),
            legend.position='bottom',
            legend.box='vertical',
            # legend.justification=c(0.975,0.975),
            legend.spacing=unit(c(0,0,0,0),"null"),
            legend.background = element_rect(fill = "white", color = NULL),
            # legend.title=element_blank(),
            axis.title.y = element_text(color='black'),
            axis.text.x = element_text(angle=90)
      ) # final plot object to then have 
    
    
    for (key in names(meta.grp)) {
      plt.meta.clst <- comb.plt[comb.plt$wy %in% meta.grp[[key]],]
      gplt <- plt +
        geom_line(data=plt.meta.clst, mapping=aes(x=dt,y=value,color=wy), linewidth=2)
      
      ggsave(paste0(plot.dir,'/metacluster_',key,'.png'), plot=gplt, width=12, height=5)
    } # end cluster loop
    
    for (uwy in as.vector(uncat.wys)) {
      plt.uncat.wy <- comb.plt[comb.plt$wy==uwy,]
      gplt <- plt +
        geom_line(data=plt.uncat.wy, mapping=aes(x=dt,y=value,color=wy), linewidth=2) +
        ggtitle(uwy)
      
      ggsave(paste0(plot.dir,'/uncatwy_',uwy,'.png'), plot=gplt, width=12, height=5)
    }
  } # end if plot
  
  # write out dataframe showing how many members belong to each cluster and how many members are in the metacluster
  meta.grp.md
}


# Plot single WY emphasis against cluster ---------------------------------

#### Debugging
# df.plt <- exp.dec.dry.rshp[c(1:29),]
# km.cb <- exp.dec.dry.km
# wy <- 1993

plt.single.wy.clust <- function(df.plt, km.cb, wy, plt.out.dir) {
  jds <- names(df.plt)[!(names(df.plt) %in% c('wy', 'clust'))]
  
  df.plt$clust <- km.cb$partition # add cluster number
  wy.to.clust <- df.plt[,c('wy','clust')]
  plt.wy <- data.frame(wyday=jds, value=NA)
  for (jd in jds) {
    plt.wy$value[plt.wy$wyday==jd] <- df.plt[df.plt$wy==wy,names(df.plt) %in% jd]
  } # end jd loop
  
  plt.wy.rshp <- melt(plt.wy, id.vars='wyday')
  # add dates for plotting
  plt.wy.rshp$dt <- lubridate::ymd('1990/1/1')
  for (var in unique(plt.wy.rshp$variable)) {
    plt.wy.rshp$dt[plt.wy.rshp$variable==var] <- lubridate::ymd(refdate) + 
      lubridate::days(plt.wy.rshp$wyday[plt.wy.rshp$variable==var])
  }
  plt.wy.rshp$color <- as.character(wy)
  
  df.plt.rshp <- melt(df.plt, id.vars=c('wy','clust'))
  df.plt.rshp <- reshape2::dcast(df.plt.rshp, variable~wy, value.var='value') # dataframe now with WYs as columns and julian days as columns, with flow values as the variable
  df.plt.rshp$variable <- as.numeric(as.character(df.plt.rshp$variable))
  df.plt.rshp <- melt(df.plt.rshp, id.vars=c('variable'))
  names(df.plt.rshp) <- c('wyday','variable','value')
  # add dates for plotting
  df.plt.rshp$dt <- lubridate::ymd('1990/1/1')
  for (var in unique(df.plt.rshp$variable)) {
    df.plt.rshp$dt[df.plt.rshp$variable==var] <- lubridate::ymd(refdate) + 
      lubridate::days(df.plt.rshp$wyday[df.plt.rshp$variable==var])
  }
  
  # plot each cluster with just the WYs highlighted
  
  c <- df.plt$clust[df.plt$wy==wy] # get cluster to assign to plot
  
  max.y.lim <- max(df.plt.rshp$value)*1.2
  
  wy.cl.subst <- df.plt.rshp[df.plt.rshp$variable %in% wy.to.clust$wy[wy.to.clust$clust==c],]
  wys.inclust <- unlist(as.list(as.vector(unique(wy.cl.subst$variable))))
  wys.notinclust <- levels(df.plt.rshp$variable)[!(levels(df.plt.rshp$variable) %in% wys.inclust)]
  c.colordf <- data.frame(breaks=c(wys.notinclust, 
                                   wys.inclust),
                          values=c(rep('grey',length(wys.notinclust)),
                                   rep(cols[c], length(wys.inclust))))
  c.colordf <- c.colordf[order(c.colordf$breaks),]
  
  ggplot() +
    geom_line(data=df.plt.rshp, mapping=aes(x=dt,y=value,color=variable), linewidth=0.5, alpha=0.5, show.legend=FALSE) +
    geom_line(data=wy.cl.subst, mapping=aes(x=dt,y=value,color=variable), linewidth=1, alpha=0.3) +
    geom_line(data=plt.wy.rshp, mapping=aes(x=dt,y=value,color=color), linewidth=2) +
    scale_color_manual(breaks= c.colordf$breaks,
                       values=c.colordf$values,
                       # guide=guide_legend(override.aes=list(linetype=0)),
                       name="Cluster") +
    # scale_y_continuous(limits=c(0,500000), breaks=seq(0,500000,100000), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,max.y.lim), expand=c(0,0)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand=c(0,0)) +
    labs(x='',y='Discharge (cfs)') +
    guides(color=guide_legend(nrow=2)) +
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
          axis.text.x = element_text(angle=90)
    ) # final plot object to then have 
  
  ggsave(paste0(plt.out.dir,'/singleWY_',wy,'.png'), width=12, height=5)
}



# Separate bivariate dataframe --------------------------------------------

seperate.bv.df <- function(df.rshp, bivars) {
  library(hash)
  
  # this is a bi-variate cluster timeseries analysis
  colnms <- sapply(strsplit(names(df.rshp), "[.]"),`[`,1)
  colnms <- colnms[!(colnms %in% c('wy'))]
  base.df <- data.frame(wyday=unique(colnms))
  wys <- df.rshp$wy
  bv.list <- hash()
  for (b in c(1:length(bivars))) {
    bvar <- bivars[b]
    bv.df <- base.df
    for (wy in wys) {
      for (jd in colnms) {
        jdcol <- names(df.plt)[grepl(paste0("^",jd,"."), names(df.plt))][b]
        bv.df[bv.df$wyday==jd,as.character(wy)] <- df.rshp[df.rshp$wy==wy,jdcol]
      } # end jd loop
    } # end wy loop
    bv.list[bvar] <- bv.df
  } # end bivariate loop

  bv.list
}