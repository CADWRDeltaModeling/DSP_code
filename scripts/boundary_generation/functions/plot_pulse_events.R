# Create plots (need to debug 2a_perturb_historical_ts.R into pulse_events
# to pulse_events() and use the variables there to create the plot here.

# vars pulled from pulse_events():
  # ts.df
  # rshift
  # z
  # stoch

library(ggplot2)
library(lubridate)
library(scales)

plot_pulses <- function(ts.df,u,z,rshift,stoch) {
  ts.df$Original <- ts.df$Value
  ts.df$Pulse <- u
  ts.df$`Attenuated Pulse` <- z
  ts.df$`Attenuated Pulse`[is.na(ts.df$`Attenuated Pulse`)] <- 0
  ts.df$`Random Shifts` <- rshift
  ts.df$`Stochastic Cycles` <- stoch
  ts.df$`Perturbed Final` <- ts.df$Original + ts.df$`Attenuated Pulse` + ts.df$rshift + ts.df$stoch
  
  plt.df <- ts.df[,names(ts.df) %in% c('Time','Original','Pulse', 'Attenuated Pulse',
                                       'Random Shifts','Stochastic Cycles',
                                       'Perturbed Final')]
  
  plt.df <- melt(plt.df, id.vars='Time')
  
  plt.df$type <- plt.df$variable
  plt.df$type <- as.character(plt.df$type)
  plt.df$type[plt.df$variable %in% c('Original','Perturbed Final')] <- 'Timeseries'
  plt.df$type <- factor(plt.df$type, 
                        levels=c('Pulse', 'Attenuated Pulse','Random Shifts',
                                 'Stochastic Cycles','Timeseries'))
  
  plt <- ggplot(data=plt.df, aes(x=Time, y=value, color=variable, linewidth=variable)) +
    geom_line() +
    facet_wrap(~type, ncol=1, scales='free_y') +
    scale_color_manual(breaks=c('Pulse', 'Attenuated Pulse','Random Shifts','Stochastic Cycles',
                                'Original','Perturbed Final'),
                       values=c('black','black','black','black',
                                'blue','red')) +
    scale_linewidth_manual(breaks=c('Pulse', 'Attenuated Pulse','Random Shifts','Stochastic Cycles',
                                    'Original','Perturbed Final'),
                       values=c(1,1,1,1,
                                1.5,1)) +
    ylab('Flow (cfs)') +
    xlab('')+
    # ggh4x::scale_y_facet(
    #   type == "Timeseries",
    #   trans  = "log10",
    #   breaks = breaks_log(),
    #   labels = label_log()
    # ) +
    theme(text=element_text(size=12),
          panel.border=element_rect(colour = "black", fill=NA, linewidth=0.5),
          panel.background = element_blank(),
          legend.key=element_blank(),
          panel.grid.major.y = element_line(linewidth=.25, colour='grey80', linetype = 'dashed'),
          axis.line = element_line(colour = "black"),
          legend.position=c(0.975,0.975),
          legend.justification=c(0.975,0.975),
          legend.spacing=unit(c(0,0,0,0),"null"),
          legend.background = element_rect(fill = "white", color = NULL),
          legend.title=element_blank(),
          axis.title.y = element_text(color='black')
    )  
  
  return(plt)
}
