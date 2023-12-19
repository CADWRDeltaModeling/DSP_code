# function: pulse_events -------------------------------------------------


## Debugging
# key <- 'Northern Flow'
# yr.slct <- 2008
# ts.df <- data.frame(Time=delta_state$Time, Value=delta_state[,key])
# ts.df <- ts.df[lubridate::year(ts.df$Time)==yr.slct,]
# scale <- pulse_params[[key]][1] # through some testing this seems to be the scaling ratio where the first number is the approximate
# nstep <- nrow(ts.df)
# rescale <- scale
# sh.mean <- pulse_params[[key]][2]
# sh.sd <- pulse_params[[key]][3]

# nstep <- 365 #1200
# rescale <- 0.0000000000001
# sh.mean <- 500
# sh.sd <- 1500

pulse_events <- function(ts.df, nstep=365, rescale=0.015, sh.mean=500, sh.sd=1500) {
  
  ## Peak events
  ## Number of events
  ns <- rpois(1, 2*as.integer(nstep/365)) # enough events for about 2 per year
  if (ns==0) {ns <- 1}
  # print(ns)
  
  ## The length of events
  ls <- rpois(ns, 2.5)
  
  ## The positions of the events. Conditional on the number of events and length
  # of series, these are just spread uniformaly
  ps <- sort(as.integer(runif(ns,min=150,max=nstep-50)))
  
  ## Populate the stimulus for each event with a postive number generally in the
  # range 40-60,000 cfs
  # Create the response
  u <- rep(0,nstep)
  for (isn in 1:ns){
    u[ps[isn]:(ls[isn]+ps[isn])] <- rgamma(1,70,rate=1/80000) # populate the positions with a random massive number
  }
  # print(max(u))
  # ts.plot(u,xlim=c(1,nstep))
  
  # Now create a response function to that stimulus that has mean
  # at shape/rate and peak mode (shape-1)/rate days after the stimulus
  f<-dgamma(1:30,shape=4.5,rate=0.5)
  f <- f/sum(f)
  # print(max(u)*rescale) # Rescale to desired range, including duration of stimulus, value of stimulus and the response
  
  f <- f*rescale
  z <- stats::filter(u,f) # shape the events in u to create z
  # ts.plot(z,xlim=c(1,nstep))
  
  # Random shifts
  ns <- rpois(1, 12*as.integer(nstep/365)) # enough for about 12 per year
  ## Position
  ps <- sort(as.integer(runif(ns,min=10,max=nstep-50)))
  rshift <- rep(0,nstep)
  
  ## The magnitudes of the shifts
  ibase <- 1
  for (isn in 1:ns){
    end <- ps[isn]
    rshift[ibase:end] <- rnorm(1, sh.mean, sh.sd) # Magnitude
    ibase <- end
    # print(ibase)
  }
  # ts.plot(rshift)
  # ts.plot(rshift+z)
  
  ts.df$Edit <- ts.df$Value + z + rshift
  ts.df$Edit[is.na(ts.df$Edit)] <- ts.df$Value[is.na(ts.df$Edit)]
  
  return(ts.df)
}


# function: perturb_all ---------------------------------------------------

perturb_all <- function(delta_df, pulse_params) {
  
  ndo_pass <- TRUE
  min.ndo <- 300
  
  while (ndo_pass) {
    
    edit.df <- data.frame(matrix(nrow=nrow(delta_df),ncol=length(keys(pulse_params))))
    names(edit.df) <- keys(pulse_params)
    edit.df$NF_nonSac <- delta_df$NF_nonSac
    edit.df$Time <- delta_df$Time
    # Create edited fields
    for (key in keys(pulse_params)) {
      # print(key)
      ts.df <- data.frame(Time=delta_df$Time, Value=delta_df[,key])
      # ts.df <- ts.df[lubridate::year(ts.df$Time)==yr.slct,]
      scale <- pulse_params[[key]][1] # through some testing this seems to be the scaling ratio where the first number is the approximate increase in flow on the peaks
      
      # create randomly edited flow field
      ts.df.e <- pulse_events(ts.df, nstep=nrow(ts.df), rescale=scale, sh.mean=pulse_params[[key]][2], sh.sd=pulse_params[[key]][3])
      
      # check for min max criteria
      if (!is.na(pulse_params[[key]][4]) | !is.na(pulse_params[[key]][5])) {
        pass <- TRUE
        while (pass) {
          ts.df.e <- pulse_events(ts.df, nstep=nrow(ts.df), rescale=scale, sh.mean=pulse_params[[key]][2], sh.sd=pulse_params[[key]][3])
          
          # check for min/max criteria
          if (is.na(pulse_params[[key]][5])) {
            if (min(ts.df.e$Edit)<pulse_params[[key]][4]) {
              pass <- TRUE } else { pass <- FALSE}
          } else if (min(ts.df.e$Edit)<pulse_params[[key]][4] | max(ts.df.e$Edit)>pulse_params[[key]][5]) {
            pass <- TRUE } else { pass <- FALSE}
        } # end while loop for min/max criteria
      } # end min/max def check
      edit.df[,key] <- ts.df.e$Edit
    } # end key loop 
    
    # Check for net delta outflow requirements
    net_delta_outflow <- edit.df$NF_nonSac + edit.df$Sacramento + edit.df$`SJR Flow` - edit.df$Exports - 
      delta_df$`Consump. Use`
    
    if (min(net_delta_outflow)>min.ndo) {
      ndo_pass <- FALSE
    } #else {
    # print(paste0("Minimum net delta outflow is ", min(net_delta_outflow)))
    # print(paste0("Maximum Exports is ", max(edit.df$Exports)))
    # print(paste0("Maximum Northern Flow is ", max(edit.df$`Northern Flow`)))
    # print(paste0("Maximum SJR Flow is ", max(edit.df$`SJR Flow`)))
    # } # end ndo check
    
    return(edit.df)
  } # end ndo loop
  
  
}
