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

# pulse_events(ts.df, nstep=nrow(ts.df), rescale=scale, sh.mean, sh.sd, num.pulse, stochit, stochscale)
# nstep <- nrow(ts.df)
# rescale <- scale
# sh.mean <- pulse_params[[key]][2]
# sh.sd <- pulse_params[[key]][3]

pulse_events <- function(ts.df, nstep=365, rescale=0.015, sh.mean=500, sh.sd=1500, num.pulse=2, stochit=FALSE, stochscale=1) {
  
  ## Peak events
  ## Number of events
  ns <- rpois(1, num.pulse*as.integer(nstep/365)) # enough events for about num.pulse per year
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
  # par(mfrow=c(5,1))
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
  
  if (stochit) {
    stoch <- stochify(nstep, scale=stochscale)
    # ts.plot(stoch)
    
    # ts.plot(rshift+z+stoch)
    ts.df$Edit <- ts.df$Value + z + rshift
  } else {
    ts.df$Edit <- ts.df$Value + z + rshift
  }
  
  ts.df$Edit[is.na(ts.df$Edit)] <- ts.df$Value[is.na(ts.df$Edit)]
  
  return(ts.df)
}


# function: stochify ------------------------------------------------------

stochify <- function(nstep, scale=1, nstep_buf=3500, base_period_days=128) {
  ########## This is a perturbation that includes drift plus periodicities
  # of 14,28 and 56 days. The first couple lines are where those periods are
  # basically you set some frequencies and it randomly modulates those as well as creating a trend.
  
  
  lambdabase = 2*pi/base_period_days
  # You could add or remove periods here by changing c(1,2,4)
  # The lowest frequency corresponds to the 1 element (128 days divided by 1)
  lambda = c(1,2,4,8)*2*pi/base_period_days
  
  order_trend <- 2
  order_cycle <- 4
  freqs <- lambda
  rho <- 0.994
  sigma2_eps <- 1e-4
  sigma2_zeta <- 2e-4
  sigma2_kappa <- rep(9e-11,length(lambda))
  sigma2_diffuse <- 1.
  
  mod <- sc_model_build(order_trend,order_cycle,freqs,rho,
                        sigma2_eps=sigma2_eps,sigma2_zeta=sigma2_zeta,
                        sigma2_kappa=sigma2_kappa,sigma2_diffuse=sigma2_diffuse)
  mod$GG[1,2] <- 0.99
  mod$GG[2,2] <- 0.99
  mod$FF[2] <- 0.99
  nstate <- dim(mod$GG)[1]
  sigma_eps <- sqrt(sigma2_eps)
  sigma <- sqrt(diag(mod$W))
  y = rep(NA,nstep_buf)
  m = matrix(nrow=nstate,ncol=nstep_buf)
  m[,]<-0.
  
  
  for (i in 2:nstep_buf){
    m[,i] = mod$GG%*%m[,i-1] + rnorm(nstate,sd=sigma)
    y[i] = mod$FF%*%m[,i] + rnorm(1,sigma_eps)
  }
  ndx <- c(1,order_trend+1+order_cycle*2*(0:(length(freqs)-1)))
  m <- m[ndx,1001:(1000+nstep)]*4
  y <- colSums(m)+rnorm(nstep)*2
  
  y <- y * scale
  
  return(y)
  # par(mfrow=c(4,1))
  # for(jplot in 2:4){ts.plot(m[jplot,])}
  # ts.plot(y)
}

# function: perturb_all ---------------------------------------------------

# debug
# key <- 'Sacramento'

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
      scale <- pulse_params[[key]][1]
      sh.mean <- pulse_params[[key]][2]
      sh.sd <- pulse_params[[key]][3]
      min.crit <- pulse_params[[key]][4]
      max.crit <- pulse_params[[key]][5]
      num.pulse <- pulse_params[[key]][6]
      stochit <- pulse_params[[key]][7]
      stochscale <- pulse_params[[key]][8]
      
      # create randomly edited flow field
      ts.df.e <- pulse_events(ts.df, nstep=nrow(ts.df), rescale=scale, sh.mean, sh.sd, num.pulse, stochit, stochscale)
      
      # debug check
      # ggplot(data=ts.df.e,aes(x=Time)) +
      #   geom_line(aes(y=Value, color='orig')) +
      #   geom_line(aes(y=Edit, color='edit'))
      
      # check for min max criteria
      if (!is.na(min.crit) | !is.na(max.crit)) {
        pass <- TRUE
        while (pass) {
          ts.df.e <- pulse_events(ts.df, nstep=nrow(ts.df), rescale=scale, sh.mean, sh.sd, num.pulse)
          
          # check for min/max criteria
          if (is.na(max.crit)) {
            if (min(ts.df.e$Edit)<min.crit) {
              pass <- TRUE } else { pass <- FALSE}
          } else if (min(ts.df.e$Edit)<min.crit | max(ts.df.e$Edit)>max.crit) {
            pass <- TRUE } else { pass <- FALSE}
        } # end while loop for min/max criteria
      } # end min/max def check
      edit.df[,key] <- ts.df.e$Edit
    } # end key loop 
    
    # Check for net delta outflow requirements
    net_delta_outflow <- edit.df$NF_nonSac + edit.df$Sacramento + edit.df$`SJR Flow` - edit.df$Exports - 
      delta_df$`Consump Use`
    
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
