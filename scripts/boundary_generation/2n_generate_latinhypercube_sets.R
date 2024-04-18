
# Libraries to be loaded -------------------------------------------------

library(lhs)

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory


pert_df <- data.frame(year=c(0,.25,.5,.75,99),
                      tide=c(0,.5,99,99,99),
                      dcc=c(0,.4,.6,.8,99),
                      flows=c(0,0.3,0.8,99,99),
                      dcd=c(0,1/3,2/3,99,99),
                      suis=c(0,1/3,2/3,99,99))
pert_desc <- data.frame(year=c(2008,2010,2012,2014),
                        tide=c('Regular','Shifted + 100d',NA,NA),
                        dcc=c('Regular','Perturbed','Always Open','Always Closed'),
                        flows=c('Regular','Perturbed','Sacramento + 10%',NA),
                        dcd=c('Regular','Perturbed - v1','Perturbed - v2',NA),
                        suis=c('Regular','Perturbed - v1', 'Perturbed -v2',NA))

# Check setup
names(pert_df) == names(pert_desc)

seed <- 1
nLHScol <- 15
nLHSsamples <- 7

version <- 'v3'


# Run ---------------------------------------------------------------------
match <- FALSE
set.seed(seed)
num.df <- randomLHS(nLHSsamples,nLHScol)

while (!match) {
  res_df <- data.frame(matrix(NA,
                              nrow <- nLHSsamples,
                              ncol <- ncol(pert_desc)))
  names(res_df) <- names(pert_desc)
  rownames(res_df) <- seq(1,nLHSsamples)# rownames(target.df)
  
  # Determine result from LHC
  for (c in seq(1,ncol(res_df))) {
    for (r in seq(1,nLHSsamples)) {
      cell <- num.df[r,c]
      pert.col <- pert_df[,c]
      for (n in seq(1,length(pert.col))) {
        if (cell>=pert.col[n] & cell<pert.col[n+1]) {
          res_df[r,c] <- pert_desc[n,c]
          break
        } # end if
      } # end pert check
    } # end row loop
  } # end col loop
  
  res_df
  
  # check for one row of 2008 with all regular
  for (rdf in seq(1,nLHSsamples)) {
    if (res_df$year[rdf]==2008 && all(res_df[rdf,2:ncol(pert_desc)]=='Regular')) {
      match <- TRUE
    }
  }
  
  if (match) {
    break
  } else {
    seed <- seed + 1
    set.seed(seed)
    num.df <- randomLHS(nLHSsamples,nLHScol)
  }
  
  if (seed%%100==0) {
    print(seed)
  }
}

seed # ends up being 113
res_df

write.csv(res_df, paste0("./data_out/lhc_",version,".csv"), row.names=FALSE)

# for adding rows to existing LHS use:
# b <- augmentLHS(num.df,m=2)


