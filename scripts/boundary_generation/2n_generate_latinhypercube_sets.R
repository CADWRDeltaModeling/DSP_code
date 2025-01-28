
# Libraries to be loaded -------------------------------------------------

library(lhs)

# Set wd to current file --------------------------------------------------

rm(list=ls(all=TRUE)) #start with empty workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Setworking directory to this file's directory

# define convert_to_numeric fxn -------------------------------------------

# Function to convert to a numeric scale
convert_to_numeric <- function(column) {
  # Assign numeric codes (starting from 0 for the first level)
  numeric_column <- seq(1,length(column))
  numeric_column[is.na(column)] <- NA
  
  # Rescale to 0-1 if there are more than one unique values
  if (length(unique(na.omit(column))) > 1) {
    numeric_column <- (numeric_column - min(numeric_column, na.rm = TRUE)) / 
      max(numeric_column, na.rm = TRUE)
  }
  
  numeric_column[numeric_column>.99 | is.na(numeric_column)] <- 99
  
  return(numeric_column)
}

# input data --------------------------------------------------------------

seed <- 916
nLHScol <- 15
nLHSsamples <- 100

version <- 'v4'

years <- c('1991','1992','1994','1995','2000','2002','2005','2013','2015','2017')

mod_pers <- data.frame(
  year= sort(years), # 1991 1992 1994 1995 2000 2002 2005 2013 2015 2017
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

pert_desc <- data.frame(year=years,
                        tide=c('Regular','Shifted + 100d', 'Shifted - 50d','Subtidal Pert',
                               rep(NA, length(years)-4)),
                        dcc=c('Regular','Perturbed','Always Open','Always Closed',
                              rep(NA, length(years)-4)),
                        flows=c('Regular','Perturbed - v1','Perturbed - v2',
                                'Sacramento + 10%','Sacramento - 10%',
                                rep(NA, length(years)-5)),
                        dcd=c('Regular','Perturbed - v1','Perturbed - v2',
                              rep(NA, length(years)-3)),
                        suis=c('Regular','Perturbed - v1', 'Perturbed - v2',
                               rep(NA, length(years)-3)))
pert_df <- as.data.frame(lapply(pert_desc, convert_to_numeric))
pert_df <- rbind(pert_df, rep(99, ncol(pert_df)))

# Check setup
names(pert_df) == names(pert_desc)



# Run with no matching criteria -------------------------------------------------------


set.seed(seed)
num.df <- randomLHS(nLHSsamples,nLHScol)
dup_rows <- TRUE

while (dup_rows) {

  res_df <- data.frame(matrix(NA,
                              nrow <- nLHSsamples,
                              ncol <- ncol(pert_desc)))
  names(res_df) <- names(pert_desc)
  rownames(res_df) <- seq(1,nLHSsamples)# rownames(target.df)
  
  # Determine result from LHC in english (not numbers)
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
  
  # Check for duplicated rows
  any_duplicates <- any(duplicated(res_df))
  
  # Print the result
  if (any_duplicates) {
    print("There are duplicated rows.")
    seed <- seed + 1
    num.df <- randomLHS(nLHSsamples,nLHScol)
    dup_rows <- TRUE
    # Optionally, show which rows are duplicates
    # duplicated_rows <- duplicated(res_df) | duplicated(res_df, fromLast = TRUE)
    # print(res_df[duplicated_rows, ])
  } else {
    print(paste0("No duplicated rows found. Seed = ", seed))
    dup_rows <- FALSE
  }
  

} # Seed ends up being 918

# Calculate percentages for each column
percentage_list <- lapply(res_df, function(column) {
  value_counts <- table(column)                # Count occurrences
  percentages <- (value_counts / nrow(res_df)) * 100  # Calculate percentages
  percentages
})

# Display percentages for each column
# percentage_list

df_out <- merge(res_df, mod_pers, by='year')
df_out$case <- paste0('lhc_',seq(1,nLHSsamples))
df_out$rndays <- as.integer(lubridate::ymd(df_out$end) - lubridate::ymd(df_out$start))

write.csv(df_out, paste0("./data_out/lhc_",version,".csv"), row.names=FALSE)

# # Run (look for matching criteria) -------------------------------------------------------
# match <- FALSE
# set.seed(seed)
# num.df <- randomLHS(nLHSsamples,nLHScol)
# 
# while (!match) {
#   res_df <- data.frame(matrix(NA,
#                               nrow <- nLHSsamples,
#                               ncol <- ncol(pert_desc)))
#   names(res_df) <- names(pert_desc)
#   rownames(res_df) <- seq(1,nLHSsamples)# rownames(target.df)
#   
#   # Determine result from LHC
#   for (c in seq(1,ncol(res_df))) {
#     for (r in seq(1,nLHSsamples)) {
#       cell <- num.df[r,c]
#       pert.col <- pert_df[,c]
#       for (n in seq(1,length(pert.col))) {
#         if (cell>=pert.col[n] & cell<pert.col[n+1]) {
#           res_df[r,c] <- pert_desc[n,c]
#           break
#         } # end if
#       } # end pert check
#     } # end row loop
#   } # end col loop
#   
#   res_df
#   
#   # check for one row of 2008 with all regular
#   for (rdf in seq(1,nLHSsamples)) {
#     if (res_df$year[rdf]==2008 && all(res_df[rdf,2:ncol(pert_desc)]=='Regular')) {
#       match <- TRUE
#     }
#   }
#   
#   if (match) {
#     break
#   } else {
#     seed <- seed + 1
#     set.seed(seed)
#     num.df <- randomLHS(nLHSsamples,nLHScol)
#   }
#   
#   if (seed%%100==0) {
#     print(seed)
#   }
# }
# 
# seed # ends up being 113
# res_df

# write.csv(res_df, paste0("./data_out/lhc_",version,".csv"), row.names=FALSE)

# for adding rows to existing LHS use:
# b <- augmentLHS(num.df,m=2)


