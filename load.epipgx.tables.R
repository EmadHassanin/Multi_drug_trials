
# local config - not in git-repo
source("config.R")
### load packages
library(readxl)


### Load tables
table.names <- c("ADR", "AED", "INVESTIGATION", "PREGN_DRUGS", "PREGNANCY", "SEIZ_REM",
                 "SUBJECT")
load.table <- function(spec){
  return(read_excel(paste0( dir, tpref,spec ,allsuf), guess_max = 10000) )}
  
# e <- map(table.names, load.table) 
for (d in table.names){
  tdf <- load.table(d)
   assign(d, tdf)
}


#AED$D108_START_date = as.Date(as.numeric(AED$D108_START_date), origin = "1900-01-01")
#AED$D113_END_date = as.Date(as.numeric(AED$D113_END_date), origin = "1900-01-01")
