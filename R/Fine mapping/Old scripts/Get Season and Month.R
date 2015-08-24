# Get the season and month function
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES), ".", fixed=TRUE), "[", 2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5]                  <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8]                  <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11]                 <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2]  <- "Winter"
  return(ret)
}

getMonth <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES), ".", fixed=TRUE), "[", 2)))
  return(mmonths)
}


getDate <- function(DATES) {
  Date <- rep(NA, length(DATES))
  for(x in 1:length(DATES)){
  Date[x] <- paste0(strsplit(as.character(DATES[x]), ".", fixed=TRUE)[[1]][c(1,2)], collapse = "/")
  } 
  return(Date)
}

QTLdataAll <- cbind(QTLdataAll, Date = getDate(QTLdataAll[,"Schlupf"]))
