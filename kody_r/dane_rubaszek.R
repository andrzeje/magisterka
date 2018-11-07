#załadowanie danych Rubaszka

home <- "C:/Users/apalu/Desktop/magisterka"
setwd(paste0(home,"/dane"))
temp <- read.csv2("R8.csv")
yields <- temp[,-1]
maturities <- c(1,3,6,12,24,36,48,60,72,84,96,108,120)
