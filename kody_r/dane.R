#czytamy irs z partia.cz

library(data.table)
home <- "C:/Users/apalu/Desktop/magisterka"
setwd(paste0(home,"/dane"))

partia <-read.csv("partiacsv.csv",header = F)

partia <- setDT(partia)
partia <- dcast(partia,V1~V3,value.var = "V4")
setcolorder(partia, c("V1", "1 rok", "2 roky", "3 roky", "4 roky", "5 let", "6 let", "7 let", "8 let", "9 let", "10 let", "12 let", "15 let", "20 let"))

partia$V1 <- as.Date(partia$V1,"%d.%m.%Y")
partia <- partia[order(V1)]
partia<-partia[,-1]    
partia<-partia[,-((ncol(partia)-2):ncol(partia))]
partia <- apply(partia, 2, function(x) as.numeric(sub(',', '.', x, fixed = TRUE)))

yields <- partia
maturities <- c()