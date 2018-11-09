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

### wibor
setwd(paste0(home,"/dane/wibor"))
wibor_file_names <- c("1m","3m","6m","9m","1y")
wibor_var_names <- c("1m","3m","6m","9m","12m")
for (i in 1:length(wibor_file_names)){
  name <- paste("plopln",wibor_file_names[i],"_m.csv",sep="")
  temp <- read.csv(name)
  temp <- temp[-nrow(temp),]
  #temp <- rev(temp[,5])
  assign(wibor_var_names[i],temp[,5])
}
n_months <- nrow(partia)

data <- matrix(`1m`[(length(`1m`)-n_months+1):length(`1m`)])
data <- cbind(data, matrix(`3m`[(length(`3m`)-n_months+1):length(`3m`)]))
data <- cbind(data, matrix(`6m`[(length(`6m`)-n_months+1):length(`6m`)]))
data <- cbind(data, partia)

yields <- data
maturities <- c(1,3,6,12,24,36,48,60,72,84,96,108,120)
