home <- "C:/Users/apalu/Desktop/magisterka"

setwd(paste0(home,"/dane"))
bonds_file_names <- c("2-letnie","3-letnie","4-letnie","5-letnie","6-letnie","8-letnie","9-letnie","10-letnie","12-letnie")
bonds_var_names <- c("24m","36m","48m","60m","72m","96m","108m","120m","132m")
for (i in 1:length(bonds_file_names)){
  name <- paste("Dane historyczne dla dochodow z obligacji Polska ",bonds_file_names[i],".csv",sep="")
  temp <- read.csv(name)
  
  assign(bonds_var_names[i], as.numeric(sub(',', '.', rev(temp[,2]), fixed = TRUE)))
}

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

used_bonds <- c("3m","6m","9m","12m","24m","36m","48m","60m","120m")
maturities <- matrix(c(1,as.numeric(sub('m', '', used_bonds, fixed = TRUE))))


data <- matrix(`1m`)

for (bond in used_bonds){
  temp <- get(bond)
  n.row<-dim(data)[1]
  if (n.row>length(temp)){
    temp <- rev(temp)
    length(temp)<-n.row
    temp <- rev(temp)
  }
  data <- cbind(data,temp)
}

yields <- as.data.frame(data)

write.csv(yields,"yieldspl.csv")
yields <- yields[-(1:30),]
  
