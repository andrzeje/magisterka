getwd()
setwd("dane/wibor")
wibor_file_names <- c("1m","3m","6m","9m","1y")

for (n_mon in wibor_file_names){
  name <- paste("plopln",n_mon,"_m.csv",sep="")
  assign(name, read.csv(name))
}
