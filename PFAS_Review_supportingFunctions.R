# PFAS_Review_supportingFunctions.R

collect.frames<-function(file_names,path){
  X<-read.csv(paste0(path,file_names[1]),header = TRUE,row.names = NULL)
  for (i in 2:length(file_names)){
    Y<-read.csv(paste0(path,file_names[i]),header = TRUE,row.names = NULL)
    X<-rbind.fill(X,Y)
  }
  return(X)
}

exp10<-function(x){
  10^x
}
