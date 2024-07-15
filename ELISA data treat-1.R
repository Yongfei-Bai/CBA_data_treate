library(readxl)
OD450 <-na.omit(
  read_excel('EXP-FV-20220048端点 Abs @ 450 2022-03-18 16-00-44.xlsx',sheet = 2))
OD540 <- na.omit(
  read_excel('EXP-FV-20220048端点 Abs @ 540 2022-03-18 16-01-56.xlsx',sheet = 2))
plate1 <- read.csv(file = 'sample detail.csv')
ELISA_data_treat <- function(OD450,OD540,plate){
  colnames(OD450)<- OD450[1,]
  OD450 <- OD450[-1,]
  row.names(OD450) <- OD450$孔
  colnames(OD540)<- OD540[1,]
  OD540 <- OD540[-1,]
  row.names(OD540) <- OD540$孔
  OD450 <- as.data.frame(OD450)
  OD540 <- as.data.frame(OD540)
  #numeric data
  OD450$原始吸收 <- as.numeric(OD450$原始吸收)
  OD540$原始吸收 <- as.numeric(OD540$原始吸收)
  a <- as.data.frame(OD450$原始吸收 - OD540$原始吸收)
  a['Well'] <- OD450$孔
  colnames(a) <- c('OD','Well')
  #merge two data set to one
  OD <- na.omit(merge(plate,a,by = 'Well',all= T))
  OD['OD1'] <- (OD$OD)-min(OD$OD)
  n <- grep('S',OD$TCR)
  standard <- OD[n,]
  #standard curve
  standard$Conc <- as.numeric(as.character(standard$Conc))
  linerfit <-lm(OD1~Conc,data = standard) 
  ##
  result <- as.data.frame(OD[,-6])
  for (i in 1:nrow(OD)) {
    result[i,7]=(((OD[i,7]-linerfit$coefficients[1])/linerfit$coefficients[2])*OD[i,5])
  }
  result['IFN-g'] <- result$V7
  result <- result[,-7]
  result [1,8]<-paste('R_squared') 
  result [2,8]<-summary(linerfit)$r.squared 
  return(result)
}
standardplot <- function(OD450,OD540,plate){
  colnames(OD450)<- OD450[1,]
  OD450 <- OD450[-1,]
  row.names(OD450) <- OD450$孔
  colnames(OD540)<- OD540[1,]
  OD540 <- OD540[-1,]
  row.names(OD540) <- OD540$孔
  OD450 <- as.data.frame(OD450)
  OD540 <- as.data.frame(OD540)
  #numeric data
  OD450$原始吸收 <- as.numeric(OD450$原始吸收)
  OD540$原始吸收 <- as.numeric(OD540$原始吸收)
  a <- as.data.frame(OD450$原始吸收 - OD540$原始吸收)
  a['Well'] <- OD450$孔
  colnames(a) <- c('OD','Well')
  #merge two data set to one
  OD <- na.omit(merge(plate,a,by = 'Well',all= T))
  OD['OD1'] <- (OD$OD)-min(OD$OD)
  n <- grep('S',OD$TCR)
  standard <- OD[n,]
  standard$Conc <- as.numeric(as.character(standard$Conc))
  linerfit <-lm(OD1~Conc,data = standard)
  require(ggplot2)
  p <- ggplot(data = standard,aes(x = Conc,y = OD1))+geom_point()+geom_smooth(method = 'lm')+
    labs(title = paste("R2 = ",summary(linerfit)$r.squared))
  return(p)
}
result <- ELISA_data_treat(OD450,OD540,plate =plate1)
p <- standardplot(OD450,OD540 = OD540,plate = plate1)
p
write.csv(file = 'ELISA data.csv',result)