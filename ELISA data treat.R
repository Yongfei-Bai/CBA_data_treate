#import data use readxl package
library(readxl)
OD450 <-na.omit(
  read_excel('EXP-FV-20220048端点 Abs @ 450 2022-03-18 16-00-44.xlsx',sheet = 2))
OD540 <- na.omit(
  read_excel('EXP-FV-20220048端点 Abs @ 540 2022-03-18 16-01-56.xlsx',sheet = 2))
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
sample <- read.csv(file = 'sample detail.csv')
#merge two data set to one
OD <- merge(sample,a,by = 'Well',all= T)
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
write.csv(file = 'ELISA data.csv',result)
library(ggplot2)

