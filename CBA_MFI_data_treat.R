# load required packages
if(!require(readxl)) install.packages("readxl", repos = "http://cran.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.r-project.org")
if(!require(reshape2)) install.packages("reshape2", repos = "http://cran.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.r-project.org")
if(!require(writexl)) install.packages("writexl", repos = "http://cran.r-project.org")
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(writexl)
# 弹出文件选择对话框
print('请选择导出流式数据文件')

chosen_file_MFI <- file.choose()

# 检查是否选择了文件
if (!is.null(chosen_file_MFI)) {
  print(paste("选择的文件路径是:", chosen_file_MFI))
} else {
  print("没有选择文件")
}
#input MFI data,gate detail
MFI <- read.csv(chosen_file_MFI,header = T)

print('请选择流式门控数据')

chosen_file_gate <- file.choose()
# 检查是否选择了文件
if (!is.null(chosen_file_gate)) {
  print(paste("选择的文件路径是:", chosen_file_gate))
} else {
  print("没有选择文件")
}
print('选择标准品定值文件')

chosen_file_standard <- file.choose()

if (!is.null(chosen_file_standard)) {
  print(paste("选择的文件路径是:", chosen_file_standard))
} else {
  print("没有选择文件")
}
gate <- read.csv(chosen_file_gate,header = T)

names(MFI) <- c('sample','P3','P4','P5','P6','P7','P8','P9','P10',
                'P11','P12','P13','P14','P15','P16')
t_MFI <- as.data.frame(t(MFI))
colnames(t_MFI) <- t_MFI[1,]
t_MFI['gate'] <- colnames(MFI)
t_MFI <- merge(t_MFI,gate,by = 'gate',all = T)
MFI_1<- as.data.frame(t(t_MFI))
MFI_1 <- MFI_1[-1,]
MFI_1[nrow(MFI_1),ncol(MFI_1)] <- 'sample'
colnames(MFI_1) <- MFI_1[nrow(MFI_1),]
MFI_1 <- MFI_1[-nrow(MFI_1),]
#MFI_1[] <- lapply(MFI_1, as.numeric)
#standard building
standard <- read.csv(chosen_file_standard,header = T)
all_standard <- na.omit(merge(standard,MFI_1,by = 'sample',all = T))
row.names(all_standard) <- all_standard$sample
all_standard <- lapply(all_standard[,-1],as.numeric)
all_standard_min <- data.frame(sapply(all_standard, function(x) x - x[length(x)]))
log_all_standard <- log10(all_standard_min)
log_all_standard <- log_all_standard[-nrow(log_all_standard),]
linerfit_IL.12p70 <- lm(log_all_standard$IL.12p70 ~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.8 <- lm(log_all_standard$IL.8~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.4 <- lm(log_all_standard$IL.4~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.5<- lm(log_all_standard$IL.5~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.17A <- lm(log_all_standard$IL.17A~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.22 <- lm(log_all_standard$IL.22~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.17F <- lm(log_all_standard$IL.17F~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.1b <- lm(log_all_standard$IL.1b~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.10 <- lm(log_all_standard$IL.10~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.6 <- lm(log_all_standard$IL.6~log_all_standard$Conc,data = log_all_standard)
linerfit_TNF.a <- lm(log_all_standard$TNF.a ~log_all_standard$Conc,data = log_all_standard)
linerfit_IL.2 <- lm(log_all_standard$IL.2~log_all_standard$Conc,data = log_all_standard)
linerfit_IFN.g <- lm(log_all_standard$IFN.g~log_all_standard$Conc,data = log_all_standard)
linerfit_TNF.b <- lm(log_all_standard$TNF.b~log_all_standard$Conc,data = log_all_standard)
#every cytokine standard plot
require(ggplot2)
plot_IL.12p70 <- ggplot(data = log_all_standard,aes(x =IL.12p70,y = Conc))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.12p70)$r.squared),digits = 2)))
plot_IL.8  <- ggplot(data = log_all_standard,aes(x =IL.8,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.8 )$r.squared),digits = 2)))
plot_IL.4  <- ggplot(data = log_all_standard,aes(x =IL.4 ,y =  Conc))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.4 )$r.squared),digits = 2)))
plot_IL.5  <- ggplot(data = log_all_standard,aes(x =IL.5,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.5 )$r.squared),digits = 2)))
plot_IL.17A  <- ggplot(data = log_all_standard,aes(x =IL.17A,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.17A )$r.squared),digits = 2)))
plot_IL.22  <- ggplot(data = log_all_standard,aes(x =IL.2 ,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.22 )$r.squared),digits = 2)))
plot_IL.17F  <- ggplot(data = log_all_standard,aes(x =IL.17F,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.17F )$r.squared),digits = 2)))
plot_IL.1b  <- ggplot(data = log_all_standard,aes(x =IL.1b,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.1b )$r.squared),digits = 2)))
plot_IL.10  <- ggplot(data = log_all_standard,aes(x =IL.10,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.10 )$r.squared),digits = 2)))
plot_IL.6  <- ggplot(data = log_all_standard,aes(x =IL.6,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.6 )$r.squared),digits = 2)))
plot_TNF.a  <- ggplot(data = log_all_standard,aes(x =TNF.a,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_TNF.a )$r.squared),digits = 2)))
plot_IL.2  <- ggplot(data = log_all_standard,aes(x =IL.2,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IL.2 )$r.squared),digits = 2)))
plot_IFN.g  <- ggplot(data = log_all_standard,aes(x =IFN.g,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_IFN.g )$r.squared),digits = 2)))
plot_TNF.b  <- ggplot(data = log_all_standard,aes(x =TNF.b,y = Conc ))+geom_point()+geom_smooth(method = 'lm')+
  labs(title = paste("R2 = ",format((summary(linerfit_TNF.b )$r.squared),digits = 2)))
library(gridExtra)
library(gridGraphics)
plot_all <-grid.arrange(plot_IL.1b ,plot_IL.2,plot_IL.4,plot_IL.5,plot_IL.6,
                                            plot_IL.8 ,plot_IL.10,plot_IL.12p70 ,plot_IL.17A ,plot_IL.17F, 
                                            plot_IL.22,plot_TNF.a ,plot_TNF.b ,plot_IFN.g ,nrow=4,ncol=4)
ggsave(filename = "C:/Users/meixi/Desktop/standard_curve.png",plot_all,width=14,height=12)
# calculate every cytokine conc by standard

MFI_1.1 <-as.data.frame(lapply(MFI_1[,-ncol(MFI_1)],as.numeric)) 
row.names(MFI_1.1) <- row.names(MFI_1)
n <- grep('S8',row.names(MFI_1.1))
for(i in 1:nrow(MFI_1.1)){
  for(j in 1:ncol(MFI_1.1)){
    MFI_1.1[i,j] <- MFI_1.1[i,j] - MFI_1.1[n,j]
  }}
MFI_1.2 <- MFI_1.1[-n,]
#
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,1]=( 10^((log10(MFI_1.2[a,1])-linerfit_IL.12p70$coefficients[1])/linerfit_IL.12p70$coefficients[2]))
}

for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,2]=( 10^((log10(MFI_1.2[a,2])-linerfit_IL.8$coefficients[1])/linerfit_IL.8$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,3]=( 10^((log10(MFI_1.2[a,3])-linerfit_IL.4$coefficients[1])/linerfit_IL.8$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,4]=( 10^((log10(MFI_1.2[a,4])-linerfit_IL.5$coefficients[1])/linerfit_IL.5$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,5]=( 10^((log10(MFI_1.2[a,5])-linerfit_IL.17A$coefficients[1])/linerfit_IL.17A$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,6]=( 10^((log10(MFI_1.2[a,6])-linerfit_IL.22$coefficients[1])/linerfit_IL.22$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,7]=( 10^((log10(MFI_1.2[a,7])-linerfit_IL.17F$coefficients[1])/linerfit_IL.17F$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,8]=( 10^((log10(MFI_1.2[a,8])-linerfit_IL.1b$coefficients[1])/linerfit_IL.1b$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,9]=( 10^((log10(MFI_1.2[a,9])-linerfit_IL.10$coefficients[1])/linerfit_IL.10$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,10]=( 10^((log10(MFI_1.2[a,10])-linerfit_IL.6$coefficients[1])/linerfit_IL.6$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,11]=( 10^((log10(MFI_1.2[a,11])-linerfit_TNF.a$coefficients[1])/linerfit_TNF.a$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,12]=( 10^((log10(MFI_1.2[a,12])-linerfit_IL.2$coefficients[1])/linerfit_IL.2$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,13]=( 10^((log10(MFI_1.2[a,13])-linerfit_IFN.g$coefficients[1])/linerfit_IFN.g$coefficients[2]))
}
for (a in 1:nrow(MFI_1.2)) {
  MFI_1.2[a,14]=( 10^((log10(MFI_1.2[a,14])-linerfit_TNF.b$coefficients[1])/linerfit_TNF.b$coefficients[2]))
}
#数值型数据转换为小数点后两位
MFI_1.2[] <- lapply(MFI_1.2, function(x) round(x, 2))
MFI_1.2['sample'] <- row.names(MFI_1.2)
#df[] <- lapply(df, function(x) round(x, 2))
write_xlsx(MFI_1.2,'C:/Users/meixi/Desktop/CBA data.xlsx')

