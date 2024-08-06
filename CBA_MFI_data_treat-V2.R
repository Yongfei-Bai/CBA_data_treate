# load required packages
if(!require(readxl)) install.packages("readxl", repos = "http://cran.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.r-project.org")
if(!require(reshape2)) install.packages("reshape2", repos = "http://cran.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.r-project.org")
if(!require(writexl)) install.packages("writexl", repos = "http://cran.r-project.org")
if(!require(tcltk)) install.packages("tcltk", repos = "http://cran.r-project.org")
library(readxl)
library(ggplot2)
library(reshape2)
library(dplyr)
library(writexl)
library(tcltk)
#设定检测项目
set_value <- 14
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
#检查导入文件列数是否符合
if ((ncol(MFI)-1)!= set_value) {
  tcltk::tk_messageBox(
    message = paste("行数不符合设定值！实际列数为", (ncol(MFI)-1)),
    type = "ok"
  )
}

print('请选择流式门控数据')

chosen_file_gate <- file.choose()
# 检查是否选择了文件
if (!is.null(chosen_file_gate)) {
  print(paste("选择的文件路径是:", chosen_file_gate))
} else {
  print("没有选择文件")
}
gate <- read.csv(chosen_file_gate,header = T)
#检查导入文件列数是否符合
if ((nrow(gate))!= set_value) {
  tcltk::tk_messageBox(
    message = paste("行数不符合设定值！实际行数为", (nrow(gate))),
    type = "ok"
  )
}
print('选择标准品定值文件')

chosen_file_standard <- file.choose()

if (!is.null(chosen_file_standard)) {
  print(paste("选择的文件路径是:", chosen_file_standard))
} else {
  print("没有选择文件")
}
standard <- read.csv(chosen_file_standard,header = T)

#names(MFI) <- c('sample','P3','P4','P5','P6','P7','P8','P9','P10',
 #               'P11','P12','P13','P14','P15','P16')

t_MFI <- as.data.frame(t(MFI))
colnames(t_MFI) <- t_MFI[1,]
t_MFI['gate'] <- colnames(MFI)
t_MFI <- merge(t_MFI,gate,by = 'gate',all = T)
MFI_1<- as.data.frame(t(t_MFI))
MFI_1 <- MFI_1[-1,]
MFI_1[nrow(MFI_1),ncol(MFI_1)] <- 'sample'
colnames(MFI_1) <- MFI_1[nrow(MFI_1),]
MFI_1 <- MFI_1[-nrow(MFI_1),]
#standard building
standard <- read.csv(chosen_file_standard,header = T)
all_standard <- na.omit(merge(standard,MFI_1,by = 'sample',all = T))
row.names(all_standard) <- all_standard$sample
all_standard <- lapply(all_standard[,-1],as.numeric)
all_standard_min <- data.frame(sapply(all_standard, function(x) x - x[length(x)]))
log_all_standard <- log10(all_standard_min)
log_all_standard <- as.data.frame(log_all_standard[-nrow(log_all_standard),])
for (i in 1:set_value) {
  var_name <- paste0("linearfit_", i)
  assign(var_name,lm(log_all_standard[,i+1] ~log_all_standard[,1],data = log_all_standard) )
}
require(ggplot2)

#for (i in 1:set_value) {
 # linerfit <- paste0("linerfit_", i)
#  assign(linerfit,lm(log_all_standard[,i+1] ~log_all_standard[,1],data = log_all_standard) )
#   for (j in 1:set_value) { 
#      var_name <- paste0("plot_", j)
#      assign(var_name,ggplot(data = log_all_standard,aes(x =log_all_standard[,j+1],y = log_all_standard[,1]))+geom_point()+geom_smooth(method = 'lm')+
#           labs(title = paste("R2 = ",format((summary(linerfit)[8]),digits = 2)),y=paste('Log_Conc'),
#                x=paste(names(log_all_standard[j+1]))))
#           }
# 
#  }

library(ggplot2)

# 定义函数来生成单个图
generate_plot <- function(i, log_all_standard) {
  var_name <- paste0("plot_", i)  
  x_data <- log_all_standard[[i + 1]]
  y_data <- log_all_standard[[1]]
  current_linerfit <- get(paste0("linearfit_", i))
  
  p <- ggplot(data = NULL, aes(x = x_data, y = y_data)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    labs(title = paste("R2 = ", format((summary(current_linerfit)[8]), digits = 2)),
         y = paste('Log_Conc'),  
         x = paste(names(log_all_standard[i + 1])))
  return(p)
}

# 生成图的列表
plots_list <- list() 
for (i in 1:set_value) {
  plots_list[[paste0("plot_", i)]] <- generate_plot(i, log_all_standard)
}

library(cowplot)
plots_all <- plot_grid(plotlist = plots_list,nrow=3,ncol=5)
plots_all


#计算每组浓度
MFI_1.1 <-as.data.frame(lapply(MFI_1[,-ncol(MFI_1)],as.numeric)) 
row.names(MFI_1.1) <- row.names(MFI_1)
min_standard <- c('S8')
n <- grep(min_standard,row.names(MFI_1.1))
for(i in 1:nrow(MFI_1.1)){
  for(j in 1:ncol(MFI_1.1)){
    MFI_1.1[i,j] <- MFI_1.1[i,j] - MFI_1.1[n,j]
  }}
MFI_1.2 <- MFI_1.1[-n,]
for (a in 1:nrow(MFI_1.2)) {
  for (i in 1:set_value) {  
  current_linerfit <- get(paste0("linearfit_", i))
  MFI_1.2[a,i]=( 10^((log10(MFI_1.2[a,i])-current_linerfit$coefficients[1])/current_linerfit$coefficients[2]))
  }
}  
#数值型数据转换为小数点后两位
MFI_1.2[] <- lapply(MFI_1.2, function(x) round(x, 2))
#MFI_1.2['sample'] <- row.names(MFI_1.2)
#df[] <- lapply(df, function(x) round(x, 2))

csv_file  <- file.choose()
if (!substring(csv_file, nchar(csv_file)-3, nchar(csv_file)) %in% '.csv') {
  csv_file <- paste0(csv_file, '.csv')
}

write.csv(MFI_1.2,csv_file,fileEncoding = "GB18030")















