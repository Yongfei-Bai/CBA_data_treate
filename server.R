#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
if(!require(shiny)) install.packages("shiny", repos = "http://cran.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.r-project.org")
if(!require(reshape2)) install.packages("reshape2", repos = "http://cran.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.r-project.org")
if(!require(writexl)) install.packages("writexl", repos = "http://cran.r-project.org")
if(!require(cowplot)) install.packages("cowplot", repos = "http://cran.r-project.org")
library(shiny)
library(readxl)
library(ggplot2) 
library(reshape2) 
library(dplyr) 
library(writexl)
library(cowplot)
library(gridGraphics)
#ntext <- eventReactive(input$Submit, {     input$n   })
#set up shiny server
function(input, output, session) {
  observeEvent(input$Submit, {
                MFI <- input$MFI_file
                if (is.null(MFI)){return(NULL)}
                MFI_input <- read.csv(MFI$datapath,header = T)
                gate <- input$gate_file
                if (is.null(gate)){return(NULL)}
                gate_input <- read.csv(gate$datapath,header = T)
                standard <- input$standard_file
               if (is.null(standard)){return(NULL)}
               standard_input <- read.csv(standard$datapath,header = T)
               t_MFI <- as.data.frame(t(MFI_input))
               colnames(t_MFI) <- t_MFI[1,]
               t_MFI['gate'] <- colnames(MFI_input)
               t_MFI <- merge(t_MFI,gate_input,by = 'gate',all = T)
               MFI_1<- as.data.frame(t(t_MFI))
               MFI_1 <- MFI_1[-1,]
               MFI_1[nrow(MFI_1),ncol(MFI_1)] <- 'sample'
               colnames(MFI_1) <- MFI_1[nrow(MFI_1),]
               MFI_1 <- MFI_1[-nrow(MFI_1),]
               all_standard <- na.omit(merge(standard_input,MFI_1,by = 'sample',all = T))
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
               require(gridExtra)
               require(gridGraphics)
               plot_all <-plot_grid(plot_IL.1b ,plot_IL.2,plot_IL.4,plot_IL.5,plot_IL.6,
                                       plot_IL.8 ,plot_IL.10,plot_IL.12p70 ,plot_IL.17A ,plot_IL.17F, 
                                       plot_IL.22,plot_TNF.a ,plot_TNF.b ,plot_IFN.g ,nrow=3,ncol=5)
               MFI_1.1 <-as.data.frame(lapply(MFI_1[,-ncol(MFI_1)],as.numeric)) 
               MFI_1.1['sample'] <- MFI_1[,ncol(MFI_1)]
               n <- grep('S8',MFI_1.1[,ncol(MFI_1.1)])
               for(i in 1:nrow(MFI_1.1)){
                 for(j in 1:(ncol(MFI_1.1)-1)){
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
               MFI_1.2[nrow(MFI_1.2),(ncol(MFI_1.2)-1)] <- lapply(MFI_1.2[nrow(MFI_1.2),(ncol(MFI_1.2)-1)] , function(x) round(x, 2))
               sample <- MFI_1.2[, ncol(MFI_1.2)]
               MFI_1.2 <- cbind(sample, MFI_1.2[, -ncol(MFI_1.2)])
               output$result = renderTable(MFI_1.2)
               output$plot = renderPlot(plot_all)
               output$downloadresult <- downloadHandler(
                 filename = function() {
                   paste('null', ".csv", sep = "")
                 },
                 content = function(file) {
                   write.csv(MFI_1.2, file, row.names = F,fileEncoding = "GB18030")})
               output$downloadplot <- downloadHandler(
                 filename = function() {
                   paste('null', ".jpg", sep = "")
                 },
                 content = function(file) {
                  ggsave(file,plot_all,width=16,height=12)})
                  # write.csv(MFI_1.2, file, row.names = F)})
               })
}