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
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.r-project.org")
library(shiny)
library(readxl)
library(ggplot2) 
library(reshape2) 
library(dplyr) 
library(writexl)
#ntext <- eventReactive(input$Submit, {     input$n   })
#set up shiny server
function(input, output, session) {
  MFI_file <- eventReactive(input$Submit,
                       {MFI <- input$MFI_file
                       if (is.null(MFI)){return(NULL)}
                       read.csv(MFI$datapath,header = T)})
  gate_file <- eventReactive(input$Submit,
                        {gate <- input$gate_file
                        if (is.null(gate)){return(NULL)}
                        read.csv(gate$datapath,header = T)})
  
  standard_file <- eventReactive(input$Submit,
                        {standard <- input$standard_file
                        if (is.null(standard)){return(NULL)}
                        read.csv(standard$datapath,header = T)})

  output$result = renderTable({     
    df <- MFI_file()     
    return(df)        
    })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #output$result = renderTable({
   # df <- standard_file()
    #return(df)   
    #})

    
    
    
    
    
    
    
    
}
