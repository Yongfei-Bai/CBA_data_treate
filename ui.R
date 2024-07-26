#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

if(!require(shiny)) install.packages("shiny", repos = "http://cran.r-project.org")
#if(!require(shinydashboard)) install.packages("shinydashboard", repos = "http://cran.r-project.org")

library(shiny)
#library(shinydashboard)

# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("14 Cytokine Detection Calculator (适配智研医康14因子试剂)"),

    sidebarPanel(
        fileInput("MFI_file", label = h4("MFI file input"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        hr(),
        fluidRow(column(4, verbatimTextOutput("value"))),
        fileInput("gate_file", label = h4("Gate detail input"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        hr(),
        fluidRow(column(4, verbatimTextOutput("value"))),
        
        fileInput("standard_file", label = h4("Standard file input"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        hr(),
        fluidRow(column(4, verbatimTextOutput("value"))),
        
        helpText("注：MFI文件为流式软件导出的PE.Mean数据；
        Gate数据为流式软件中门控代表细胞因子类型；
        Standard数据为标准品定量文件。建议除MFI文件外，其余文件请勿删除或改动，否则导致程序无法运行！
                 祝您使用愉快"),
        actionButton("Submit", "Submit"),
        downloadButton("downloadresult", "Download_result"),
        downloadButton("downloadplot", "Download_Plot")),
    
    mainPanel(
        tabsetPanel(id = "inTabset",
                tabPanel("Result",tableOutput("result")),
                tabPanel("Standard_plot",plotOutput("plot"))
   # tabsetPanel(type = "tabs",
    #                tabPanel("Result",tableOutput("result"),
     #               tabPanel("Standard_plot",plotOutput("plot")))
                )
                
 )  )          


