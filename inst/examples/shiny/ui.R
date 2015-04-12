library(shiny)

shinyUI(fluidPage(
  titlePanel("ALD: Asymmetric Linkage Disequilibrium"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose Frequency File',
        accept=c('text/csv', 
          'text/comma-separated-values,text/plain', 
          '.csv')),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),      
      radioButtons('sep', 'Separator',
        c(Comma=',',
          Semicolon=';',
          Tab='\t'),
        ','),
      checkboxInput('values', 'Show ALD values', TRUE)      
    ),

    mainPanel(
      includeHTML("d3Script.js"),
      tabsetPanel(type = "tabs", 
        tabPanel("ALD Plot", plotOutput('heatmap')), 
        #tabPanel("ASF Table", dataTableOutput("asf_table")),
        tabPanel("Testing!", uiOutput("testingAsf")),
        tabPanel("Data", dataTableOutput('contents'))
      )
    )
    
    
  )
))