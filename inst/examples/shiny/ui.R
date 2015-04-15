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
      checkboxInput('values', 'Show ALD values', TRUE),
      numericInput('tol', label = "tolerance (for sum of haplo.freqs)", value = 0.01, min = 0.01, max = 0.1, step = 0.01)
      
    ),

    mainPanel(
      tags$head(tags$script(src="d3.js")),
      tabsetPanel(type = "tabs", 
        tabPanel("ALD Plot", plotOutput('heatmap')), 
        tabPanel("ASF Table", uiOutput("asf_display")),
        tabPanel("ALD Data", textOutput("text_ALDdata"),
                 dataTableOutput('plot_data')),
        tabPanel("raw Data", dataTableOutput('raw_data'))
      )
    )
    
    
  )
))