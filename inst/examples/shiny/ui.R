library(shiny)

shinyUI(fluidPage(
  titlePanel("ALD: Asymmetric Linkage Disequilibrium"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose Frequency File',
        accept=c('text/csv', 
          'text/comma-separated-values,text/plain', 
          '.csv')),
      
      tags$p("Need to see file formats?"),
      tags$a(href="shiny.rstudio.com/tutorial", "Click Here!"),
      tags$hr(),
      
     #checkboxInput('header', 'Header', TRUE),      
      radioButtons('sep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
      numericInput('tol', label = "tolerance (for sum of haplo.freqs)", value = 0.01, min = 0.01, max = 0.1, step = 0.01),
      tags$hr(),
      paste("ALD Plot options"),
      checkboxInput('values', 'Show ALD values', TRUE),
      tags$hr(),
      paste("Allele Specific Homozygosity options"),
      uiOutput('choose_locus_pair'), #input$ var: selected_pair
      uiOutput('choose_locus')       #input$ var: selected_locus
    ),

    mainPanel(
      tags$head(
        tags$script(src="d3.js"),
        tags$script(src="lodash.js")
        ),
      tabsetPanel(type = "tabs", 
        tabPanel("ALD Plot", plotOutput('heatmap')), 
        tabPanel("Allele Specific Homozygosity Table", uiOutput("asf_display")),
        tabPanel("ALD Table", 
                 textOutput("text_ALDdata1"),
                 textOutput("text_ALDdata2"),
                 textOutput("text_ALDdata3"),
                 textOutput("text_ALDdata4"),
                 textOutput("text_ALDdata5"),
                 textOutput("text_ALDdata6"),
                 dataTableOutput('plot_data')),
        tabPanel("raw Data", dataTableOutput('raw_data'))
      )
    )
    
    
  )
))