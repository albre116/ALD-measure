library(shiny)

shinyUI(fluidPage(
  titlePanel("ALD: Asymmetric Linkage Disequilibrium"),
  sidebarLayout(
    sidebarPanel(
      tags$a(href="http://www.uvm.edu/~rsingle/software/shiny_ALD/file_formats.txt", target="_blank", "Click here to see the file formats"),
      fileInput('file1', 'Upload a Haplotype Frequency File',
        accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
     #checkboxInput('header', 'Header', TRUE),      
      radioButtons('sep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
      numericInput('tol', label = "tolerance (for sum of haplo.freqs)", value = 0.01, min = 0.01, max = 0.1, step = 0.01),
      tags$p("(sum of haplo.freqs + tolerance) must be >= 1"),
      tags$hr(),
      tags$h4("ALD Plot options"),
      checkboxInput('values', 'Show ALD values', TRUE),
      tags$hr(),
      tags$h4("Allele Specific Homozygosity options"),
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
                 htmlOutput("text_ALDdata1"),
                 dataTableOutput('plot_data')),
        tabPanel("raw Data", dataTableOutput('raw_data'))
      )
    )
    
    
  )
))