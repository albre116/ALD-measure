library(shiny)

shinyUI(navbarPage("Asymmetric LD (ALD)",
                   
  tabPanel("Upload Data", 
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Upload a Haplotype Frequency File',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        #checkboxInput('header', 'Header', TRUE),      
        tags$a(href="http://www.uvm.edu/~rsingle/software/shiny_ALD/file_formats.txt", target="_blank", "Click here to see the file formats"),

        radioButtons('sep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
        numericInput('tol', label = "tolerance (for sum of haplo.freqs)", value = 0.01, min = 0.01, max = 0.1, step = 0.01),
        tags$p("(sum of haplo.freqs + tolerance) must be >= 1")
      ),
      mainPanel(
        titlePanel("Asymmetric Linkage Disequilibrium (ALD)"),
        helpText("Note: something about tolerance ..."),
        textOutput("text_rawdata1"),
        
        dataTableOutput('raw_data')
      )
    )
  ),

  tabPanel("ALD Plot", 
    sidebarLayout( 
     sidebarPanel(
       tags$h4("ALD plot options"),
       checkboxInput('values', 'Show ALD values', TRUE)
     ),
     mainPanel(
       plotOutput('heatmap')
     )
    )
  ),
  
  tabPanel("ALD Table", 
           htmlOutput("text_ALDdata1"),
           dataTableOutput('plot_data')
  ),
  
  tabPanel("Allele Specific Homozygosity Table", 
    sidebarLayout( 
      sidebarPanel(
        tags$h4("Allele Specific Homozygosity options"),
        uiOutput('choose_locus_pair'), #input$ var: selected_pair
        uiOutput('choose_locus')       #input$ var: selected_locus
      ),
      mainPanel(
        uiOutput("asf_display")
      )
    )
  )
))

  