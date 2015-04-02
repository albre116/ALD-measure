library(shiny)

shinyServer(function(input, output) {
 
  # display data
   output$contents <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    data <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
    # remove haplotypes with freq == 0 
    rm.inds <- data[, dim(data)[2]] == 0
    data <- data[!rm.inds, ]
    data
  })
  
   
   # plot of asymetric LD
  output$heatmap <- renderPlot({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    data <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
    nloci <- dim(data)[2]-1
    loci <- names(data)[1:nloci]
    
    ald.allpairs <- NULL
    for (i in 1:(nloci-1)){
      for (j in (i+1):nloci){
        bi.data <- get_bilocus_data(data, i, j)
        ald.allpairs <- rbind(ald.allpairs, compute.ALD(bi.data))
      }
    }
    
    # prepare the data for plotting function
    ald.allpairs$pop <- "ALD" 
    pop.use <- "ALD"
    map.order <- c("A","B","DRB1")
    ld.matrix.plot.2vars(dat=ald.allpairs, pop.name=pop.use,
      ld.varnames=c("ALD.1.2","ALD.2.1"), map.order=map.order,
      ld.labnames=c("ALD",""), bw=T, xlab.shift=1, ylab.shift=-0, values=T)
    title(sub=paste(pop.use,": Asymmetric LD\n row gene conditional on
      column gene",sep=""),font.sub=2,cex.sub=1.2)
    
  })
  
  
  # --------- HELPER FUNCTIONS --------------
  
  # Get the frequency data for two loci from the frequency file and format it 
  # for the compute.ALD() function.

  get_bilocus_data <- function(data, i, j){
    data.2 <- data[,c(i, j, nloci+1)]
    data.2$pair <- paste(data.2[, 1], data.2[, 2], sep='-')
    aggregate.freqs <- by(data.2, data.2$pair, function(x) sum(x[, 3]))
    allele_combos <- names(aggregate.freqs)
    haplo.freq <- as.numeric(aggregate.freqs)
    alleles <- t(data.frame(strsplit(allele_combos, split='-')))
    colnames(alleles) <- c('allele1', 'allele2')
    bi.data <- cbind(data.frame(haplo.freq, locus1=rep(loci[i], length(haplo.freq)), 
      locus2=rep(loci[j], length(haplo.freq))), alleles)
    row.names(bi.data) <- NULL
    # remove freq == 0 rows 
    rminds <- bi.data$haplo.freq == 0
    bi.data <- bi.data[!rminds,]
    return(bi.data)
  }
  
  })

