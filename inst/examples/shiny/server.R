library(shiny)
library(asymLD)

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
    data <- data[data[, dim(data)[2]] != 0, ] # remove freq == 0 rows 
    pop.use <- names(data)[dim(data)[2]]
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
    ald.allpairs$pop <- pop.use
    ld.matrix.plot.2vars(dat=ald.allpairs, pop.name=pop.use,
      ld.varnames=c("ALD.1.2","ALD.2.1"), map.order=loci,
      ld.labnames=c("ALD","ALD"), bw=T, xlab.shift=1, ylab.shift=-0, values=T)
    title(sub=paste(pop.use,": Asymmetric LD\n row gene conditional on
      column gene",sep=""),font.sub=2,cex.sub=1.2)
    
  })
  
  
  # -------------------------- HELPER FUNCTIONS -------------------------------
  
  # ---------------------------------------------------------------------------
  # Get the frequency data for two loci from the frequency file and format it 
  # for the compute.ALD() function.
  # ---------------------------------------------------------------------------
  get_bilocus_data <- function(data, i, j){
    nloci <- dim(data)[2]-1
    loci <- names(data)[1:nloci]
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
    return(bi.data)
  }
  
  
  # ---------------------------------------------------------------------------
  # Function for plotting the ALD heatmap
  # ---------------------------------------------------------------------------
  
  ld.matrix.plot.2vars <- function(dat, pop.name, ld.varnames, map.order, 
    ld.labnames=ld.varnames, bw=FALSE, xlab.shift=1, ylab.shift=-5, values=FALSE, cex.ld=.9) {
    dat2 <- dat[dat$pop==pop.name,]  
    
    # swap order of locus1 & locus2 in cases where they are not listed in 
    # map.order with locus1 coming before locus2
    dat2 <- dat2[dat2$locus1 %in% map.order & dat2$locus2 %in% map.order,]
    dat2$locus1 <- factor(dat2$locus1, levels=map.order)
    dat2$locus2 <- factor(dat2$locus2, levels=map.order)
    dat2$swap <- rep(0,dim(dat2)[1]) 
    dat2$swap[as.numeric(dat2$locus1) > as.numeric(dat2$locus2)] <- 1
    dat2$swap.locus1 <- dat2$locus1
    dat2$locus1[dat2$swap==1] <- dat2$locus2[dat2$swap==1]
    dat2$locus2[dat2$swap==1] <- dat2$swap.locus1[dat2$swap==1]  
    
    # make local copies of the vars on the dataframe to use from now on
    pval.color1 <- dat2[,names(dat2)==ld.varnames[1]]
    pval.color2 <- dat2[,names(dat2)==ld.varnames[2]]
    locus1 <- dat2$locus1
    locus2 <- dat2$locus2
    
    n.loci <- length(map.order)
    loci <- 1:length(map.order)
    names(loci) <- map.order
    
    x <- 1:n.loci
    y <- 1:n.loci
    names(x) <- map.order
    names(y) <- map.order
    
    # set up a matrix w/ rownames & colnames = loci.ordered
    z <- outer(x, y, FUN="+") 
    
    for (i in 1:(n.loci-1)) 
    {
      for (j in i:(n.loci))
      {
        if (i == j) 
        {
          z[i,j] <- -1 # a value outside of the zlim range for the diagonal
        } else 
        {
          # ld.varnames[1] on the upper triangle
          if (sum(locus1==names(x)[i] & locus2==names(y)[j]) > 0)  
          {
            #1st entry from pval.color1
            z[i,j] <- pval.color1[locus1==names(x)[i] & locus2==names(y)[j]] 
          } else
          {
            z[i,j] <- -1
          }
          # ld.varnames[2] on the lower triangle
          if (sum(locus1==names(x)[i] & locus2==names(y)[j]) > 0)  
          {
            #2nd entry from pval.color2
            z[j,i] <- pval.color2[locus1==names(x)[i] & locus2==names(y)[j]] 
          } else
          {
            z[j,i] <- -1
          }
        }  
      }
    }
    z[n.loci,n.loci] <- -1 # take care of bottom right element
    
    zz <- NULL
    for (i in (n.loci):1) 
    {
      zz <- cbind(zz, z[i,]) # rotate matrix z by +90 degrees
    }
    
    # NB: can change zz, but not z, since z is the return value
    zz[zz==-1] <- 9 # code untyped locus pairs as 9 for coloring
    
    # plot the LD values
    z.lim <- c(0,1) #zlim: plot only considers values between these 2 values
    stat.breaks <- c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
    if (bw)  stat.colors <- rev(gray((0:(length(stat.breaks)-2))/(length(stat.breaks)-2)))
    if (!bw) stat.colors <- c("gray",rev(rainbow((length(stat.breaks)-1)))[-(length(stat.breaks)-10)])
    # zlim: only use colors in this range
    fields::image.plot(x,y,zz,xaxt="n",yaxt="n",axes=F,xlab="",ylab="",zlim=z.lim,
      horizontal=F,col=stat.colors,breaks=stat.breaks) 
    
    # Add a light grid with dashed lines
    # NB: abline() must be called before mtext() and axis()
    abline(h=0.5:(max(y)+.5),v=0.5:(max(x)+.5),lty=2,col=gray(.3))
    
    
    # add 1st var name as a text string rotated 90 degrees on the right margin
    mtext(side=4,ld.labnames[1],srt=90,font=2,line=ylab.shift) 
    # add 2nd var name as a text string on the bottom margin
    mtext(side=1,ld.labnames[2], font=2,line=xlab.shift) 
    
    # Add axes with locus names  
    axis(3, at=x, labels=names(x), tick=F,las=2)
    axis(2, at=y, labels=rev(names(y)), tick=F,las=2) # rev() since rotated matrix 90 degrees 
    
    
    if (values)
    {
      for (x in 1:n.loci)
      {
        for (y in n.loci:1)
        {
          if (x != (n.loci-y+1))
          {
            yyz <- as.character(round(zz[x,y],2))
            if (!is.na(zz[x,y]))
            {
              if (nchar(as.character(yyz))==3) yyz <- paste(as.character(yyz),"0",sep="")
              yyz <- substr(yyz,2,4)
              if (zz[x,y]>.5 & bw==T)
                # font=2 for bold, NB: cex was 1.2
                text(x,y,yyz,cex=cex.ld,col="white",font=1) 
              else 
                # font=2 for bold, NB: cex was 1.2
                text(x,y,yyz,cex=cex.ld,col="black",font=1) 
            }
          }
        }
      }
    }
    return(z)
  }
  
  
  })

