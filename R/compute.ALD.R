



compute.ALD <-
function(dat, tolerance=.005)
{
  names.req <- c("locus1","locus2","allele1","allele2","haplo.freq")
  check.names <- names.req %in% names(dat)
  if ( sum(!check.names)>0 )             stop("The following required variables are missing from the dataset: ", paste(names.req[!check.names], collapse=" ")) 
  if ( sum(dat$haplo.freq)>1+tolerance ) stop("Sum of haplo freqs > 1; sum=", sum(dat$haplo.freq) )
  if ( sum(dat$haplo.freq)<1-tolerance ) stop("Sum of haplo freqs < 1; sum=", sum(dat$haplo.freq) )
  if ( abs(1-sum(dat$haplo.freq))> .01 ) warning("Sum of haplo freqs is not 1; sum=", sum(dat$haplo.freq) )
  locus1 <- unique(dat$locus1)
  locus2 <- unique(dat$locus2)
  
 #normalize haplotype freqs to sum to 1.0
  dat$haplo.freq <- dat$haplo.freq/sum(dat$haplo.freq)

 #generate allele freqs based on haplo freqs
  by.vars1 <- list(dat$allele1) 
  by.vars2 <- list(dat$allele2)
    names(by.vars1) <- c("allele1")
    names(by.vars2) <- c("allele2")
  af1 <- aggregate(dat$haplo.freq, by=by.vars1, FUN=sum)
  af2 <- aggregate(dat$haplo.freq, by=by.vars2, FUN=sum)
    names(af1)[length(names(af1))] <- "allele.freq1"
    names(af2)[length(names(af2))] <- "allele.freq2"
  mrg1 <- merge(dat,  af1, by.x=c("allele1"), by.y=c("allele1"), all.x=T, all.y=F) 
  mrg2 <- merge(mrg1, af2, by.x=c("allele2"), by.y=c("allele2"), all.x=T, all.y=F) 
  dat <- mrg2

  F.1 <- 0; F.2.1 <- 0
  F.2 <- 0; F.1.2 <- 0
  for (a1 in unique(dat$allele1))
  {
    af.1 <- unique(dat$allele.freq1[dat$allele1==a1])
    F.2.1 <- sum(F.2.1, sum( dat$haplo.freq[dat$allele1==a1]^2 ) / af.1)    
    F.1 <- F.1 + af.1^2
  }
  for (a2 in unique(dat$allele2))
  {
    af.2 <- unique(dat$allele.freq2[dat$allele2==a2])
    F.1.2 <- sum(F.1.2, sum( dat$haplo.freq[dat$allele2==a2]^2 ) / af.2)
    F.2 <- F.2 + af.2^2
  } 
  if (F.2==1) { F.2.1.prime <- NA; ALD.2.1 <- NA }
  else 
  {
    F.2.1.prime <- (F.2.1 - F.2) / (1 - F.2)
    ALD.2.1 <- sqrt(F.2.1.prime)
  }
  if (F.1==1) { F.1.2.prime <- NA; ALD.1.2 <- NA }
  else
  {
    F.1.2.prime <- (F.1.2 - F.1) / (1 - F.1)
    ALD.1.2 <- sqrt(F.1.2.prime)
  }
 ALD <- data.frame(locus1, locus2, F.1, F.1.2, F.2, F.2.1, ALD.1.2, ALD.2.1) 
 cat("   ALD.x.y is asymmetric LD for locus x conditioned on locus y\n\n")
 return(ALD)
}
