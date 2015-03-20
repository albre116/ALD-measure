#' Compute ALD.
#'
#' A function to compute asymmetric Linkage Disequilibrium measures (ALD) for polymorphic genetic data. These measures are identical to the correlation measure (r) for bi-allelic data.
#' 
#' @param dat	A data.frame with 5 required variables:
#'      \tabular{ll}{
#'           \code{haplo.freq} \tab A numeric vector of haplotype frequencies.\cr
#'            \code{locus1} \tab A character vector indentifying the first locus.\cr
#'            \code{locus2} \tab A character vector indentifying the second locus.\cr
#'            \code{allele1} \tab A character vector indentifying the allele at locus 1.\cr
#'            \code{allele2} \tab A character vector indentifying the allele at locus 1.\cr
#'            }
#' @param tolerance A threshold for the sum of the haplotype frequencies.
#'    If the sum of the haplotype frequencies is greater than 1+tolerance or less
#'    than 1-tolerance an error is returned.
#'
#' @return The return value is a dataframe with the following components:
#'  \tabular{ll}{
#'  \code{locus1}	\tab The name of the first locus.\cr
#'  \code{locus2}	\tab The name of the second locus.\cr
#'  \code{F.1}	\tab Homozygosity (expected under HWP) for locus 1.\cr
#'  \code{F.1.2}	\tab Conditional homozygosity* for locus1 given locus2.\cr
#'  \code{F.2}	\tab Homozygosity (expected under HWP) for locus 1.\cr
#'  \code{F.2.1}	\tab Conditional homozygosity* for locus2 given locus1.\cr
#'  \code{ALD.1.2} \tab Asymmetric LD for locus1 given locus2.\cr
#'  \code{ALD.2.1} \tab  Asymmetric LD for locus2 given locus1.\cr
#'  }
#'  *Overall weighted haplotype-specific homozygosity for the first locus given the second locus.
#'
#' @section Details:
#' A warning message is given if the sum of the haplotype frequencies is greater than 1.01 or less than 0.99. The haplotype frequencies that are passed to the function are normalized within the function to sum to 1.0 by dividing each frequency by the sum of the passed frequencies.
#'
#' @examples
#' data(hla.freqs)
#' hla.a_b <- hla.freqs[hla.freqs$locus1=="A" & hla.freqs$locus2=="B",]
#' compute.ALD(hla.a_b)
#' hla.freqs$locus <- paste(hla.freqs$locus1, hla.freqs$locus2, sep="-")
#' compute.ALD(hla.freqs[hla.freqs$locus=="C-B",])

compute.ALD <- function(dat, tolerance = 0.005) {
  names.req <- c("locus1", "locus2", "allele1", "allele2", "haplo.freq")
  check.names <- names.req %in% names(dat)
  if (sum(!check.names) > 0)
    stop("The following required variables are missing from the dataset: ",
      paste(names.req[!check.names], collapse = " "))
  if (sum(dat$haplo.freq) > 1 + tolerance)
    stop("Sum of haplo freqs > 1; sum=", sum(dat$haplo.freq))
  if (sum(dat$haplo.freq) < 1 - tolerance)
    stop("Sum of haplo freqs < 1; sum=", sum(dat$haplo.freq))
  if (abs(1 - sum(dat$haplo.freq)) > 0.01)
    warning("Sum of haplo freqs is not 1; sum=", sum(dat$haplo.freq))
  locus1 <- unique(dat$locus1)
  locus2 <- unique(dat$locus2)
  
  # normalize haplotype freqs to sum to 1.0
  dat$haplo.freq <- dat$haplo.freq/sum(dat$haplo.freq)
  
  # generate allele freqs based on haplo freqs
  by.vars1 <- list(dat$allele1)
  by.vars2 <- list(dat$allele2)
  names(by.vars1) <- c("allele1")
  names(by.vars2) <- c("allele2")
  af1 <- aggregate(dat$haplo.freq, by = by.vars1, FUN = sum)
  af2 <- aggregate(dat$haplo.freq, by = by.vars2, FUN = sum)
  names(af1)[length(names(af1))] <- "allele.freq1"
  names(af2)[length(names(af2))] <- "allele.freq2"
  mrg1 <- merge(dat, af1, by.x = c("allele1"), by.y = c("allele1"), all.x = T, all.y = F)
  mrg2 <- merge(mrg1, af2, by.x = c("allele2"), by.y = c("allele2"), all.x = T, all.y = F)
  dat <- mrg2
  
  F.1 <- 0
  F.2.1 <- 0
  F.2 <- 0
  F.1.2 <- 0
  for (a1 in unique(dat$allele1)) {
    af.1 <- unique(dat$allele.freq1[dat$allele1 == a1])
    F.2.1 <- sum(F.2.1, sum(dat$haplo.freq[dat$allele1 == a1]^2)/af.1)
    F.1 <- F.1 + af.1^2
  }
  for (a2 in unique(dat$allele2)) {
    af.2 <- unique(dat$allele.freq2[dat$allele2 == a2])
    F.1.2 <- sum(F.1.2, sum(dat$haplo.freq[dat$allele2 == a2]^2)/af.2)
    F.2 <- F.2 + af.2^2
  }
  if (F.2 == 1) {
    F.2.1.prime <- NA
    ALD.2.1 <- NA
  } else {
    F.2.1.prime <- (F.2.1 - F.2)/(1 - F.2)
    ALD.2.1 <- sqrt(F.2.1.prime)
  }
  if (F.1 == 1) {
    F.1.2.prime <- NA
    ALD.1.2 <- NA
  } else {
    F.1.2.prime <- (F.1.2 - F.1)/(1 - F.1)
    ALD.1.2 <- sqrt(F.1.2.prime)
  }
  ALD <- data.frame(locus1, locus2, F.1, F.1.2, F.2, F.2.1, ALD.1.2, ALD.2.1)
  cat("   ALD.x.y is asymmetric LD for locus x conditioned on locus y\n\n")
  return(ALD)
}

