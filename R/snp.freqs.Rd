\name{snp.freqs}
\alias{snp.freqs}
\docType{data}
\title{
Pairwise haplotype frequencies for 4 SNP loci
}
\description{
HLA haplotype frequencies for 6 pairs of SNP loci from the de Bakker et al. 2006 data for 90 unrelated individuals with European ancestry (CEU) from the Centre d'Etude du Polymorphisme Humain (CEPH)
collection obtained from the Tagger/MHC webpage.
}
\usage{data("snp.freqs")}
\format{
  A data frame with 20 observations on the following 5 variables.
  \describe{
    \item{\code{locus1}}{a character vector}
    \item{\code{locus2}}{a character vector}
    \item{\code{allele1}}{a character vector}
    \item{\code{allele2}}{a character vector}
    \item{\code{haplo.freq}}{a numeric vector}
  }
}
\details{
For bi-allelic SNP data the ALD measures are symmetric and equivalent to the r correlation measure of LD.
}
\source{
de Bakker, P. I., G. McVean, P. C. Sabeti, M. M. Miretti, T. Green et al., 2006 A high-resolution HLA and SNP haplotype map for disease association studies in the extended human MHC. Nat.Genet. 38: 1166-1172.
}
\references{
http://www.broadinstitute.org/mpg/tagger/mhc.html
}
\examples{
 data(snp.freqs)
 snps <- c("rs1548306", "rs6923504", "rs4434496", "rs7766854")
 compute.ALD(snp.freqs[snp.freqs$locus1==snps[2] & snp.freqs$locus2==snps[3],])

 snp.freqs$locus <- paste(snp.freqs$locus1, snp.freqs$locus2, sep="-")
 by(snp.freqs,list(snp.freqs$locus),compute.ALD)

#SNP1 & SNP2 : the r correlation & ALD measures are equivalent due to symmetry for bi-allelic SNPs
 p.AB <- snp.freqs$haplo.freq[1]
 p.Ab <- snp.freqs$haplo.freq[2]
 p.aB <- snp.freqs$haplo.freq[3]
 p.ab <- snp.freqs$haplo.freq[4]
 p.A <- p.AB + p.Ab
 p.B <- p.AB + p.aB 
 r.squared <- (p.AB - p.A*p.B)^2 / (p.A*(1-p.A)*p.B*(1-p.B))
 sqrt(r.squared) #the r correlation measure
 compute.ALD(snp.freqs[snp.freqs$locus1==snps[1] & snp.freqs$locus2==snps[2],])

}
\keyword{datasets}
