\name{hla.freqs}
\alias{hla.freqs}
\docType{data}
\title{
Pairwise haplotype frequencies for 7 HLA loci
}
\description{
HLA haplotype frequencies for 21 pairs of HLA loci in a set of 300 controls from a study of myopericarditis incidence following smallpox vaccination. 
}
\usage{data("hla.freqs")}
\format{
  A data frame with 3063 observations on the following 5 variables.
  \describe{
    \item{\code{haplo.freq}}{a numeric vector}
    \item{\code{locus1}}{a character vector}
    \item{\code{locus2}}{a character vector}
    \item{\code{allele1}}{a character vector}
    \item{\code{allele2}}{a character vector}
  }
}
\source{
Wilson, C., 2010 Identifying polymorphisms associated with risk for the development of myopericarditis following smallpox vaccine. The Immunology Database and Analysis Portal (ImmPort), Study #26.
}
\references{
https://immport.niaid.nih.gov/immportWeb/clinical/study/displayStudyDetails.do?itemList=SDY26
}
\examples{
data(hla.freqs)
hla.a_b <- hla.freqs[hla.freqs$locus1=="A" & hla.freqs$locus2=="B",]
compute.ALD(hla.a_b)
hla.freqs$locus <- paste(hla.freqs$locus1, hla.freqs$locus2, sep="-")
compute.ALD(hla.freqs[hla.freqs$locus=="C-B",])

}
\keyword{datasets}
