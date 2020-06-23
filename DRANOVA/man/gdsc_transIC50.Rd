\name{gdsc_transIC50}
\docType{data}
\alias{gdsc_transIC50}
\title{gdsc drug response data measured in -log10(IC50/10^6)}
\description{
Contains the drug response from GDSC for 292 cell lines (rows) and 15 drugs (columns), which are the overlapping cell lines and drugs with CCLE. Row names containing the COSMIC cell identifiers and column names are the drug names. Drug response is reported in IC50, the half maximal inhibitory concentration of a substance with respect to cell viability. Those IC50s are transfromed as in Haibe-Kains et al, which is -log10(IC50/10^6).
}
\usage{data(gdsc_transIC50)}
\format{A matrix containing 292 cell lines and 15 drugs}
\references{
Garnett, M.J. et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature 483, 570-575 (2012).

Barretina, J. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603-607 (2012).

Haibe-Kains, B. et al. Inconsistency in large pharmacogenomic studies. Nature 504, 389-93 (2013).
}
\keyword{datasets}