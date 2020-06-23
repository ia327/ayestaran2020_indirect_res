\name{gdsc_ActArea}
\docType{data}
\alias{gdsc_ActArea}
\title{gdsc drug response data measured in 1-AUC}
\description{
Contains the drug response from GDSC for 292 cell lines (rows) and 15 drugs (columns), which are the overlapping cell lines and drugs with CCLE. Row names containing the COSMIC cell identifiers and column names are the drug names. Drug response is reported in one minus normalised area under fitted drug response curve. 
}
\usage{data(gdsc_ActArea)}
\format{A matrix containing 292 cell lines and 15 drugs}
\references{
Garnett, M.J. et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature 483, 570-575 (2012).

Barretina, J. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603-607 (2012).
}
\keyword{datasets}