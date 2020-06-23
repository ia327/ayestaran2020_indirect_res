\name{DRANOVA-package}
\alias{DRANOVA-package}
\alias{DRANOVA}
\docType{package}
\title{
\strong{DRANOVA} - \strong{D}rug \strong{R}esistance \strong{ANOVA}
}

\description{
This package runs an ANOVA for estimating the drug response variance determined by mutational status of oncogenes. The tissue-of-origin is used as a covariate. Furthermore, this package contains a Volcano plot, comparison plot of two drug response ANOVA outputs, and methods for saving the output. 
}

\details{
\tabular{ll}{
Package: \tab DRANOVA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-18\cr
License: \tab GPL v3\cr
Depends: \tab methods\cr
}
}
\author{
Michael P Menden, Francesco Iorio, Michael Schubert, Julio Saez-Rodriguez

Maintainer: Michael P Menden <\email{menden@ebi.ac.uk}>
}
\references{
Garnett, M.J. et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature 483, 570-575 (2012).

Barretina, J. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603-607 (2012).
}
\keyword{ package }
\seealso{
  \code{\link{runANOVA}},
  \code{\link{plotVolcano}},
  \code{\link{saveDRANOVAout}},
  \code{\link{plotComparison}},
  \code{\link{saveContigencyTable}}
}
\examples{
# load data
data(ccle_ActArea)  # Drug response in ActiveArea from CCLE
data(gdsc_ActArea)  # Drug response in 1-AUC from GDSC
data(bem)           # Binary event matrix (71 oncogenes)
data(tissue)        # Cell line tissue-of-origin

# run anova with drug response from CCLE
out_ccle = runANOVA(new("DRANOVA", 
                        res=ccle_ActArea, 
                        bem=bem, 
                        tissue=tissue))

# Volcano plot of CCLE
plotVolcano(out_ccle, 
            main="CCLE",
            withLabel=TRUE,
            labelSig=0.01)

# run anova with drug response from GDSC
out_gdsc = runANOVA(new("DRANOVA", 
                       res=gdsc_ActArea, 
                       bem=bem, 
                       tissue=tissue))

# Volcano plot of GDSC
plotVolcano(out_gdsc, 
            main="GDSC",
            withLabel=TRUE,
            labelSig=0.01)

# Comparison plot between CCLE and GDSC
plotComparison(out_gdsc, out_ccle,
               main="CCLE & GDSC comparison with Active Area",
               xlab="GDSC ANOVA signed -log10(q-value)",
               ylab="CCLE ANOVA signed -log10(q-value)")
               
# Saves the contigency table of 
saveContigencyTable(out_gdsc, out_ccle, 
                    file="ccle_gdsc_comparison.csv")

# Saves the ANOVA output in CSV files
saveDRANOVAout(out_gdsc, file="gdsc_ANOVA.csv")
saveDRANOVAout(out_ccle, file="ccle_ANOVA.csv")
}
