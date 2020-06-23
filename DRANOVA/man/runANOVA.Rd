\name{runANOVA}
\alias{runANOVA}
\title{Runs a drug ANOVA}
\description{
Runs for each drug an ANOVA, while considering the mutational status as input and tissue-of-origin as covariate. 
}
\usage{
runANOVA(X)
}
\arguments{
  \item{X}{
     Object of \link{DRANOVA-class}, which contains the drug response, mutational status (bem) and tissue-of-origin for each cell line.}
}
\details{
Internally we use the \link{aov} function from R package \link{stats} with the formula "res ~ tissue + gene". "res" is the drug response, tissue is a covariate (e.g. breast, lung, etc.) and gene is a binary value (wild type or mutant). \link{aov} performs a sequential ANOVA and first fits tissue, and than the mutational status. Tissue on its own is a powerful but trivial predictor, e.g. treating a lung cancer patient with a lung cancer drug. We are focusing on molecular descriptors which are stratifying variance in drug response with previously removing the impact of tissue specificity.   
}
\references{
Garnett, M.J. et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature 483, 570-575 (2012).
}
\seealso{
  \code{\link{DRANOVA-class}},
  \code{\link{plotVolcano}},
  \code{\link{saveDRANOVAout}},
  \code{\link{plotComparison}},
  \code{\link{saveContigencyTable}}
}
\examples{
# load data
data(bem)             # Binary event matrix (71 oncogenes)
data(gdsc_transIC50)  # Drug response in -log10(IC50/10^6)
data(tissue)          # Cell line tissue-of-origin

# run anova
output <- runANOVA(new("DRANOVA", 
                       res=gdsc_transIC50, 
                       bem=bem, 
                       tissue=tissue))

# Volcano plot of ANOVA
plotVolcano(output, 
            main="ANOVA with IC50s from GDSC",
            withLabel=TRUE,
            labelSig=0.01)
}
\keyword{ANOVA}
