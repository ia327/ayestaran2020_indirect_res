\name{runTissueMultiComp}
\alias{runTissueMultiComp}
\title{Runs a tissue multi-comparison}
\description{
Evaluates which tissue-of-origin is significantly associated with a drug. A high correlation between tissue and drug responses are expected. 
}
\usage{
runTissueMultiComp(X, Y, FDR=0.2)
}
\arguments{
  \item{X}{
     Object of \link{DRANOVA-class}, which contains the drug response, mutational status (bem) and tissue-of-origin for each cell line.}
  \item{Y}{
     Object of \link{DRANOVA_MUT_OUT-class}, which contains the output from \link{runANOVA}.}
  \item{FDR}{False discovery rate threshold (default = 0.2).}
}
\references{
Garnett, M.J. et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. Nature 483, 570-575 (2012).

Barretina, J. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603-607 (2012).
}
\seealso{
  \code{\link{DRANOVA-class}},
  \code{\link{DRANOVA_MUT_OUT-class}},
  \code{\link{plotVolcano}},
  \code{\link{saveDRANOVAout}},
  \code{\link{plotComparison}},
  \code{\link{saveContigencyTable}}
}
\examples{
# load data
data(bem)           # Binary event matrix (71 oncogenes)
data(gdsc_ActArea)  # Drug response from GDSC
data(ccle_ActArea)  # Drug response from CCLE
data(tissue)        # Cell line tissue-of-origin

# for GDSC: (1) generate DRANOVA input, (2) run anova and (3) run tissue multi comparison.
gdsc_actArea = new("DRANOVA", res=gdsc_ActArea, bem=bem, tissue=tissue)
gdsc_actArea_out = runANOVA(gdsc_actArea)
gdsc_tissue_Act_out<-runTissueMultiComp(gdsc_actArea, gdsc_actArea_out)

# for CCLE: same steps as outlined for GDSC above
ccle_actArea = new("DRANOVA", res=ccle_ActArea, bem=bem, tissue=tissue)
ccle_actArea_out = runANOVA(ccle_actArea)
ccle_tissue_Act_out<-runTissueMultiComp(ccle_actArea, ccle_actArea_out)

# Plotting concordance from CCLE and GDSC by focusing on the tissue-of-origin
plotComparison(gdsc_tissue_Act_out, ccle_tissue_Act_out, 
               main="CCLE & GDSC tissue associations comparison with Activity curve",
               xlab="GDSC tissue assoc signed -log10(q-value)",
               ylab="CCLE tissue assic signed -log10(q-value)")
}
\keyword{ANOVA, tissue multi-comparison }
