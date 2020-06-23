\name{saveContigencyTable}
\alias{saveContigencyTable}
\title{Saves the contigency table}
\description{
Saves the contigency table, showing the concordance of two different datasets, using the ANOVA analysis. The Fisher exact test p-value for independence of columns and rows is also stored in the CSV file. 
}
\usage{
saveContigencyTable(X, Y, file, ...)
}
\arguments{
   \item{X}{
  Object of \link{DRANOVA_MUT_OUT-class} or \link{DRANOVA_OUT-class} object, which can be generated with \link{runANOVA} and \link{runTissueMultiComp}, respectively. It will be compared to \code{"Y"}.}
  
  \item{Y}{
  Object of \link{DRANOVA_MUT_OUT-class} or \link{DRANOVA_OUT-class} object, which will be compared to \code{"X"}. }
  
  \item{file}{
  Absolute of relative path + name of the CSV file.}
  
  \item{overwrite}{Boolean input has to be either \code{"TRUE"} or \code{"FALSE"}. If set to \code{"TRUE"} it will create or replace \code{"file"}. In case of \code{"FALSE"} it will create the file only if it doesn't exists, otherwise throws an error. Default is \code{"FALSE"}.}
  
  \item{FDR}{False discovery rate (FDR) threshold, which is equivalent to Benjamini & Hochberg correction. Default is 0.2 (=20\% FDR).}
  
  \item{xlab}{Name of object \code{"X"}.}
  
  \item{ylab}{Name of object \code{"Y"}.}
 
 \item{\dots}{...} 
}
\seealso{
  \code{\link{runANOVA}},
  \code{\link{DRANOVA_MUT_OUT-class}},
  \code{\link{runTissueMultiComp}},
  \code{\link{DRANOVA_OUT-class}},
  \code{\link{saveDRANOVAout}}
}
\examples{
# load data
data(gdsc_ActArea)  # Drug response in 1-AUC from GDSC
data(ccle_ActArea)  # Drug response in ActiveArea from CCLE
data(bem)           # Binary event matrix (71 oncogenes)
data(tissue)        # Cell line tissue-of-origin

# run anova with drug response from GDSC
out_gdsc = runANOVA(new("DRANOVA", 
                       res=gdsc_ActArea, 
                       bem=bem, 
                       tissue=tissue))
                       
# run anova with drug response from GDSC
out_ccle = runANOVA(new("DRANOVA", 
                       res=ccle_ActArea, 
                       bem=bem, 
                       tissue=tissue))

# saves the contigency table in ccle_gdsc_comparison.csv
saveContigencyTable(out_gdsc, 
                    out_ccle, 
                    file="ccle_gdsc_comparison.csv", 
                    overwrite=TRUE)
}
\keyword{I/O}