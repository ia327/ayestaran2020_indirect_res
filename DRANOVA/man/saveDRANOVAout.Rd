\name{saveDRANOVAout}
\alias{saveDRANOVAout}
\title{Saves the ANOVA output}
\description{
Saves the ANOVA output in a CSV file
}
\usage{
saveDRANOVAout(X, file, ...)
}
\arguments{
  \item{X}{Object of \link{DRANOVA_MUT_OUT-class} or \link{DRANOVA_OUT-class} object, which can be generated with \link{runANOVA} and \link{runTissueMultiComp}, respectively.
}
  
  \item{file}{Absolute of relative path + name of the CSV file.}
  
  \item{overwrite}{Boolean input has to be either \code{"TRUE"} or \code{"FALSE"}. If set to \code{"TRUE"} it will create or replace \code{"file"}. In case of \code{"FALSE"} it will create the file only if it doesn't exists, otherwise throws an error. Default is \code{"FALSE"}.}
  
}
\seealso{
  \code{\link{runANOVA}},
  \code{\link{DRANOVA_MUT_OUT-class}},
  \code{\link{runTissueMultiComp}},
  \code{\link{DRANOVA_OUT-class}},
  \code{\link{saveContigencyTable}}
}
\examples{
# load data
data(gdsc_ActArea)  # Drug response in 1-AUC from GDSC
data(bem)           # Binary event matrix (71 oncogenes)
data(tissue)        # Cell line tissue-of-origin

# run anova with drug response from GDSC
out_gdsc = runANOVA(new("DRANOVA", 
                       res=gdsc_ActArea, 
                       bem=bem, 
                       tissue=tissue))

saveDRANOVAout(out_gdsc, file="gdsc_ANOVA.csv", overwrite=TRUE)
}
\keyword{I/O}