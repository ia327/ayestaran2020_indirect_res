\name{DRANOVA-class}
\Rdversion{1.1}
\docType{class}
\alias{DRANOVA-class}
\alias{runANOVA,DRANOVA-method}
\alias{runTissueMultiComp,DRANOVA-method}
\alias{show,DRANOVA-method}

\title{Class \code{"DRANOVA"}}
\description{
This class specifies how the input for \code{\link{runANOVA}} has to look like. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DRANOVA", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{res}:}{Drug response \code{"matrix"} of cell lines (rows) and drugs (columns). \code{"rownames"} need to be cell IDs and \code{"colnames"} should contain drug names / IDs}
    \item{\code{bem}:}{Binary event \code{"matrix"} of cell lines (rows) and mutational status of genes (columns). \code{"rownames"} need to be cell IDs and \code{"colnames"} should contain gene names / IDs}
    \item{\code{tissue}:}{A vector containing the tissue names / IDs. \code{"names"} of \code{"tissue"} has to fit with cell line IDs}
  }
}
\section{Methods}{
  \describe{
    \item{runANOVA}{\code{signature(X = "DRANOVA")}: ... }
    \item{runTissueMultiComp}{\code{signature(X = "DRANOVA")}: ... }
    \item{show}{\code{signature(object = "DRANOVA")}: ... }
	 }
}
\seealso{
 \code{\link{runANOVA}}, 
 \code{\link{show}}
}
\examples{
# load data
data(ccle_ActArea)  # Drug response in 1-AUC from CCLE
data(bem)           # Binary event matrix (71 oncogenes)
data(tissue)        # Cell line tissue-of-origin

obj <- new("DRANOVA", 
            res=ccle_ActArea, 
            bem=bem, 
            tissue=tissue)
show(obj)
}
\keyword{classes}
