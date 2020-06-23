\name{DRANOVA_MUT_OUT-class}
\Rdversion{1.1}
\docType{class}
\alias{DRANOVA_MUT_OUT-class}
\alias{plotComparison,DRANOVA_MUT_OUT-method}
\alias{plotVolcano,DRANOVA_MUT_OUT-method}
\alias{saveDRANOVAout,DRANOVA_MUT_OUT-method}
\alias{saveContigencyTable,DRANOVA_MUT_OUT-method}
\alias{show,DRANOVA_MUT_OUT-method}

\title{Class \code{"DRANOVA_MUT_OUT"}}
\description{
Output from \link{runANOVA}, which can be visualised with \link{plotVolcano} and stored with \link{saveDRANOVAout}. For comparing two \code{"DRANOVA_MUT_OUT"} objects the methods  \link{plotComparison} and \link{saveContigencyTable} could be applied.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DRANOVA_MUT_OUT", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{label}:}{Vector conataining the drug-to-oncogene label.}
    \item{\code{pVal}:}{Vector of the associated p-values for each drug-to-oncogene association.}
    \item{\code{qVal}:}{\code{"pVal"} correct for multiple hypothesis with Benjamini & Hochberg.}
    \item{\code{pValTissue}:}{Vector of the associated p-values for each drug-to-tissue association.}
    \item{\code{qValTissue}:}{\code{"qValTissue"} correct for multiple hypothesis with Benjamini & Hochberg the tissue p-value.}
    \item{\code{effect}:}{Effect size of drug-to-oncogene association with signed Cohen's D.}
    \item{\code{size}:}{Vector containing the number of cell lines carrying the oncogenes of above outlined drug-to-oncogene associations.}
  }
}
\section{Methods}{
  \describe{
    \item{plotComparison}{\code{signature(X = "DRANOVA_MUT_OUT")}: ... }
    \item{plotVolcano}{\code{signature(X = "DRANOVA_MUT_OUT")}: ... }
    \item{saveDRANOVAout}{\code{signature(X = "DRANOVA_MUT_OUT")}: ... }
    \item{saveContigencyTable}{\code{signature(X = "DRANOVA_MUT_OUT")}: ... }
    \item{show}{\code{signature(object = "DRANOVA_MUT_OUT")}: ... }
	 }
}

\seealso{
  \code{\link{runANOVA}},
  \code{\link{plotVolcano}},
  \code{\link{saveDRANOVAout}},
  \code{\link{plotComparison}},
  \code{\link{saveContigencyTable}},
  \code{\link{show}}
}
\examples{
# load data
data(ccle_ActArea)  # Drug response in 1-AUC from CCLE
data(bem)           # Binary event matrix (71 oncogenes)
data(tissue)        # Cell line tissue-of-origin

obj <- runANOVA(new("DRANOVA", 
                    res=ccle_ActArea, 
                    bem=bem, 
                    tissue=tissue))
show(obj)
}
\keyword{classes}
