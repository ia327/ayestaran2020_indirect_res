# MODIFIED VERSION
## TODO: merge CFEs with identical patterns
##       filter out cell lines that weren't screened with a single drug

# ---------------------------------------------------------------------------------
# DRANOVA class definition:
# ---------------------------------------------------------------------------------
setClass("DRANOVA", 
         representation=representation(
           res="matrix", 
           bem="matrix",
           covariates="data.frame"
         ),
         validity=function(object) {
           if (is.null(rownames(object@res))) {
             stop("ERROR: Rownames in drug response matrix need to be cell IDs.")
             return(FALSE)
           }
           if (is.null(rownames(object@bem))) {
             stop("ERROR: Rownames in binary event matrix (bem) need to be cell IDs.")
             return(FALSE)
           }
           if (any(rownames(object@res) != rownames(object@bem))) {
             cat("ERROR: cell IDs in drug response matrix (rownames) and binary event matrix (bem) rownames are not mapping together.")
             return(FALSE)
           }
           if (is.null(names(object@covariates))) {
             stop("ERROR: Rowames in the covariate matrix need to be cell IDs.")
             return(FALSE)
           }
           if (any(rownames(object@covariates) != rownames(object@bem))) {
             cat("ERROR: cell IDs of covariates (rownames) and binary event matrix (=bem) rownames are not mapping together.")
             return(FALSE)
           }
           return(TRUE)
         }
)

setGeneric("runANOVA", function(X, ...) standardGeneric("runANOVA"))
setGeneric("runTissueMultiComp",function(X, ...) standardGeneric("runTissueMultiComp"))

setMethod("show", "DRANOVA",
          function(object) {
            cat("Object of class", class(object),"\n")
            cat(" Number of drugs:", ncol(object@bem),"\n")
            cat(" Number of cell lines:", nrow(object@bem),"\n")
            cat(" Number of covariates: \n")
            sapply(colnames(object@covariates), 
                   function(x) cat('\t', x, ': ',
                                   length(unique(object@covariates[,x])),
                                   '\n'))
          }
)

# ---------------------------------------------------------------------------------
# DRANOVA_OUT abstract class definition, which is the template for DRANOVA_MUT_OUT
# and DRANOVA_TISUE_OUT:
# ---------------------------------------------------------------------------------
setClass("DRANOVA_OUT", 
         representation=representation(
           label="character",
           pVal='numeric',
           qVal="numeric",
           effect="numeric",
           size="numeric")
)

setGeneric("plotVolcano", function(X, ...) standardGeneric("plotVolcano"))
setGeneric("plotComparison", function(X, Y, ...) standardGeneric("plotComparison"))
setGeneric("saveDRANOVAout", function(X, file, ...) standardGeneric("saveDRANOVAout"))
setGeneric("saveContigencyTable", function(X, Y, file, ...) standardGeneric("saveContigencyTable"))

setMethod("show", "DRANOVA_OUT",
          function(object) {
            cat("Object of class", class(object),"\n")
            cat(" Number of drug-to-oncogene tested associations:", length(object@qVal),"\n")
            qVal <- summary(object@qVal)
            cat(" q-value:\n")
            cat(" ", names(qVal), "\n")
            cat(" ", qVal, "\n")
            eff <- summary(object@effect)
            cat(" effect:\n")
            cat(" ", names(eff), "\n")
            cat(" ", eff, "\n")
            size <- summary(object@size)
            cat(" size:\n")
            cat(" ", names(size), "\n")
            cat(" ", size, "\n")
          }
)

# ---------------------------------------------------------------------------------
# DRANOVA_MUT_OUT, which inherits all slots and methods from abstract class 
# DRANOVA_OUT:
# ---------------------------------------------------------------------------------
setClass("DRANOVA_MUT_OUT", 
         representation=representation(
           pValCov='list',
           qValCov='list'),
         contains="DRANOVA_OUT"   
)



# ---------------------------------------------------------------------------------
# Method: runANOVA
# ---------------------------------------------------------------------------------
setMethod("runANOVA", "DRANOVA", function(X, cov.names=NA, n.min.mutants = 3) {
  res <- X@res
  bem <- X@bem
  covariates <- X@covariates
  
  # To avoid blank spaces and punctuation in covariate names
  cov.names <- str_replace_all(cov.names, "[[:space:]]", '.')
  cov.names <- str_replace_all(cov.names, "[[:punct:]]", '.')
  colnames(covariates) <- str_replace_all(colnames(covariates), "[[:space:]]", '.')
  colnames(covariates) <- str_replace_all(colnames(covariates), "[[:punct:]]", '.')
  
  # If cov.names is missing, use all present covariates
  if (all(!is.na(cov.names))){
    covariates <- covariates[,cov.names, drop = FALSE]
  }
  
  pvalues <- matrix(nrow=ncol(res), ncol=ncol(bem), 
                    dimnames=list(colnames(res), colnames(bem)))
  
  Cov_pvalues <- lapply(1:ncol(covariates), function(x) matrix(nrow=ncol(res), ncol=ncol(bem), 
                                                               dimnames=list(colnames(res), colnames(bem))))
  names(Cov_pvalues) <- colnames(covariates)
  
  populationSize <- matrix(nrow=ncol(res), ncol=ncol(bem), 
                           dimnames=list(colnames(res), colnames(bem)))
  effect <- matrix(nrow=ncol(res), ncol=ncol(bem), 
                   dimnames=list(colnames(res), colnames(bem)))
  
  getCohensD <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation  
    return(md/csd)
  }
  
  
  
  for(res_i in 1:ncol(res)) { # loop through drugs
    cat("drug#:",res_i,'\n')
    # Check which of the covariates have more than one level
    multiLev <- apply(covariates, 2, 
                      function(x) length(unique(x[!is.na(res[,res_i]) &
                                                    !is.na(x)])) > 1)
    used.covs <- colnames(covariates)[multiLev]
    
    for(bem_i in 1:ncol(bem)) { # loop through bem
      # Filter for minimum number of mutations
      if (sum(bem[!is.na(res[,res_i]), bem_i], na.rm=TRUE) >= n.min.mutants &
          length(bem[!is.na(res[,res_i]), bem_i]) - sum(bem[!is.na(res[,res_i]), bem_i], na.rm=TRUE) >= n.min.mutants){
        # Write the command and run it
        if (length(used.covs != 0)){
          aov.command <- paste0('aov(res[,res_i] ~ ', 
                                paste(paste0('covariates$',used.covs,''), collapse = ' + '),
                                ' + bem[,bem_i])')
        } else {
          aov.command <- 'aov(res[,res_i] ~ bem[,bem_i])'
        }
        
        fit <- eval(parse(text = aov.command))
        
        anovaResult <- anova(fit)  
        pvalues[res_i, bem_i] <- anovaResult['bem[, bem_i]','Pr(>F)']
        
        if (length(used.covs != 0)){
          for (i in 1:length(used.covs)){
            Cov_pvalues[[used.covs[i]]][res_i, bem_i] <- anovaResult[paste0('covariates$',used.covs[i]),'Pr(>F)']
          }
        }
        
        if (!is.na(pvalues[res_i, bem_i])) {
          sign <- -1
          if (mean(res[bem[,bem_i]==1, res_i], na.rm=T) > mean(res[bem[,bem_i]==0, res_i], na.rm=T)){
            sign <- 1
          }
          effect[res_i, bem_i] <- sign * getCohensD(na.omit(res[bem[,bem_i]==1, res_i]), 
                                                    na.omit(res[bem[,bem_i]==0, res_i]))
          populationSize[res_i, bem_i] <- sum(bem[!is.na(res[,res_i]), bem_i]==1, na.rm=TRUE)
        }
      }
    }
  }
  
  # remove drug-to-gene associations, 
  # where the mutant gene was never treated
  mask <- !is.na(effect)
  
  new("DRANOVA_MUT_OUT", 
      label=apply(expand.grid(rownames(pvalues), colnames(pvalues)), 1, 
                  function(x) paste(x,collapse=":"))[mask],
      pVal=as.numeric(as.vector(pvalues)[mask]),
      qVal=p.adjust(as.numeric(as.vector(pvalues)[mask]), method = 'BH'),
      pValCov=lapply(Cov_pvalues, 
                     function(x) as.numeric(as.vector(x)[mask])),
      qValCov=lapply(names(Cov_pvalues), function(x) {
        if (multiLev[x]){
          p.adjust(as.numeric(as.vector(Cov_pvalues[[x]])[mask]), method = 'BH')
        } else {
          NA[mask]
        }
      }),      
      effect=as.numeric(as.vector(effect)[mask]), 
      size=as.numeric(as.vector(populationSize)[mask]))
})

# ---------------------------------------------------------------------------------
# Method: plotVolcano
# ---------------------------------------------------------------------------------
setMethod("plotVolcano", "DRANOVA_OUT", 
          function(X, 
                   FDR=0.2,
                   main="", 
                   withLabel=TRUE, 
                   labelSig=FDR
          ) {
            if (class(X)[1] != "DRANOVA_OUT" &
                class(X)[1] != "DRANOVA_MUT_OUT")
              stop("ERROR: input X must be of class \"DRANOVA_OUT\" or \"DRANOVA_MUT_OUT\"")
            
            # format anova output
            size <- X@size+20
            labels <- X@label
            effect <- X@effect
            qVal <- X@qVal
            
            # colours of the the drug-to-oncogene associations
            col_doa <- rep(rgb(0, 151, 30, 120, maxColorValue=255), length(qVal)) # green
            col_doa[effect>=0] <- rgb(225, 0, 0, 120, maxColorValue=255) # red
            col_doa[qVal>FDR] <- rgb(200, 200, 200, 120, maxColorValue=255) # grey
            
            # Labels for the drug-to-oncogene associations  
            if (withLabel) { 
              labels[qVal>labelSig] <- ""
            } else { 
              labels <- "" 
            }
            
            # function for transforming y-axis
            reverselog_trans <- function (base = exp(1)) {
              trans <- function(x) -log(x, base)
              inv <- function(x) base^-x
              trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base), 
                        domain = c(1e-100, Inf))
            }
            
            # function for relabeling the y-axis
            scientific_10 <- function(x) {
              parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
            }
            
            # Volcano plot
            ggplot(data=data.frame(effect, qVal, col_doa, labels, size), 
                   aes(x=effect, y=qVal)) + 
              scale_y_continuous(trans=reverselog_trans(10),label=scientific_10,limits=c(1, min(qVal)/10)) +
              scale_x_continuous(breaks =seq(floor(min(effect)), ceiling(max(effect)+2), 2), 
                                 limits=c(floor(min(effect)), ceiling(max(effect)+2))) +
              geom_point(size=size/10, colour=col_doa) +
              geom_vline(xintercept=0, lwd=0.3) +
              geom_hline(yintercept=FDR, lwd=0.3, linetype=2) +
              annotate("text", x=min(effect) + 0.1, y=FDR-0.07, size=3.5, 
                       label=paste((FDR*100), "% FDR", sep=""), colour="black") +
              xlab("Drug effect (signed Cohen's d)") + 
              ylab("Adjusted P-value") +
              ggtitle(main) +
              theme_bw() +
              geom_text(mapping=aes(x=effect, 
                                    y=qVal,
                                    label=labels), colour="darkgrey", size=4, vjust=-1) 
          })


# ---------------------------------------------------------------------------------
# Method: plotComparison
# ---------------------------------------------------------------------------------
setMethod("plotComparison", "DRANOVA_OUT", 
          function(X, Y, FDR=0.2,
                   main="",
                   xlab="X signed -log10(q-value)",
                   ylab="Y signed -log10(q-value)",
                   withLabel=TRUE, labelSig=FDR) {
            if (class(X)[1] != "DRANOVA_OUT" &
                class(X)[1] != "DRANOVA_MUT_OUT")
              stop("ERROR: input X must be of class \"DRANOVA_OUT\" or \"DRANOVA_MUT_OUT\"")
            if (class(Y)[1] != "DRANOVA_OUT" &
                class(Y)[1] != "DRANOVA_MUT_OUT")
              stop("ERROR: input Y must be of class \"DRANOVA_OUT\" or \"DRANOVA_MUT_OUT\"")
            
            labels <- intersect(X@label, Y@label)
            X_idx <- match(labels, X@label)
            Y_idx <- match(labels, Y@label)
            
            col<-rep(rgb(0,0,0,20,maxColorValue=255),length(labels))
            allGreen<-rgb(0,255,100,180,maxColorValue=255)
            allRed<-rgb(255,100,0,180,maxColorValue=255)
            
            lightGreen<-rgb(0,255,100,80,maxColorValue=255)
            lightRed<-rgb(255,100,0,80,maxColorValue=255)
            
            col[which(X@qVal[X_idx]<FDR & X@effect[X_idx] < 0)]<-lightGreen
            col[which(Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] < 0)]<-lightGreen
            
            col[which(X@qVal[X_idx]<FDR & X@effect[X_idx] > 0)]<-lightRed
            col[which(Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] > 0)]<-lightRed
            
            col[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR & X@effect[X_idx] < 0 & Y@effect[Y_idx] <0)]<-allGreen
            col[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR & X@effect[X_idx] > 0 & Y@effect[Y_idx] >0)]<-allRed
            
            bcol<-rep(NA,length(labels))
            bcol[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR & X@effect[X_idx] < 0 & Y@effect[Y_idx] <0)]<-'black'
            bcol[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR & X@effect[X_idx] > 0 & Y@effect[Y_idx] >0)]<-'black'
            
            x <- -log10(X@qVal[X_idx]) * sign(X@effect[X_idx])
            y <- -log10(Y@qVal[Y_idx]) * sign(Y@effect[Y_idx])
            plot(x, 
                 y,
                 main=main,
                 cex=1.4, bg=col, pch=21,col=bcol,
                 xlab=xlab,ylab=ylab,frame.plot=F, 
                 xlim=c(min(x)-3, max(x)+3),
                 ylim=c(min(y)-5, max(y)+5))
            
            abline(h=0,col='gray')
            abline(v=0,col='gray')
            
            abline(h=-log10(FDR),col='gray',lty=2)
            abline(h=log10(FDR),col='gray',lty=2)
            
            abline(v=-log10(FDR),col='gray',lty=2)
            abline(v=log10(FDR),col='gray',lty=2)
            
            # Labels for the drug-to-oncogene associations  
            labels <-X@label[X_idx]
            if (withLabel) { 
              labels[X@qVal[X_idx]>labelSig] <- ""
              labels[Y@qVal[Y_idx]>labelSig] <- ""
            } else { 
              labels <- "" 
            }  
            text(x = -log10(X@qVal[X_idx]) * sign(X@effect[X_idx]),cex=0.7,
                 y = -log10(Y@qVal[Y_idx]) * sign(Y@effect[Y_idx])+0.5,labels=labels)
            
            # legend in bottom right
            X_Sidx<-X_idx[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR)]
            Y_Sidx<-Y_idx[which(X@qVal[X_idx]<FDR & Y@qVal[Y_idx]<FDR)]
            CT<-cor.test(-log10(X@pVal[X_Sidx]) * sign(X@effect[X_Sidx]),method = 'pearson',
                         -log10(Y@pVal[Y_Sidx]) * sign(Y@effect[Y_Sidx]))
            legend('bottomright',c(paste('pc =',format(CT$estimate,digits = 2)),
                                   paste('p =',format(CT$p.value,digits=2,scientific=TRUE))), cex=0.7)
          })

# ---------------------------------------------------------------------------------
# Method: saveDRANOVAout
# ---------------------------------------------------------------------------------
setMethod("saveDRANOVAout", "DRANOVA_OUT", 
          function(X, file, overwrite=F) {
            if (!overwrite && file.exists(file))
              stop(paste("ERROR: \"", file, "\" already exists", sep=""))
            # Add covariate p values and q values to the table
            covs <- names(X@pValCov)
            cov.res.names <- as.vector(sapply(covs, function(x) paste0(c('p.value.',
                                                                         'q.value.'),
                                                                       x)))
            head <- t(c("association", 
                        "p.value",
                        "corrected p.value (=q.value)", 
                        "effect", 
                        "n.samples",
                        cov.res.names))
            cov.matrix <- matrix(nrow = length(X@label), ncol = length(covs) * 2)
            for (i in 1:length(covs)){
              index <- 2*(i-1) + 1
              cov.matrix[,index] <- X@pValCov[[i]]
              cov.matrix[,(index + 1)] <- X@qValCov[[i]]
            }
            write.table(head, file=file, sep = ",", row.names=F, col.names=F)
            write.table(cbind(X@label,X@pVal, X@qVal, X@effect, X@size, cov.matrix), 
                        file=file, append=TRUE, sep = ",", row.names=F, col.names=F)
          })

# ---------------------------------------------------------------------------------
# Method: saveContigencyTable
# ---------------------------------------------------------------------------------
setMethod("saveContigencyTable", "DRANOVA_OUT", 
          function(X, Y, file,
                   xlab="X signed -log10(q-value)",
                   ylab="Y signed -log10(q-value)",
                   FDR=0.2, overwrite=F) {
            if (!overwrite && file.exists(file))
              stop(paste("ERROR: \"", file, "\" already exists", sep=""))
            
            labels <- intersect(X@label, Y@label)
            X_idx <- match(labels, X@label)
            Y_idx <- match(labels, Y@label)
            
            Ftable<-matrix(0,3,3)
            Ftable[1,1]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx]<0 &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx]>0))
            Ftable[1,2]<-length(which(X@qVal[X_idx]>=FDR &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx]>0))
            Ftable[1,3]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx] > 0 &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] > 0))
            Ftable[2,1]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx] < 0 &
                                        Y@qVal[Y_idx]>=FDR))
            Ftable[2,2]<-length(which(X@qVal[X_idx]>=FDR &
                                        Y@qVal[Y_idx]>=FDR))
            Ftable[2,3]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx] > 0 &
                                        Y@qVal[Y_idx]>=FDR))
            Ftable[3,1]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx] < 0 &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] < 0))
            Ftable[3,2]<-length(which(X@qVal[X_idx]>=FDR &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] < 0))
            Ftable[3,3]<-length(which(X@qVal[X_idx]<FDR & X@effect[X_idx] > 0 &
                                        Y@qVal[Y_idx]<FDR & Y@effect[Y_idx] < 0))
            rownames(Ftable) <- c("resistant", "non-significant", "sensitive")
            colnames(Ftable) <- c("sensitive", "non-significant", "resistant")
            
            head <- paste("# Contingency Table\n", 
                          "# Fisher's exact test p-value = ", fisher.test(Ftable)$p.value, "\n",
                          "# Columns containing ", xlab, "\n",
                          "# Rows containing ", ylab, "\n\n", 
                          ",", paste(colnames(Ftable), collapse=","), sep="")
            
            write(head, file=file)
            write.table(Ftable, file=file, sep = ",", row.names=T, col.names=F, append=T)
          })

# # ---------------------------------------------------------------------------------
# # Method: tissueMultiComparison
# # ---------------------------------------------------------------------------------
# setMethod("runTissueMultiComp", "DRANOVA", 
# function(X, Y, FDR=0.2) {
#   if (class(X)[1] != "DRANOVA")
#     stop("ERROR: input X must be of class \"DRANOVA\"")
#   if (class(Y)[1] != "DRANOVA_MUT_OUT")
#     stop("ERROR: input Y must be of class \"DRANOVA_MUT_OUT\"")
#   
#   id<-which(Y@qValTissue<FDR)
#   splittedLabels<-unlist(str_split(Y@label[id],':'))
#   DRUGS<-unique(splittedLabels[seq(1,length(splittedLabels),2)])
#   
#   print(paste(length(DRUGS),'drugs have a tissue qval <', FDR, 'in at least one MANOVA model'))
#   
#   TISSUE_factor<-unique(tissue)
#   TISSUE_factor<-setdiff(TISSUE_factor,'other')
#   
#   pvalues <- matrix(nrow=length(DRUGS), ncol=length(TISSUE_factor), 
#                              dimnames=list(DRUGS, TISSUE_factor))
#   
#   populationSizes <- matrix(nrow=length(DRUGS), ncol=length(TISSUE_factor), 
#                             dimnames=list(DRUGS, TISSUE_factor))
#   effectSizes <- matrix(nrow=length(DRUGS), ncol=length(TISSUE_factor), 
#                         dimnames=list(DRUGS, TISSUE_factor))
#   
#   getCohensD <- function(x, y) {
#     lx <- length(x)- 1
#     ly <- length(y)- 1
#     md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
#     csd <- lx * var(x) + ly * var(y)
#     csd <- csd/(lx + ly)
#     csd <- sqrt(csd)                     ## common sd computation  
#     return(md/csd)
#   }
#   
#   for(drug in DRUGS) { # loop through drugs
#     print(paste("drug#:",drug))
#     for(TISSUE in TISSUE_factor) { # loop through tissues   
#         Val_from_tissue<-X@res[tissue==TISSUE,drug]
#         Val_from_other<-X@res[tissue!=TISSUE,drug]
#         n_from_tissue<-length(which(!is.na(Val_from_tissue)))
#         n_from_other<-length(which(!is.na(Val_from_other)))
#         
#         if(n_from_tissue>1 & n_from_other>1){
#           fit <- t.test(X@res[,drug]~(tissue==TISSUE))
#           pvalues[drug, TISSUE] <- fit$p.value
#           populationSizes[drug,TISSUE]<-n_from_tissue
#           effectSizes[drug,TISSUE]<-sign(fit$estimate[1]-fit$estimate[2])*
#             getCohensD(Val_from_tissue[!is.na(Val_from_tissue)],Val_from_other[!is.na(Val_from_other)])
#         }    
#     }
#   }
#   mask <- !is.na(pvalues)
#   
#   new("DRANOVA_OUT", 
#       label=apply(expand.grid(rownames(pvalues), colnames(pvalues)), 1, 
#                   function(x) paste(x,collapse=":"))[mask],
#       
#       pVal=as.vector(pvalues)[mask],
#       qVal=p.adjust(as.vector(pvalues)[mask],method = 'BH'),
#       size=as.vector(populationSizes)[mask],
#       effect=as.vector(effectSizes)[mask]
#   )
# })
