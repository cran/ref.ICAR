#' OLM and ICAR model probabilities for areal data
#' @description Performs simultaneous selection of covariates and spatial model structure for areal data.
#' @references
#' \insertRef{Porter_2023}{ref.ICAR}
#'
#' @author Erica M. Porter, Christopher T. Franck, and Marco A.R. Ferreira
#'
#' @importFrom spdep sf hier.part pracma stats
#' @param X A matrix of covariates, which should include a column of 1's for models with a non-zero intercept
#' @param Y A vector of responses.
#' @param b Training fraction for the fractional Bayes factor (FBF) approach.
#' @param H Neighborhood matrix for spatial subregions.
#' @param H.spectral Spectral decomposition of neighborhood matrix, if user wants to pre-compute it to save time.
#' @param Sig_phi Pseudo inverse of the neighborhood matrix, if user wants to pre-compute it to save time.
#' @param verbose If FALSE, marginal likelihood progress is not printed.
#'
#' @return A list containing a data frame with all posterior model probabilities and other selection information.
#'     \item{probs.mat}{Data frame containing posterior model probabilities for all candidate OLMs and ICAR models from the data.}
#'     \item{mod.prior}{Vector of model priors used to obtain the posterior model probabilities.}
#'     \item{logmargin.all}{Vector of all (log) fractional integrated likelihoods.}
#'     \item{base.model}{Maximum (log) fractional integrated likelihood among all candidate models.  All fractional Bayes factors are obtained with respect to this model.}
#'     \item{BF.vec}{Vector of fractional Bayes factors for all candidate models.}
#' @export

probs.icar <- function(Y,X,H,H.spectral=NULL,Sig_phi=NULL,b=0.05,verbose=FALSE) {

  ## Predefine tauc grid for obtaining logmax value
  tauc.grid <- c(seq(0.001,10,by=.01),seq(10.1,100,by=10))
  X <- as.matrix(X)
  n <- length(Y)
  p <- ncol(X)-1
  num.mods <- 2^(p+1)

  ind.logmargin <- rep(NA,(num.mods/2))
  spatial.logmargin <- rep(NA,(num.mods/2))
  spatial.logmargin.1 <- rep(NA,(num.mods/2))
  spatial.logmargin.b <- rep(NA,(num.mods/2))
  ind.BF <- rep(NA,(num.mods/2))
  spatial.BF <- rep(NA,(num.mods/2))
  ind.mod.probs <- rep(NA,(num.mods/2))
  spatial.mod.probs <- rep(NA,(num.mods/2))

  if(p>1){
    ## dX is a matrix of 0's and 1's indicating all possible covariate combinations (no interactions)
    x <- rbind(rep(0,p),combos(p)$binary)
    dX <- cbind(rep(1,2^p),x)

    ## Create a vector of formula labels for each candidate model
    predictor <- lapply(1:2^p, matrix, data=NA)
    predictor.labels <- lapply(1:2^p,matrix, data=NA)
    predictor.formulas <- lapply(1:2^p,matrix, data=NA)
    formula.vec <- c(rep(NA, 2^p))
    labels <- c(rep(NA,p))
    for(i in 1:(p)){labels[i] <- paste("X",i,sep="")}
    labels <- c("Intercept",labels)

    for(i in 1:2^p){
      predictor[[i]] <- diag(dX[i,])
      colnames(predictor[[i]]) <- labels
      predictor.labels[[i]]<- c(rownames(as.matrix((which(colSums(predictor[[i]]) != 0)))))
      predictor.formulas[[i]] <- formula(paste("Y ~ ", paste(predictor.labels[[i]], collapse=" + ")))
      formula.vec[i] <- Reduce(paste,deparse(predictor.formulas[[i]]))
    }

    mod.prior <- rep(NA,2^p)
    for(i in 1:2^p){
      k <- (length(which(colSums(predictor[[i]]) != 0)) - 1)
      mod.prior[i] <- ((1/(p+1))*(nchoosek(p,k))^(-1))/2
    }

    mod.prior <- c(mod.prior,mod.prior)

    ## combo.list is a list of matrices with the actual covariate values for each combination in dX
    ## Consider the case of p>1, p=1, and p=0 because the combos function only accommodates p>1
    combo.list <- lapply(1:2^p, matrix, data= NA)
    for(i in 1:2^p){
      combo.list[[i]] <- X%*%diag(dX[i,])
    }
    for(i in 1:((2^p)-1)) {combo.list[[i]] <- as.matrix(combo.list[[i]][,-(which(colSums(combo.list[[i]]) == 0))])}
  }

  else if(p==1) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))
    formula.vec[2] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept + X1", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[2]] <- X
    combo.list[[1]] <- as.matrix(X[,1])

    mod.prior <- rep(NA,2^p)
    for(i in 1:2^p){
      k <- 1
      mod.prior[i] <- ((1/(p+1))*(nchoosek(p,k))^(-1))/2
    }

    mod.prior <- c(mod.prior,mod.prior)
  } else if(p==0) {
    formula.vec <- c(rep(NA, 2^p))
    formula.vec[1] <- Reduce(paste,deparse(formula(paste("Y ~ ", paste("Intercept", collapse=" + ")))))

    combo.list <- lapply(1:2^p, matrix, data= NA)
    combo.list[[1]] <- X
    mod.prior <- c(0.5,0.5)
  }

  spatial.joint <- array(dim=c(2^p,2,length(tauc.grid)))
  logmax.value <- array(dim=c(2^p,2))

    for(i in 1:(num.mods/2)){
      ind.logmargin[i] <- vague.logmargin.independent.b(b=b, X=combo.list[[i]], Y=Y)

      spatial.joint[i,1,] <- ref.logmax.spatial.b(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                                  H.spectral=H.spectral,Sig_phi=Sig_phi,b=1)

      spatial.joint[i,2,] <- ref.logmax.spatial.b(tauc.grid,Y=Y,X=combo.list[[i]],H=H,
                                                  H.spectral=H.spectral,Sig_phi=Sig_phi,b=b)

      logmax.value[i,1] <- max(spatial.joint[i,1,])
      logmax.value[i,2] <- max(spatial.joint[i,2,])

      spatial.logmargin.1[i] <- log(integrate(ref.integrand.spatial.b,lower=0,upper=Inf,
                                              Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                              Sig_phi=Sig_phi,b=1,logmax=logmax.value[i,1])$value) + logmax.value[i,1]

      spatial.logmargin.b[i] <- log(integrate(ref.integrand.spatial.b,lower=0,upper=Inf,
                                              Y=Y,X=combo.list[[i]],H=H,H.spectral=H.spectral,
                                              Sig_phi=Sig_phi,b=b,logmax=logmax.value[i,2])$value) + logmax.value[i,2]

      spatial.logmargin[i] <- spatial.logmargin.1[i] - spatial.logmargin.b[i]
      if(verbose==TRUE){
      print(paste('Ref Integrated Likelihood ', i,' of ', num.mods/2, ' complete',' at ',date(),sep=''))
        }
    }

  logmargin.all <- c(ind.logmargin,spatial.logmargin)
  #base.model <- logmargin.all[which.max( abs(logmargin.all) )]
  base.model <- max(logmargin.all)

  for(i in 1:(num.mods/2)){
    ind.BF[i] <- exp(ind.logmargin[i]-base.model)
    spatial.BF[i] <- exp(spatial.logmargin[i]-base.model)
  }

  BF.vec <- c(ind.BF,spatial.BF)
  BF.sum.adj <- BF.vec %*% mod.prior

  mod.probs.all <- rep(NA,num.mods)
  for (i in 1:num.mods){
    mod.probs.all[i] <- (BF.sum.adj)^(-1)*BF.vec[i]*mod.prior[i]
  }

  probs.mat <- data.frame(mod.probs.all,c(rep("Independent",num.mods/2),rep("Spatial",num.mods/2)),c(formula.vec,formula.vec))
  names(probs.mat) <- c("model prob","model type","model form")

  return(list(probs.mat=probs.mat,
              mod.prior=mod.prior,
              logmargin.all=logmargin.all,
              base.model=base.model,
              BF.vec=BF.vec))
}
