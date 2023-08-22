## ----setup, include=FALSE,warning=F-------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#knitr::opts_chunk$set(fig.width=8, fig.height=6)
library(rcrossref)

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=58)----

library(sf)

system.path <- system.file("extdata", "us.shape48.shp", package = "ref.ICAR", mustWork = TRUE)
shp.layer <- gsub('.shp','',basename(system.path))
shp.path <- dirname(system.path)

us.shape48 <- st_read(dsn = path.expand(shp.path), layer = shp.layer, quiet = TRUE)

## ----echo = T, eval = T, warning=F, message=F, tidy=TRUE, tidy.opts=list(width.cutoff=58)----

library(utils)

data.path <- system.file("extdata", "states-sats48.txt", package = "ref.ICAR", mustWork = TRUE)

sats48 <- read.table(data.path, header = T)
us.shape48$verbal <- sats48$VERBAL
us.shape48$percent <- sats48$PERCENT

## ----echo = T, eval = T,fig.cap="Figure 1: Observed Verbal SAT Scores", warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.align="center",fig.pos="h",fig.width=7,fig.height=5.5----

library(ggplot2)
library(classInt)
library(dplyr)

breaks_qt <- classIntervals(c(min(us.shape48$verbal) - .00001, us.shape48$verbal), n = 7, style = "quantile")

us.shape48_sf <- mutate(us.shape48, score_cat = cut(verbal, breaks_qt$brks))
ggplot(us.shape48_sf) + 
    geom_sf(aes(fill=score_cat)) +
    scale_fill_brewer(palette = "OrRd") + 
    labs(title="Plot of observed \n verbal SAT scores") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          axis.title=element_text(size=25,face="bold"),
          plot.title = element_text(face="bold", size=25, hjust=0.5)) +
    guides(fill=guide_legend("Verbal score"))

## ----echo = T, eval = T,fig.cap="Figure 2: Percent of eligible students taking the SAT", warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.align="center",fig.pos="h",fig.width=7,fig.height=5.5----

breaks_qt <- classIntervals(c(min(us.shape48$percent) - .00001, us.shape48$percent), n = 7, style = "quantile")

us.shape48_sf <- mutate(us.shape48, pct_cat = cut(percent, breaks_qt$brks))
ggplot(us.shape48_sf) + 
    geom_sf(aes(fill=pct_cat)) +
    scale_fill_brewer(palette = "OrRd") + 
    labs(title="Plot of observed \n percent SAT takers") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          axis.title=element_text(size=25,face="bold"),
          plot.title = element_text(face="bold", size=25, hjust=0.5)) +
    guides(fill=guide_legend("Percent taking"))

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60)----

library(ref.ICAR)
library(spdep)

shp.data <- shape.H(system.path)
H <- shp.data$H

class(shp.data$map)
length(shp.data$map)

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60)----

Y <- sats48$VERBAL
x <- sats48$PERCENT
X <- cbind(1,x)

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60)----
library(mvtnorm)
set.seed(3456)

ref.SAT <- ref.MCMC(y=Y,X=X,H=H,iters=15000,burnin=5000,verbose=FALSE)

names(ref.SAT)

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.width=7,fig.height=5----
par(mfrow=c(2,2))
ref.plot(ref.SAT$MCMCchain,X,burnin=5000,num.reg=length(Y))

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60)----
library(coda)

summary.params <- ref.summary(MCMCchain=ref.SAT$MCMCchain,tauc.MCMC=ref.SAT$tauc.MCMC,sigma2.MCMC=ref.SAT$sigma2.MCMC,beta.MCMC=ref.SAT$beta.MCMC,phi.MCMC=ref.SAT$phi.MCMC,accept.phi=ref.SAT$accept.phi,accept.sigma2=ref.SAT$accept.sigma2,accept.tauc=ref.SAT$accept.tauc,iters=15000,burnin=5000)

names(summary.params)
summary.params

## ----echo = T, eval = T,fig.cap="Figure 3: Posterior Medians for Verbal SAT", warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.align="center",fig.pos="h",fig.width=7,fig.height=5.5----

summary.region <- reg.summary(ref.SAT$MCMCchain,X,Y,burnin=5000)

us.shape48$verbalfits <- summary.region$reg.medians

breaks_qt <- classIntervals(c(min(us.shape48$verbalfits) - .00001, us.shape48$verbalfits), n = 7, style = "quantile")

us.shape48_sf <- mutate(us.shape48, reg_cat = cut(verbalfits, breaks_qt$brks))
ggplot(us.shape48_sf) + 
    geom_sf(aes(fill=reg_cat)) +
    scale_fill_brewer(palette = "OrRd") + 
    labs(title="Plot of fitted \n verbal SAT scores") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          axis.title=element_text(size=25,face="bold"),
          plot.title = element_text(face="bold", size=25, hjust=0.5)) +
    guides(fill=guide_legend("Region medians"))

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.width=7,fig.height=5----

### The SAT scores and percent of students are already arranged by state alphabetically
x.reg.names <- us.shape48$NAME
y.reg.names <- us.shape48$NAME

set.seed(3456)
par(mfrow=c(2,2))
sat.analysis <- ref.analysis(system.path,X,Y,x.reg.names,y.reg.names,shp.reg.names = NULL,iters=15000,burnin=5000,verbose = FALSE,tauc.start=.1,beta.start=-1,sigma2.start=.1,step.tauc=0.5,step.sigma2=0.5)

names(sat.analysis)

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=58)----

library(pracma)
library(stats)
library(spdep)

# read in the data as contained in the spdep package
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)

## ----echo = T, eval = T,fig.cap="Figure 4: Plot of observed crime rates for Columbus, OH neighborhoods", warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=60),fig.align="center",fig.pos="h",fig.width=7,fig.height=5.5----

breaks_qt <- classIntervals(c(min(columbus$CRIME) - .00001, columbus$CRIME), n = 7, style = "quantile")

columbus_sf <- mutate(columbus, crime_cat = cut(CRIME, breaks_qt$brks))
ggplot(columbus_sf) + 
    geom_sf(aes(fill=crime_cat)) +
    scale_fill_brewer(palette = "OrRd") + 
    labs(title="Plot of observed \n neighborhood crime rates") +
    theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.y=element_blank(),
          axis.title=element_text(size=25,face="bold"),
          plot.title = element_text(face="bold", size=25, hjust=0.5)) +
    guides(fill=guide_legend("Neighborhood \n crime rate"))

## ----echo = T, eval = T, warning=F,message=F, tidy=TRUE, tidy.opts=list(width.cutoff=58)----

# create neighborhood matrix
columbus.listw <- poly2nb(columbus)
summary(columbus.listw)
W <- nb2mat(columbus.listw, style="B")
Dmat <- diag(apply(W,1,sum))
num.reg <- length(columbus$CRIME)

H <- Dmat - W
H <- (H+t(H))/2
rownames(H) <- NULL
isSymmetric(H)  # check that neighborhood matrix is symmetrix before proceeding

# spectral quantities for use in model selection
H.spectral <- eigen(H, symmetric=TRUE)
Q <- H.spectral$vectors
eigH <- H.spectral$values
phimat <- diag(1/sqrt(eigH[1:(num.reg-1)]))
Sig_phi <- matrix(0,num.reg, num.reg) #initialize
for(i in 1:(num.reg-1)){
  total <- (1/(eigH[i]))*Q[,i]%*%t(Q[,i])
  Sig_phi <- Sig_phi + total
}

# define response and design matrix
Y <- columbus$CRIME
X <- cbind(1, columbus$HOVAL, columbus$INC, columbus$OPEN, columbus$PLUMB, columbus$DISCBD)
b <- (ncol(X)+1)/num.reg  # specify the minimal training size for this example

# perform model selection
columbus.select <- probs.icar(Y=Y,X=X,H=H,
                              H.spectral=H.spectral,
                              Sig_phi=Sig_phi,
                              b=b,verbose=FALSE)

# print the model with highest posterior model probability
columbus.select$probs.mat[which.max(columbus.select$probs.mat[,1]),]

# print vector of posterior inclusion probabilities for each covariate
post.include.cov <- matrix(NA,nrow = 1, ncol=ncol(X)-1)
labels <- c(rep(NA, ncol(X)-1))
for(i in 1:(ncol(X)-1)){labels[i] <- paste("X", i, sep="")}
colnames(post.include.cov) <- labels

for(j in 1:ncol(X)-1){
  post.include.cov[,j] <- sum(columbus.select$probs.mat[grep(paste("X",j,sep=""), columbus.select$probs.mat$'model form'), 1])
}

post.include.cov

