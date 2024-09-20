## Code for RMA-Mean and RMA-75 as discussed in McGee and Chen
##
## To use the following code, copy and paste it into a plain text file.
## Open an R session, source the file into R and run once.
## Then type:
##
## library(affy)
## bgcorrect.methods=c(bgcorrect.methods,"rma001","rma002")
## These command load the affy package, and add RMA-Mean and
## RMA-75 to the list of background correction methods.
##
## b = ReadAffy()
## This obtains an AffyBatch object from your CEL files
##
## a = expresso(b, bgcorrect.method="rma001 (or rma002)",
## normalize.method="quantiles", pmcorrect.method="pmonly",
## summary.method="medianpolish")
##
## To read output into a CSV file on original scale (makes retrieval
## easier)
## aa=exprs(a)
## exprs(a)=2^aa
## aaa=exprs(a)
## write.table(data.frame(aaa,check.names=FALSE),
## "filename.csv",sep=",",col.names=NA,quote=FALSE)
## RMA-Mean
bg.correct.rma001=function (object, ...)
{
  bg.parameters001=function(pm,n.pts = 2^14)
  {
    ## mu-correction function
    mu.est.correct=function(m,s,a){
      f <- function(x) (dnorm(x-s*a)-s*a*(pnorm(x-s*a)+pnorm(m/s+s*a)-
                                            1))
      t=uniroot(f, c(-5, 10), tol = 1e-12)$root
      t=m-s*t
      return(t)
    }
    ### getting mode function
    max.density <- function(x, n.pts) {
      aux <- density(x, kernel = "epanechnikov", n = n.pts,
                     na.rm = TRUE)
      aux$x[order(-aux$y)[1]]
    }
    pmbg <- max.density(pm, n.pts)
    bg.data <- pm[pm < pmbg]
    pmbg <- max.density(bg.data, n.pts)
    mubg <- pmbg ## the mode
    bg.data <- pm[pm < pmbg]
    bg.data <- bg.data - pmbg
    bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2) ##
    #estimate sigma
    sig.data<-pm[pm > pmbg]
    sig.data <- sig.data - pmbg
    alpha1=mean(sig.data) ## by mean
    ## mode-correction
    mu1=mu.est.correct(m=mubg,s=bgsd,a=1/alpha1)
    mu1=(mu1+mubg)/2 ## take ave
    bg.data1<- pm[pm < mu1]
    bg.data1 <- bg.data1 - mu1
    bgsd1 <- sqrt(sum(bg.data1^2)/(length(bg.data1) - 1)) * sqrt(2)
    sig.data1<-pm[pm > mu1]
    sig.data1<-sig.data1-mu1
    alpha1<-mean(sig.data1)
    list(alpha = 1/alpha1, mu = mu1, sigma = bgsd1) ## be careful here
    "alpha"
  }
  ## bg.correct:
  bg.adjust001=function (pm, n.pts = 2^14, ...)
  {
    param <- bg.parameters001(pm, n.pts)
    b <- param$sigma
    pm <- pm - param$mu - param$alpha * b^2
    pm + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((pm/b)^2)))/pnorm(pm/b)
  }
  pm(object) <- apply(pm(object), 2, bg.adjust001)
  return(object)
}
## RMA-75
bg.correct.rma002=function (object, ...)
{
  bg.parameters002=function(pm,n.pts = 2^14)
  {
    ## mu-correction function
    mu.est.correct=function(m,s,a){
      f <- function(x) (dnorm(x-s*a)-s*a*(pnorm(x-s*a)+pnorm(m/s+s*a)-
                                            1))
      t=uniroot(f, c(-5, 10), tol = 1e-12)$root
      t=m-s*t
      return(t)
    }
    ### getting mode function
    max.density <- function(x, n.pts) {
      aux <- density(x, kernel = "epanechnikov", n = n.pts,
                     na.rm = TRUE)
      aux$x[order(-aux$y)[1]]
    }
    pmbg <- max.density(pm, n.pts)
    bg.data <- pm[pm < pmbg]
    pmbg <- max.density(bg.data, n.pts)
    mubg <- pmbg ## the mode
    bg.data <- pm[pm < pmbg]
    bg.data <- bg.data - pmbg
    bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2) ##
    #estimate sigma
    sig.data<-pm[pm > pmbg]
    sig.data <- sig.data - pmbg
    q75=0.75
    alpha3=-(quantile(pm,q75)-pmbg)/log(1-q75) ##quantile 75%
    estimation
    
    ## mode-correction
    mu3=mu.est.correct(m=mubg,s=bgsd,a=1/alpha3)
    mu3=(mu3+mubg)/2 ## take ave
    bg.data3<- pm[pm < mu3]
    bg.data3 <- bg.data3 - mu3
    bgsd3 <- sqrt(sum(bg.data3^2)/(length(bg.data3) - 1)) * sqrt(2)
    sig.data3<-pm[pm > mu3]
    alpha3=-(quantile(pm,q75)-mu3)/log(1-q75)
    list(alpha = 1/alpha3, mu = mu3, sigma = bgsd3)
  }
  ## bg.correct:
  bg.adjust002=function (pm, n.pts = 2^14, ...)
  {
    param <- bg.parameters002(pm, n.pts)
    b <- param$sigma
    pm <- pm - param$mu - param$alpha * b^2
    pm + b * ((1/sqrt(2 * pi)) * exp((-1/2) * ((pm/b)^2)))/pnorm(pm/b)
  }
  pm(object) <- apply(pm(object), 2, bg.adjust002)
  return(object)
}