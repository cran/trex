#$Author: sinnwell $
#$Date: 2010/04/05 14:15:13 $
#$Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/optimalDesign.q,v 1.1 2010/04/05 14:15:13 sinnwell Exp $
#$Locker:  $
#$Log: optimalDesign.q,v $
#Revision 1.1  2010/04/05 14:15:13  sinnwell
#Initial revision
#
#Revision 1.1  2010/04/02 21:38:06  sinnwell
#Initial revision
#


optimalDesign <- function(probCarrier, or, alpha, power, maxThreshold=50){
  # input probCarrier in population, or for risk of carriers,
  # alpha, power (unconditional power), and maxThreshold for stage-1.
  # Here, maxThreshold is the maximum allowable threshold over all designs.

  # the returned "designs" data.frame contains the following:
  # Ntotal = total number of cases and controls (assuming Ntotal/2 = Ncases = Ncontrols)
  # Power.cond = conditional power, conditional on continuing to stage-2
  # Ncase1 =  number of cases for stage-1
  # Threshold = continuation threshold for stage-1 (no. carriers in stage-1 >= threshold)
  # ProbStopNull = probability of stopping at stage-1 under the null (OR = 1)
  # ProbStopAlt = probabilit of stopping at stage-1 under alterntative (OR = input or)
  # ESN.Null = expected sample size under the null
  # ESN.Alt = expected sample size under the alternative specified or
  # Optimal = TRUE for design with minimum ESN.Null

  pow.uncond <- power

  pow.seq <- seq(from=pow.uncond + .01, to = .99, by = .01)
  choices <- NULL

  for(pow in pow.seq){
   n  <- findN(probCarrier, or, alpha, pow.conditional = pow)
   dat1 <- design.stage1(pCarrierNull = probCarrier, or=or, maxCases= n/2, maxThreshold=maxThreshold)
   save <- possibleDesign(pCarrier=probCarrier, or=or, nTotal= n, pow.min=pow.uncond, dat=dat1)
   choices <- rbind(choices, save)
  }

  designs <- as.data.frame(choices, row.names=1:nrow(choices))
  
  opt <- designs$esn.null == min(designs$esn.null)
  designs <- cbind(designs, optimal=opt)
  names(designs) <- c("ProbCarrier","OR","Ntotal","Power.cond","Power","Ncase1",
     "Threshold","ProbStopNull","ProbStopAlt","ESN.Null","ESN.Alt","Optimal")

  # remove constant terms from data frame
  designs$ProbCarrier=NULL
  designs$OR = NULL
  designs$Power = NULL

  return(list(probCarrier=probCarrier, or=or, alpha=alpha, power=power, designs))
}


############### utility funcs for optimalDesign ###################################


dbin <- function(x, n, p){
  # binomial density for vector x, integer n, probability p

  lnprob <-	lgamma(n+1) - lgamma(x+1) - lgamma(n-x+1) + x*log(p) + (n-x)*log(1-p)
  return(exp(lnprob))
}


probStop <- function(thresh, ncase1, probCarrier) {
  # prob stop at stage-1, for threshold (thresh), number of cases at stage-1
  # (ncase1) and probabilty of carrying any variants (probCarrier)

   pstop <- sum(dbin(0:(thresh-1), n = ncase1, p= probCarrier))
   return(pstop)
}


powerCaseControl <- function(po,or,alpha,nTotal){

 # Power for an unmatched case-control study with 
 # 1 control per case.

 # po = exposure rate among controls
 # or = relative risk (odds ratio)
 # alpha = 1-sided alpha error. Use  alpha/2 for 2-sided test
 # nTotal = total no. cases + no. controls
 # Power from eq (6.9), p 151 of Schlesselman

 qo <-  1 - po
 p1 <- po*or/(1 + po*(or-1))
 q1 <-  1 - p1
 pbar <- (p1+po)/(2)
 qbar <-  1 - pbar
 n <- nTotal / 2
 za <- qnorm(1 - alpha)
 zb <-  ( sqrt(n*(p1-po)^2) - za * sqrt((2)*pbar*qbar) )/
     sqrt(p1*q1 + po*qo)
 pow <- pnorm(zb)
 return(pow)
}


opt.min.esn <- function(dat, nTotal, pow.conditional, pow.min){

  # compute ESN under null and alt, and return design with 
  # sufficient desired power and min ESN.Null

  esn.null <- dat$design$ncase1*dat$design$pstopNull +
		 nTotal*(1 - dat$design$pstopNull)

  pow.unconditional <- (1-dat$design$pstopAlt) * pow.conditional
 
  if(max(pow.unconditional) < pow.min){
	pow.min <- max(pow.unconditional)
  }

  # choose design with desired power and min ESN.Null

  zed <- pow.unconditional >= pow.min
  minN <- min(esn.null[zed])
  zed <- zed & esn.null == minN

  esn.alt <- dat$design$ncase1*dat$design$pstopAlt +
		 nTotal*(1 - dat$design$pstopAlt)
  esn.alt <- esn.alt[zed]

  return(unlist(c(pow.uncond=pow.min, dat$design[zed,],esn.null=minN, esn.alt= esn.alt)))
}


design.stage1 <- function(pCarrierNull, or, maxCases, maxThreshold){

 pCarrierAlt <- pCarrierNull * or /(1 + pCarrierNull *(or-1))

 pstopNull <- pstopAlt <- NULL
 thresh <- ncase1 <- NULL

 increment <- 1
 case1Vec <- seq(from=10, to=maxCases, by=increment)
 for(case1 in case1Vec){

  threshVec <- 1:min(case1, maxThreshold)
 
  for(t in threshVec){

    probStopNull <- probStop(thresh=t, ncase1=case1, probCarrier=pCarrierNull)
    probStopAlt <-  probStop(thresh=t, ncase1=case1, probCarrier=pCarrierAlt)

    pstopNull <- c(pstopNull,  probStopNull)
    pstopAlt  <- c(pstopAlt,   probStopAlt)
    thresh <- c(thresh, t)
    ncase1 <- c(ncase1, case1)
    # increasting threshold increases prob stop under null and alt, so if both
    # prob's for stop = 1, no need to increase threshold further
    if(probStopNull > 0.999999999 & probStopAlt > 0.999999999) break
   }
 }

 design1 <- data.frame(ncase1, thresh, pstopNull, pstopAlt)

 return(list(design=design1, pCarrierNull = pCarrierNull, or=or, maxN=maxCases))
}


possibleDesign <- function(pCarrier, or, nTotal, pow.min, alpha=.05, dat=NULL){
  # compute a possible design that has desired power

  pow.conditional <- powerCaseControl(pCarrier, or=or, alpha=alpha, nTotal=nTotal)

  if(is.null(dat)){
    dat <- design.stage1(pCarrierNull = pCarrier, or=or, maxCases=nTotal/2)
  }
  tmp <- opt.min.esn(dat, nTotal, pow.conditional, pow.min=pow.min)

  tbl.opt <- list(pCarrier=pCarrier, or=or, nTotal=nTotal, 
	pow.cond=pow.conditional, tmp)
  return(unlist(tbl.opt))
}

bisect2 <- function(fn, pars, xL, xU)
{

  # function to compute the minimum of the input function, fn, where
  # fn is a defined function with 2  input parameters: pars and  x.
  # pars can be a list of other needed fixed parameters, and x is the
  # iterating parameter. 
  # The returned solution x* lies between xL and xU. This alogorithm is
  # based on the bisect method.

  niter <- 50
  eps <- 1e-05
  f <- fn(pars,xL)
  fmid <- fn(pars,xU)
  if(f * fmid >= 0)
  {
    stop(paste("bisect2 root not bracketed in ", xL, xU))
  }
  else{
   rtb <- xL
   dx <- xU - xL
   for(i in 1:niter) {
     iters <- i
     dx <- dx/2
     xmid <- rtb + dx
     fmid <- fn(pars,xmid)
     if(fmid <= 0) {
	 rtb <- xmid
     }
     if(abs(dx) < eps || fmid == 0)
       break
  }
  rtb
  }
 }

powBisect2.func <- function(pars, n){
  pCarrier=pars$pCarrier
  or <- pars$or
  alpha <- pars$alpha
  pow.target <- pars$pow.target
  pow <- powerCaseControl(pCarrier, or=or, alpha=alpha, nTotal=n)
  return(pow - pow.target)
}

  
findN <- function(pCarrier, or, alpha, pow.conditional){
  pars <- list(pCarrier=pCarrier, or=or, alpha=alpha, pow.target=pow.conditional)
  n <- ceiling(bisect2(powBisect2.func, pars, xL=100, xU=10000))
  if( (n - floor(n/2)*2) != 0) n <- n+1
  return(n)
}


######### examples presented in Table 4, Schaid and Sinnwell (Human Genetics)

#example1 <- optimalDesign(probCarrier=.01, or=5, alpha=.05, power=.9)

#example2 <- optimalDesign(probCarrier=.05, or=5, alpha=.05, power=.9)

