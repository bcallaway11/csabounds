## QoTT

## kernel function
#' @title k
#'
#' @description kernel function
#'
#' @param z evaluate kernel at
#' @param h bandwidth
#' @param type either "gaussian" or "epanechnikov"
#'
#' @return k(z/h)
#'
#' @keywords internal
#' @export
k <- function(z,h=1,type="gaussian") {
    u <- z/h
    if (type=="gaussian") {
        dnorm(u)*(abs(u)<=1) ## gaussian
    } else if (type=="epanechnikov") {
        0.75*(1-u^2)*(abs(u)<=1) ## epanechnikov
    }
}

#' @title G
#'
#' @description vectorized version of logit function
#'
#' @param z compute G(z)
#'
#' @return exp(z)/(1+exp(z))
#'
#' @keywords internal
#' @export
G <- function(z) {
    vapply(z, FUN=function(x) {
        if (exp(x) == Inf) {
            return(1)
        } else {
            exp(x)/(1+exp(x))
        }
    }, 1.0)
}

#' @title dg
#'
#' @description derivative of logit function
#'
#' @param z value to compute g(z)
#'
#' @return exp(z)/((1+exp(z)^2)
#'
#' @keywords internal
#' export
dg <- function(z) {
    exp(z)/((1+exp(z))^2)
}

#' @title wll
#'
#' @description weighted log likelihood function for local logit model
#'  for scalar x only
#'
#' @param bet value of parameters
#' @param y vector of y values
#' @param x vector of x values
#' @param thisx scalar value of x for which evaluating
#' @param h bandwidth
#'
#' @return negative of log likelihood function
#'
#' @keywords internal
#' @export
wll <- function(bet,y,x,thisx,h) {
    idx <- bet[1] + bet[2]*x
    ## matrix version no weights
    ##-sum(y*log(G(X%*%bet)) + log((1-G(X%*%bet)))*(1-y)
    -sum( (y*log(G(idx)) + log((1-G(idx)))*(1-y)) * k(x-thisx,h) )
}

#' @title wscore
#'
#' @description weighted score function
#'
#' @inheritParams wll
#'
#' @return weighted score
#'
#' @keywords internal
#' @export
wscore <- function(bet,y,x,thisx,h) {
    idx <- bet[1] + bet[2]*x
    X <- cbind(1,x)
    kh <- k(x-thisx,h)
    ( ( as.numeric(y*(dg(idx)/G(idx))) - (1-y)*(dg(idx)/(1-G(idx))) )*kh )*X
}


#' @title wgr
#'
#' @description weighted gradient (wscore takes care of weights)
#'
#' @inheritParams wll
#'
#' @return weighted gradient
#'
#' @keywords internal
#' @export
wgr <- function(bet,y,x,thisx,h) {
    colSums(-wscore(bet,y,x,thisx,h))
    ##colSums(-as.numeric(y*(g(X%*%bet)/G(X%*%bet    )))*X + as.numeric((1-y)*(g(X%*%bet)/(1-G(X%*%bet))))*X)
}



#' @title F.Y1
#'
#' @description calculate F(y|ytmin1), the conditional distribution
#'  of treated potential outcomes conditional on ytmin1;
#'  The order of the variables is due to the way that the function
#'  is called later on
#'
#' @param y.seq possible values for y to take
#' @param ytmin1 the value of ytmin1 to condition on
#' @param Y1t vector of outcomes for the treated group in period t
#' @param Y0tmin1 vector of outcomes for the treated group in period t-1
#' @param h optional bandwidth
#' @param method "level" or "rank" determining whether method should
#'  be used conditional on ytmin1 or the rank of ytmin1
#'
#' @return distribution F(y|ytmin1)
#'
#' @examples
#' data(displacements)
#' ytmin1 <- 10
#' Y1t <- subset(displacements, year==2011 & treat==1)$learn
#' Y0tmin1 <- subset(displacements, year==2007 & treat==1)$learn
#' y.seq <- seq(min(c(Y0tmin1,Y1t)), max(c(Y0tmin1,Y1t)), length.out=100)
#' F.Y1(ytmin1, y.seq, Y1t, Y0tmin1)
#' 
#' @export
F.Y1 <- function(ytmin1, y.seq, Y1t, Y0tmin1,  h=NULL, method="level") {
    n <- length(Y1t)
    if (method=="rank") {
        Y0tmin1 <- order(Y0tmin1)/n 
    }
    X <- cbind(1, Y0tmin1-ytmin1)
    if (is.null(h)) {
        h <- 1.06*sd(Y0tmin1)*n^(-1/4) ## check that this is right
    }
    ##h <- h/5
    K <- diag(k(Y0tmin1 - ytmin1,h), n, n)

    Fy1tcondy0tmin1.vals <- vapply(y.seq, FUN=function(y) {
            IY <- 1*(Y1t <= y)
            ## (solve(t(X)%*%K%*%X) %*% t(X)%*%K%*%IY)[1] ##local linear
            ##predict(glm(IY ~ Y0tmin1, family=binomial(link=logit),
            ##            weights=k(Y0tmin1-ytmin1,2)),
            ##        newdata=data.frame(Y0tmin1=ytmin1),
            ##        type="response") ## local logit, I think weights are treated the wrong way here
            o <- optim(c(0,0), wll, gr=wgr, y=IY, x=(Y0tmin1-ytmin1), thisx=0, h=h,##Y0tmin1, thisx=ytmin1, h=h,##thisx=ytmin1, h=h,
                       control=list(maxit=1000, reltol=1e-2),
                       method="BFGS")
            thet <- o$par
            G(thet[1])
    } , 1.0)

    ## rearrangement step
    Fy1tcondy0tmin1.vals <- Fy1tcondy0tmin1.vals[order(Fy1tcondy0tmin1.vals)]
    
    BMisc::makeDist(y.seq, Fy1tcondy0tmin1.vals, TRUE)
}


#' @title F.Y0
#'
#' @description compute F(y|ytmin1) where F is the conditional
#'  distribution of untreated potential outcomes for the treated group
#'  conditional on ytmin1;  This is computed under the copula
#'  stability assumption
#'
#' @inheritParams F.Y1
#' @param Y0tmin1 vector of outcomes for the treated group in period t-1
#' @param Y0tmin2 vector of outcomes for the treated group in period t-2
#' @param Y0tqteobj a qte object for obtaining the counterfactual distribution
#'  of untreated potential outcomes for the treated group in period t
#'
#' @return distribution F(y|ytmin1)
#'
#' @examples
#' data(displacements)
#' ytmin1 <- 10
#' Y1t <- subset(displacements, year==2011 & treat==1)$learn
#' Y0tmin1 <- subset(displacements, year==2007 & treat==1)$learn
#' Y0tmin2 <- subset(displacements, year==2003 & treat==1)$learn
#' y.seq <- seq(min(c(Y0tmin2,Y0tmin1,Y1t)), max(c(Y0tmin2,Y0tmin1,Y1t)), length.out=100)
#' cc <- qte::CiC(learn ~ treat,
#'                t=2011, tmin1=2007, tname="year",
#'                idname="id", panel=TRUE, data=displacements,
#'                probs=seq(.05,.95,.01),se=FALSE)
#' cc$F.treated.tmin2 <- ecdf(subset(displacements, year==2003 & treat==1)$learn)
#' cc$F.treated.tmin1 <- ecdf(subset(displacements, year==2007 & treat==1)$learn)
#' F.Y0(ytmin1, y.seq, Y0tmin1, Y0tmin2, cc)
#' 
#' @export
F.Y0 <- function(ytmin1, y.seq, Y0tmin1, Y0tmin2, Y0tqteobj, h=NULL,
                 method="level") {
    ddid <- Y0tqteobj
    n <- length(Y0tmin1)
    if (method=="rank") {
        Y0tmin2 <- order(Y0tmin2)/n
        xtmin1 <- ytmin1
    } else {
        xtmin1 <- quantile(ddid$F.treated.tmin2,
                           probs=ddid$F.treated.tmin1(ytmin1), type=1)
    }
    X <- cbind(1, Y0tmin2-xtmin1)
    if (is.null(h)) {
        h <- 1.06*sd(Y0tmin2)*n^(-1/4) ## check that this is right
    }
    ##h <- h/5
    K <- diag(k(Y0tmin2-xtmin1, h), n, n)

    Fy0tcondy0tmin1.vals <- vapply(y.seq, FUN=function(y) {
        Z <- 1*(Y0tmin1 <= quantile(ddid$F.treated.tmin1,
                                    probs=ddid$F.treated.t.cf(y),
                                    type=1))
        ## (solve(t(X)%*%K%*%X) %*% t(X)%*%K%*%Z)[1] ## local linear
        ##predict(glm(Z ~ Y0tmin2, family=binomial(link=logit),
        ##                weights=k(Y0tmin2-xtmin1,2)),
        ##            newdata=data.frame(Y0tmin2=xtmin1),
        ##            type="response") ## local logit
        o <- optim(c(0,0), wll, gr=wgr, y=Z, x=(Y0tmin2-xtmin1), thisx=0, h=h,##x=Y0tmin2, thisx=xtmin1, h=h,
                   control=list(maxit=1000, reltol=1e-2),
                   method="BFGS")
            thet <- o$par
            G(thet[1])##G(thet[1] + thet[2]*xtmin1)
    } , 1.0)

    ## rearrangement step
    Fy0tcondy0tmin1.vals <- Fy0tcondy0tmin1.vals[order(Fy0tcondy0tmin1.vals)]
    
    BMisc::makeDist(y.seq, Fy0tcondy0tmin1.vals, TRUE)
}




#' @title l.inner
#'
#' @description for a particular value of y, ytmin1, delt, compute the
#'  max part (before taking sup) in computing the lower bound of DoTT(delt),
#'  see Callaway (2017), Lemma 3
#'
#' @inheritParams l.ytmin1
#' @param y a value to compute the inner part of Lemma 3 for
#'
#' @return scalar that depends on the values of y, ytmin1, and delt
#'
#' @keywords internal
#' @export
l.inner <- function(y, ytmin1, delt, ytmin1.seq,
                    Y1t, Y0tmin1, Y0tmin2, Y0tqteobj, F.y1, F.y0) {
    i <- which(ytmin1.seq==ytmin1)[1]
    max(F.y1[[i]](y) - F.y0[[i]](y-delt),0)
}

#' @title l.ytmin1
#'
#' @description carry out the sup step in computing the lower bound of the
#'  DoTT; everything is still conditional on ytmin1
#'
#' @inheritParams l
#' @param ytmin1 the value of ytmin1 to condition on
#'
#' @return for any value of ytmin1, returns a scaler F^L(delt|ytmin1)
#'
#' @keywords internal
#' @export
l.ytmin1 <- function(ytmin1, delt, y.seq, ytmin1.seq, Y1t, Y0tmin1, Y0tmin2,
                     Y0tqteobj, F.y1, F.y0) {
    max(vapply(y.seq, l.inner, 1.0, ytmin1, delt, ytmin1.seq,
               Y1t, Y0tmin1, Y0tmin2,
               Y0tqteobj, F.y1, F.y0))
}

#' @title l
#'
#' @description Obtains the lower bound on the Distribution of the
#'  Treatment Effect for the Treated (DoTT), DoTT(delt)
#'  under the Copula Stability Assumption.
#'
#' @param delt the value to obtain the DoTT for
#' @param y.seq possible values of y
#' @param Y1t vector of outcomes for the treated group in period t
#' @param Y0tmin1 vector of outcomes for the treated group in period t-1
#' @param Y0tmin2 vector of outcomes for the treated group in period t-2
#' @param h optional bandwidth
#'
#' @return scalar F^L(delt)
#'
#' @keywords internal
#' @export
l <- function(delt, y.seq, ytmin1.seq,
              Y1t, Y0tmin1, Y0tmin2, Y0tqteobj, F.y1, F.y0) {
    mean(vapply(Y0tmin1, l.ytmin1, 1.0, delt, y.seq, ytmin1.seq,
                Y1t, Y0tmin1, Y0tmin2,
                Y0tqteobj, F.y1, F.y0))
}

## get the upper bound on the dte
#' @title l.inner
#'
#' @description for a particular value of y, ytmin1, delt, compute the
#'  min part (before taking sup) in computing the upper bound of DoTT(delt),
#'  see Callaway (2017), Lemma 3
#'
#' @inheritParams u.inner
#' @param y a value to compute the inner part of Lemma 3 for
#'
#' @return scalar that depends on the values of y, ytmin1, and delt
#'
#' @keywords internal
#' @export
u.inner <- function(y, ytmin1, delt, ytmin1.seq,
                    Y1t, Y0tmin1, Y0tmin2, Y0tqteobj,
                    F.y1, F.y0) {
    i <- which(ytmin1.seq==ytmin1)[1]
    1 + min((F.y1[[i]](y) - F.y0[[i]](y-delt)),0)
}

#' @title u.ytmin1
#'
#' @description carry out the 1+inf step in computing the lower bound of the
#'  DoTT; everything is still conditional on ytmin1
#'
#' @inheritParams u
#' @param ytmin1 the value of ytmin1 to condition on
#'
#' @return for any value of ytmin1, returns a scaler F^U(delt|ytmin1)
#'
#' @keywords internal
#' @export
u.ytmin1 <- function(ytmin1, delt, y.seq, ytmin1.seq,
                     Y1t, Y0tmin1, Y0tmin2, Y0tqteobj,
                     F.y1, F.y0) {
    min(vapply(y.seq, u.inner, 1.0, ytmin1, delt, ytmin1.seq,
               Y1t, Y0tmin1, Y0tmin2,
               Y0tqteobj, F.y1, F.y0))
}

#' @title u
#'
#' @description Obtains the upper bound on the Distribution of the
#'  Treatment Effect for the Treated (DoTT), DoTT(delt)
#'  under the Copula Stability Assumption.
#'
#' @inheritParams l
#'
#' @return scalar F^U(delt)
#'
#' @keywords internal
#' @export
u <- function(delt, y.seq, ytmin1.seq,
              Y1t, Y0tmin1, Y0tmin2, Y0tqteobj, F.y1, F.y0) {
    mean(vapply(Y0tmin1, u.ytmin1, 1.0, delt, y.seq, ytmin1.seq,
                Y1t, Y0tmin1, Y0tmin2,
                Y0tqteobj, F.y1, F.y0))
}


#' @title wd.l.inner
#'
#' @description Williamson-Downs bounds inner
#'
#' @inheritParams l
#'
#' @keywords internal
#' @export
wd.l.inner <- function(y, delt, Y1t, Y0tqteobj) {
    max(ecdf(Y1t)(y) - Y0tqteobj$F.treated.t.cf(y-delt),0)
}

#' @title wd.l
#'
#' @description Williamson-Downs lower bound
#'
#' @inheritParams l
#'
#' @keywords internal
#' @export
wd.l <- function(delt, y.seq, Y1t, ddid) {
    max(vapply(y.seq, wd.l.inner, 1.0, delt, Y1t, ddid))
}

#' @title wd.u.inner
#'
#' @description Williamson-Downs bounds inner
#'
#' @inheritParams l
#'
#' @keywords internal
#' @export
wd.u.inner <- function(y, delt, Y1t, ddid) {
    1 + min(ecdf(Y1t)(y) - ddid$F.treated.t.cf(y-delt),0)
}


#' @title wd.u
#'
#' @description Williamson-Downs upper bound
#'
#' @inheritParams l
#'
#' @keywords internal
#' @export
wd.u <- function(delt, y.seq, Y1t, ddid) {
    min(vapply(y.seq, wd.u.inner, 1.0, delt, Y1t, ddid))
}


#' @title csa.bounds
#'
#' @description Compute bounds on the distribution and quantile of the
#'  treatment effect as given in Callaway (2017) under the copula
#'  stability assumption and when a first step estimator of the counterfactual
#'  distribution of untreated potential outcomes for the treated group is
#'  available.
#'
#' @inheritParams F.Y0
#' @inheritParams F.Y1
#' @param formla outcomevar ~ treatmentvar
#' @param data a panel data frame
#' @param t the 3rd period
#' @param tmin1 the 2nd period
#' @param tmin2 the 1st period
#' @param tname the name of the column containing periods
#' @param idname the name of the column containing ids
#' @param delt.seq the possible values to compute bounds on the distribution
#'  of the treatment effect for
#' @param y.seq the possible values for y to take
#' @param F.y0 (optional) pre-computed distribution of counterfactual untreated outcomes for the treated group
#' @param F.y1 (optional) pre-computed distribution of treated outcomes for the treated group
#' @param cl (optional) number of multi-cores to use
#'
#' @examples
#' \dontrun{
#' data(displacements)
#' delt.seq <- seq(-4,4,length.out=50)
#' y.seq <- seq(6.5,13,length.out=50)
#' cc <- qte::CiC(learn ~ treat,
#'                t=2011, tmin1=2007, tname="year",
#'                idname="id", panel=TRUE, data=displacements,
#'                probs=seq(.05,.95,.01),se=FALSE)
#' cc$F.treated.tmin2 <- ecdf(subset(displacements, year==2003 & treat==1)$learn)
#' cc$F.treated.tmin1 <- ecdf(subset(displacements, year==2007 & treat==1)$learn)
#' cb <- csa.bounds(learn ~ treat, 2011, 2007, 2003, "year", "id",
#'         displacements, delt.seq, y.seq, cc,
#'         method="level", cl=1)
#' cb
#' ggCSABounds(cb)
#' }
#' @import stats
#' @importFrom pbapply pblapply
#'
#' @return csaboundsobj
#'
#' @export
csa.bounds <- function(formla, t, tmin1, tmin2, tname, idname,
                       data, delt.seq, y.seq, Y0tqteobj,
                       F.y0=NULL, F.y1=NULL, h=NULL,
                       method=c("level","rank"), cl=1) {

    form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data) 
    colnames(dta) <- c("y","treat")
    data <- cbind.data.frame(dta,data)
    data <- subset(data, treat==1) ## get treated group

    Y1t <- data[data[,tname]==t,]$y
    Y0tmin1 <- data[data[,tname]==tmin1,]$y
    Y0tmin2 <- data[data[,tname]==tmin2,]$y

    
    n <- length(Y1t)
    
    ## if (method=="rank") {
    ##     Y1t <- order(Y1t)/n
    ##     Y0tmin1 <- order(Y0tmin1)/n
    ##     Y0tmin2 <- order(Y0tmin2)/n
    ## }

    ytmin1.seq <- Y0tmin1
    ##ytmin1.seq <- unique(Y0tmin1)  ## if really continuous, this line
    ## doesn't do anything, if not it will add some extra computational
    ## time, but should still work, if not actually continuous, it will
    ## break the method when using rank method, that is why I am commenting
    ytmin1.seq <- ytmin1.seq[order(ytmin1.seq)]

    if (method=="rank") {
        ytmin1.seq <- order(ytmin1.seq)/length(ytmin1.seq)
    }

    print("Step 1 of 4: Calculating conditional distribution of treated potential outcomes...")
    if (is.null(F.y1)) {
        F.y1 <- pbapply::pblapply(ytmin1.seq, F.Y1, y.seq,
                                  Y1t, Y0tmin1, h=h, method=method,
                                  cl=cl)
            ##parallel::mclapply(ytmin1.seq, F.Y1, y.seq,
                ##                   Y1t, Y0tmin1, h=h, method=method,
                ##                   mc.cores=8)
    }

    print("Step 2 of 4: Calculating conditional distribution of untreated potential outcomes...")
    if (is.null(F.y0)) {
        F.y0 <- pbapply::pblapply(ytmin1.seq, F.Y0, y.seq, Y0tmin1,
                                  Y0tmin2, Y0tqteobj, h=h, method=method,
                                  cl=cl)
            ##parallel::mclapply(ytmin1.seq, F.Y0, y.seq, Y0tmin1, Y0tmin2,
                  ##                 Y0tqteobj, h=h, method=method, mc.cores=8)
    }
    ## y.seq <- unique(Y1t, Y0tmin1, Y0tmin2)
    ## y.seq <- y.seq[order(y.seq)]
    ## y.seq <- seq(min(y.seq),max(y.seq),length.out=200)


    if (method=="rank") {
        Y0tmin1r <- order(Y0tmin1)/length(Y0tmin1)
    } else {
        Y0tmin1r <- Y0tmin1
    }

    print("Step 3 of 4: Calculating lower bound")
    l.vec <- pbapply::pblapply(delt.seq, l, y.seq, ytmin1.seq,
                               Y1t, Y0tmin1r, Y0tmin2, Y0tqteobj,
                               F.y1, F.y0, cl=cl)
        ##parallel::mclapply(delt.seq, l, y.seq, ytmin1.seq,
             ##                   Y1t, Y0tmin1r, Y0tmin2, Y0tqteobj, F.y1, F.y0,
             ##                   mc.cores=8) ##TODO
    l.vec <- unlist(l.vec)

    print("Step 4 of 4: Calculating upper bound")
    u.vec <- pbapply::pblapply(delt.seq, u, y.seq, ytmin1.seq,
                               Y1t, Y0tmin1r, Y0tmin2, Y0tqteobj,
                               F.y1, F.y0, cl=cl)
        ##parallel::mclapply(delt.seq, u, y.seq, ytmin1.seq,
             ##                   Y1t, Y0tmin1r, Y0tmin2, Y0tqteobj,
             ##                   F.y1, F.y0, mc.cores=8) ##TODO
    u.vec <- unlist(u.vec)

    F.l <- BMisc::makeDist(delt.seq, l.vec)

    F.u <- BMisc::makeDist(delt.seq, u.vec)

   
    wd.l.vec <- vapply(delt.seq, wd.l, 1.0, y.seq, Y1t, Y0tqteobj)
    wd.u.vec <- vapply(delt.seq, wd.u, 1.0, y.seq, Y1t, Y0tqteobj)

    F.wd.l <- BMisc::makeDist(delt.seq, wd.l.vec)
    F.wd.u <- BMisc::makeDist(delt.seq, wd.u.vec)

    return(list(F.l=F.l, F.u=F.u, F.wd.l=F.wd.l, F.wd.u=F.wd.u))
}

#' @title ggCSABounds
#'
#' @description plot bounds on the quantile of the treatment effect
#'  using ggplot2
#'
#' @param csaboundsobj an object returned from the csa.bounds method
#' @param tau vector of values between 0 and 1 to plot quantiles for
#' @param wdbounds boolean whether or not to also plot Williamson-Downs bounds
#' @param otherdist1 optional ecdf of the distribution of the treatment effect
#'  under cross sectional rank invariance
#' @param otherdist2 optional ecdf of the distribution of the treatment effect
#'  under panel rank invariance
#'
#' @import ggplot2
#'
#' @export
ggCSABounds <- function(csaboundsobj, tau=seq(.05,.95,.05), wdbounds=FALSE,
                        otherdist1=NULL, otherdist2=NULL) {
    tau <- seq(0.05, 0.95, .05)
    c <- csaboundsobj
    
    qu <- quantile(c$F.l, tau, type=1)
    ql <- quantile(c$F.u, tau, type=1)
    qwdu <- quantile(c$F.wd.l, tau, type=1)
    qwdl <- quantile(c$F.wd.u, tau, type=1)

    cmat <- data.frame(tau=tau, qu=qu, ql=ql, group="CSA Bounds")
    cmat2 <- data.frame(tau=tau, qu=qwdu, ql=qwdl, group="WD Bounds")
    if (wdbounds) {
        cmat <- rbind.data.frame(cmat2, cmat)
    }
    if (!is.null(otherdist1)) {
        cmat3 <- data.frame(tau=tau, qu=quantile(otherdist1, tau, type=1),
                            ql=ql, group="CS PPD")
        cmat <- rbind.data.frame(cmat3, cmat)
    }
    if (!is.null(otherdist2)) {
        cmat4 <- data.frame(tau=tau, qu=quantile(otherdist2, tau, type=1),
                            ql=ql, group="Panel PPD")
        cmat <- rbind.data.frame(cmat4, cmat)
    }
    

    p <- ggplot(data=cmat) +
        geom_line(aes(x=tau, y=qu, color=factor(group)), size=1) +
        geom_line(aes(x=tau, y=ql, color=factor(group)), size=1) +
        scale_x_continuous(limits=c(0,1)) + 
        theme_bw() +
        theme(legend.title=element_blank())##%,
              ##legend.position="top",
              ##legend.direction="horizontal")

    p
}
