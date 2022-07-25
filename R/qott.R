## QoTT

#' @title F.Y1
#'
#' @description calculate F(y|ytmin1), the conditional distribution
#'  of treated potential outcomes conditional on ytmin1.  This function
#'  is typically computed internally in the \code{csabounds} package but
#'  is provided here for convenience.
#'
#' @inheritParams csa.bounds
#' @param data the data.frame.  It should contain colums called "Y1t" and
#'  "Y0tmin1" which correspond to treated potential outcomes in the third
#'  period and untreated potential outcomes in the second period
#' @param retF whether or not to return the distribution function itself;
#'  the default is \code{TRUE}.  If false, it returns a \code{distreg::DR} object
#' @param yvals Sequence of values to compute distributions over
#'  (currently not used as these are determined internally)
#' @param tvals Sequence of values to compute conditional distribution over
#'  (currently not used as these are determined internally)
#'
#' @return distribution F(y|ytmin1)
#'
#' @export
F.Y1 <- function(firststep=c("dr","qr","ll"), xformla, data, yvals, tvals, h=NULL, retF=TRUE) {

    firststep <- firststep[1]
    formla <- BMisc::toformula("Y1t", c("Y0tmin1" , BMisc::rhs.vars(xformla)))
    if (firststep=="dr") {
        dr <- distreg::distreg(formla, data=data, yvals=data$Y1t)
    } else if (firststep=="ll") {
        dr <- distreg::lldistreg(formla, xformla=xformla, data=data, yvals=data$Y1t,
                                       tvals=data$Y0tmin1, h=h)
    } else {
        dr <- quantreg::rq(formla, data=data, tau=seq(.01,.99,.01))
    }

    if (retF) {
        F.y1 <- distreg::Fycondx(dr, data$Y1t, xdf=model.frame(formla, data=data))
        F.y1 ## list with n_1 distribution functions
    } else {
        return(list(dr=dr, yvals=data$Y1t))
    }
}

#' @title F.Y0
#'
#' @description compute F(y|ytmin1) where F is the conditional
#'  distribution of untreated potential outcomes for the treated group
#'  conditional on ytmin1.  This is computed under the copula
#'  stability assumption.  This function
#'  is typically computed internally in the \code{csabounds} package but
#'  is provided here for convenience.  
#'
#' @inheritParams F.Y1
#' @param retZ whether or not to return the distribution of the transformed
#'  random variables due to the copula stability assumption.  This is mainly
#'  used in the numerical bootstrap procedure and, therefore, the default is
#'  \code{FALSE}.
#' @return distribution F(y|ytmin1)
#'
#' @export
F.Y0 <- function(firststep=c("dr","qr","ll"), xformla, data, yvals, tvals, h=NULL, retF=TRUE, retZ=FALSE) {

    ## "adjust" the outcomes under copula stability assumption
    ## here, some things are hard-coded, in particular, using first step quantile regression...
    n <- nrow(data)
    firststep <- firststep[1]
    qformla <- BMisc::toformula("y0t", BMisc::rhs.vars(xformla))
    QR0t <- quantreg::rq(qformla, tau=seq(.01,.99,.01), data=data) ## tau hard-coded...
    QR0tQ <- predict(QR0t, newdata=data, type="Qhat", stepfun=TRUE)
    qformla <- BMisc::toformula("Y0tmin1", BMisc::rhs.vars(xformla))
    QR0tmin1 <- quantreg::rq(qformla, tau=seq(.01,.99,.01), data=data)
    QR0tmin1F <- predict(QR0tmin1, newdata=data, type="Fhat", stepfun=TRUE)
    ## "adjusted" outcome in period t
    Zt <- sapply(1:n, function (i) QR0tQ[[i]](QR0tmin1F[[i]](data$Y0tmin1[i])))
    QR0tmin1Q <- predict(QR0tmin1, newdata=data, type="Qhat", stepfun=TRUE)
    QR0tmin2 <- quantreg::rq(qformla, tau=seq(.01,.99,.01), data=data)
    QR0tmin2F <- predict(QR0tmin2, newdata=data, type="Fhat", stepfun=TRUE)
    ## "adjusted" outcome in period tmin1
    Ztmin1 <- sapply(1:n, function(i) QR0tmin1Q[[i]](QR0tmin2F[[i]](data$Y0tmin2[i])))

    dd <- data.frame(Zt=Zt,Ztmin1=Ztmin1)
    dd <- cbind.data.frame(dd, data)
    formla <- BMisc::toformula("Zt", c("Ztmin1", BMisc::rhs.vars(xformla)))
    if (firststep=="dr") {
        dr <- distreg::distreg(formla, data=dd, yvals=Zt)
    } else if (firststep=="ll") {
        dr <- distreg::lldistreg(formla, data=dd, yvals=y.seq, tvals=dd$Ztmin1)
    } else {
        dr <- quantreg::rq(formla, tau=2)
    }

    if (retF & retZ) {
        F.y0 <- distreg::Fycondx(dr, Zt, xdf=model.frame(formla, data=dd))
        return(list(F.y0=F.y0, xdf=model.frame(formla, data=dd)))
    } else if (retF) {  ## default behavior; everything else is a bit hack way of getting things needed for numerical bootstrap
        F.y0 <- distreg::Fycondx(dr, Zt, xdf=model.frame(formla, data=dd))
        return(F.y0)
    } else {
        return(list(dr=dr,yvals=Zt))
    }
}




#' @title l.inner
#'
#' @description Obtains the lower bound on the Distribution of the
#'  Treatment Effect for the Treated (DoTT), DoTT(delt)
#'  under the Copula Stability Assumption.  These are generic in the
#'  sense that it is going to work without the CSA or with or without covariates.
#'  You just need to provide a vector of distributions -- if there are no
#'  covariates, all of these distributions are going to be the same
#'  and there is not going to be any tightening of the bounds.  If there are
#'  covariates or otherwise conditional distributions (e.g. coming from CSA),
#'  then same computational approach will work, but it will also lead to
#'  (potentially) tighter bounds.
#'
#' @param delt the value to obtain the DoTT for
#' @param y.seq possible values of y
#' @param F.y0 list of distributions of counterfactual untreated outcomes for the treated group
#'  (length of list should be n_1)
#' @param F.y1 list of distributions of treated outcomes for the treated group (length of list should be n_1)
#'
#' @return the value of the lower bound of DoTT(
#'
#' @export
l.inner <- function(delt, y.seq, F.y1, F.y0) {
    n <- length(F.y1)
    condqott <- sapply(1:n, function(i) {
        inner.y <- sapply(y.seq, function(y) {
            max(F.y1[[i]](y) - F.y0[[i]](y-delt), 0)
        })
        max(inner.y) ## this is sup step
    })
    mean(condqott) ## return the average
}

#' @title u.inner
#'
#' @description Obtains the upper bound on the Distribution of the
#'  Treatment Effect for the Treated (DoTT), DoTT(delt)
#'  under the Copula Stability Assumption
#'
#' @inheritParams l.inner
#'
#' @return the value of the upper bound of DoTT(delt)
#' @export
u.inner <- function(delt, y.seq, F.y1, F.y0) {
    n <- length(F.y1)
    condqott <- sapply(1:n, function(i) {
        inner.y <- sapply(y.seq, function(y) {
            1+min(F.y1[[i]](y) - F.y0[[i]](y-delt), 0)
        })
        min(inner.y)
    })
    mean(condqott)
}

#' @title wd.l.inner
#'
#' @description Williamson-Downs bounds inner
#'
#' @inheritParams l.inner
#'
#' @keywords internal
#' @export
wd.l.inner <- function(y, delt, Y1t, Y0tqteobj) {
    max(ecdf(Y1t)(y) - Y0tqteobj$F.treated.t.cf(y-delt),0)
}

#' @title wd.l
#'
#' @description Williamson-Downs lower bound.
#'
#' @inheritParams l.inner
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
#' @inheritParams l.inner
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
#' @inheritParams l.inner
#'
#' @keywords internal
#' @export
wd.u <- function(delt, y.seq, Y1t, ddid) {
    min(vapply(y.seq, wd.u.inner, 1.0, delt, Y1t, ddid))
}



setupTreatedData <- function(data, Y0tqteobj, formla, xformla, idname, tname, t, tmin1, tmin2) {
    form  <- BMisc::toformula(BMisc::lhs.vars(formla),c(BMisc::rhs.vars(formla),BMisc::rhs.vars(xformla)))
    ##form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data)
    colnames(dta) <- c("y","treat",colnames(dta)[-c(1,2)]) ## TODO: update when including covariates
    dta <- cbind.data.frame(dta,data[,c(idname,tname)])

    ## only need treated group for most computations here
    data1 <- subset(dta, treat==1)
    dta1 <- subset(dta, treat==1)

    ## shortened versions of each term
    Y1t <- dta1[dta1[,tname]==t,]$y
    Y0tmin1 <- dta1[dta1[,tname]==tmin1,]$y
    Y0tmin2 <- dta1[dta1[,tname]==tmin2,]$y

    ## this takes covariates from teh first time period
    dta1 <- dta1[dta1[,tname]==tmin2,]
    data1 <- data1[data1[,tname]==tmin2,]
    dta1$Y1t <- Y1t
    dta1$Y0tmin1 <- Y0tmin1
    dta1$Y0tmin2 <- Y0tmin2

    y0t <- Y0tqteobj$y0t
    dta1$y0t <- y0t

    ## if (is.null(xformla)) {
    ##     xformla <- ~1
    ## }

    ## dta1 <- cbind.data.frame(dta1, model.frame(xformla, data=data1))
    dta1
}

#' @title compute.csa.bounds
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
#' @param Y0tqteobj a previously computed first step estimator of the
#'  distribution of counterfactual outcomes for the treated group
#' @param F.y0 (optional) pre-computed distribution of counterfactual untreated outcomes for the treated group
#' @param F.y1 (optional) pre-computed distribution of treated outcomes for the treated group
#' @param cl (optional) number of multi-cores to use
#'
#' @import stats
#' @importFrom pbapply pblapply
#'
#' @return csaboundsobj
#'
#' @export
compute.csa.bounds <- function(formla, t, tmin1, tmin2, tname, idname,
                       data, delt.seq, y.seq, Y0tqteobj, F.y0=NULL, F.y1=NULL,
                       h=NULL,
                       xformla=~1,
                       firststep=c("dr","qr","ll"), cl=1) {


    firststep <- firststep[1]
    ## setup data
    dta1 <- setupTreatedData(data, Y0tqteobj, formla, xformla, idname, tname, t, tmin1, tmin2)

    n <- nrow(dta1)


    ## ytmin1.seq <- Y0tmin1
    ## ##ytmin1.seq <- unique(Y0tmin1)  ## if really continuous, this line
    ## ## doesn't do anything, if not it will add some extra computational
    ## ## time, but should still work, if not actually continuous, it will
    ## ## break the method when using rank method, that is why I am commenting
    ## ytmin1.seq <- ytmin1.seq[order(ytmin1.seq)]
    ##seq(min(Y1t),max(Y1t), length.out=100)

    ## calculate conditional distributions

    ##print("Step 1 of 4: Calculating conditional distribution of treated potential outcomes...")

    if (is.null(F.y1)) {
        ##F.y1 <- F.Y1(firststep, xformla, dta1, yvals=dta1$Y1t, tvals=dta1$Y0tmin1, h=h)
        F.y1 <- F.Y1(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h)
    }

    ##print("Step 2 of 4: Calculating conditional distribution of untreated potential outcomes...")

    if (is.null(F.y0)) {
        ##F.y0 <- F.Y0(firststep, xformla, dta1, yvals=dta1$Y1t, tvals=dta1$Y0tmin1, h=h)
        F.y0 <- F.Y0(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h)
    }

    out <- compute.DoTTBounds(delt.seq, y.seq, F.y1, F.y0)

    ##also compute WD Bounds
    ## compute Williamson-Downs bounds (these are fast, so can compute
    ## them whether or not we report the results)
    wd.l.vec <- vapply(delt.seq, wd.l, 1.0, y.seq, dta1$Y1t, Y0tqteobj)
    wd.u.vec <- vapply(delt.seq, wd.u, 1.0, y.seq, dta1$Y1t, Y0tqteobj)

    F.wd.l <- BMisc::makeDist(delt.seq, wd.l.vec)
    F.wd.u <- BMisc::makeDist(delt.seq, wd.u.vec)

    out$F.wd.l <- F.wd.l
    out$F.wd.u <- F.wd.u

    out
}



#' @title compute.DoTTBounds
#'
#' @description function that does the heavy lifting for computing bounds on the DoTT (and hence
#'  on the QoTT) when conditional distributions of Y1t and Y0t are available
#'
#' @inheritParams compute.csa.bounds
#'
#' @return List with bounds coming from conditional distributions and Williamson-Downs bounds
#' @export
compute.DoTTBounds <- function(delt.seq, y.seq, F.y1, F.y0) {
    ## compute lower bound;
    ##print("Step 3 of 4: Calculating lower bound")

    l.vec <- sapply(delt.seq, l.inner, y.seq, F.y1, F.y0)
    u.vec <- sapply(delt.seq, u.inner, y.seq, F.y1, F.y0)

    ## turn lower and upper bounds into distributions
    F.l <- BMisc::makeDist(delt.seq, l.vec)
    F.u <- BMisc::makeDist(delt.seq, u.vec)

    return(list(F.l=F.l, F.u=F.u))
}

simpleBoot <- function(data) {
    rownums <- sample(1:nrow(data), replace=TRUE)
    data[rownums,]
}




## This is a really ambitious attempt to write this generically for any first
## step estimator; there are some very trick issues here to get this to work --
## a main one is that the first step estimator needs to be computed again for
## each bootstrap draw which is a real headache.  I am giving up on this
## for the moment and just going to write it for the cases under selection
## on observables and under CSA.  These are main cases.  Under CSA, the first
## step estimator is hard coded to be Change-in-changes; but it would be fairly
## minor modification to the code to replace that with something else...
## #' @title DoTTBounds
## #'
## #' @description Compute bounds on the distribution of the treatment effect (DoTT) when the conditional distribution of Y1t and Y0t are available
## #'
## #' @param Y0method This should be a function that can be called with data=data and other parameters
## #'  specified in ... that returns a list of length n_1 of distribution functison for Y0
## #' @param Y1method This should be a function that can be called with data=data and other parameters
## #'  specified in ... that returns a list of length n_1 of distribution functions for Y1
## #' @param simpleBoot A method to bootstrap the data that can be called with data=data and other
## #'  parameters specified in ...; the default is simpleBoot which should work for cross sectional
## #'  data but is going to fail in other cases; in the case with panel data, one can use BMisc::blockBootSample
## #'  but also user-supplied functions can be used here
## #' @inheritParams compute.DoTTBounds
## #'
## #' @return a list containing estimates and standard errors for the DoTT
## #' @export
## DoTTBounds <- function(delt.seq, y.seq, Y0method, Y1method, data, bootmethod=simpleBoot,
##                        se=FALSE, bootiters=100, cl=1, ...) {

##     stop("Not implemented yet...still some bugs...")
##     F.y0 <- Y0method(data=data, ...)
##     F.y1 <- Y1method(data=data, ...)

##     dott <- compute.DoTTBounds(delt.seq, y.seq, F.y1, F.y0)

##     if (se) {
##         n <- nrow(ata)
##         en <- n^(-1/2.1) ## update this...
##         boot.out  <- pbapply::pblapply(1:bootiters, function(b) {
##             boot.data <- bootmethod(data,...)
##             boot.F.y0 <- Y0method(data=boot.data,...)
##             boot.F.y1 <- Y1method(data=boot.data,...)

##             ## ***left off, how to get hat.F...,etc. is not clear...

##             hl.F.y0t <- lapply(1:nrow(boot.dta1), function(i) makeDist(y.seq, hat.F.y0t[[i]](y.seq) - en*sqrt(n)*(hat.F.y0t[[i]](y.seq) - boot.F.y0t[[i]](y.seq)), rearrange=TRUE, force01=TRUE))
##             hl.F.y1t <- lapply(1:nrow(boot.dta1), function(i) makeDist(y.seq, hat.F.y1t[[i]](y.seq) - en*sqrt(n)*(hat.F.y1t[[i]](y.seq) - boot.F.y0t[[i]](y.seq))))

##             csa.boot <- compute.csa.bounds(formla, t, tmin1, tmin2, tname,
##                                            idname, boot.data, delt.seq,
##                                            y.seq, boot.Y0tqteobj, F.y0=hl.F.y0t, F.y1=hl.F.y1t, h,
##                                            xformla, firststep)
##             csa.boot
##         }, cl=cl)


##         ## ***this should work...
##         ## get pointwise confidence intervals for DoTT and QoTT
##         ## CSA lower
##         boot.F.l <- t(sapply(1:bootiters, function(b) {
##             boot.out[[b]]$F.l(delt.seq)
##         }))
##         F.l.upper <- apply(boot.F.l, 2, quantile, (1-alp))
##         F.l.lower <- apply(boot.F.l, 2, quantile, alp)
##         ## CSA upper
##         boot.F.u <- t(sapply(1:bootiters, function(b) {
##             boot.out[[b]]$F.u(delt.seq)
##         }))
##         F.u.upper <- apply(boot.F.u, 2, quantile, (1-alp))
##         F.u.lower <- apply(boot.F.u, 2, quantile, alp)
##         ## WD lower
##         boot.F.wd.l <- t(sapply(1:bootiters, function(b) {
##             boot.out[[b]]$F.wd.l(delt.seq)
##         }))
##         F.wd.l.upper <- apply(boot.F.wd.l, 2, quantile, (1-alp))
##         F.wd.l.lower <- apply(boot.F.wd.l, 2, quantile, alp)
##         ## WD upper
##         boot.F.wd.u <- t(sapply(1:bootiters, function(b) {
##             boot.out[[b]]$F.wd.u(delt.seq)
##         }))
##         F.wd.u.upper <- apply(boot.F.wd.u, 2, quantile, (1-alp))
##         F.wd.u.lower <- apply(boot.F.wd.u, 2, quantile, alp)

##         ## add results to object to return
##         csa.res$F.u.upper <- F.u.upper
##         csa.res$F.u.lower <- F.u.lower
##         csa.res$F.l.upper <- F.l.upper
##         csa.res$F.l.lower <- F.l.lower
##         csa.res$F.wd.u.upper <- F.wd.u.upper
##         csa.res$F.wd.u.lower <- F.wd.u.lower
##         csa.res$F.wd.l.upper <- F.wd.l.upper
##         csa.res$F.wd.l.lower <- F.wd.l.lower
##     }


## }

#' @title csa.bounds
#'
#' @description The main function of the \code{csabounds} package.
#'  It computes bounds on the distribution of the treatment effect
#'  when panel data is available and under the Copula Stability Assumption.
#'  The function can also compute tighter bounds when other covariates are
#'  available.
#' 
#' @param formla A formula of the form: outcome ~ treatment
#' @param t the value for the third time period
#' @param tmin1 the value of the second time period
#' @param tmin2 the value of the first time period
#' @param tname the name of the variable in \code{data} that contains the
#'  time period
#' @param idname the name of the variable in \code{data} that contains the
#'  id variable
#' @param data the name of the \code{data.frame}.  It should be in "long"
#'  format rather than "wide" format.
#' @param delt.seq a vector of values to compute the distribution of the
#'  treatment effect for
#' @param y.seq a vectof of values to compute first-step distributions over
#'  (this is currently not used as it is computed internally)
#' @param Y0tmethod the name of a function to estimate the distribution
#'  of Y0t in a first step; for example \code{qte::panel.qtet},
#'  \code{qte::ddid2}, or \code{qte::CiC}
#' @param h optional bandwidth when using local linear regression
#' @param firststep whether to use distribution regression ("dr"),
#'  quantile regression ("qr"), or local linear distribution regression ("ll")
#'  for the first step estimation of condtional distributions
#' @param xformla a formula for which covariates to use
#' @param se whether or not to compute standard errors (if \code{TRUE},
#'  they are computed using the bootstrap
#' @param bootiters if computing standard errors using the bootstrap, how
#'  many bootstrap iterations to use
#' @param cl if computing standard errors using the bootstrap, how many
#'  cores to use in parallel computation (default is 1)
#' @param alp significance level for confidence intervals
#' @param ... whatever extra arguments need to be passed to Y0tmethod
#'
#' @export
csa.bounds <- function(formla,
                       t,
                       tmin1,
                       tmin2,
                       tname,
                       idname,
                       data,
                       delt.seq,
                       y.seq=NULL,
                       Y0tmethod,
                       h=NULL,
                       xformla=~1,
                       firststep=c("dr","qr","ll"),
                       se=FALSE,
                       bootiters=100,
                       cl=1,
                       alp=.05,...) {

    ## this is not robust at all, probably only works with condcic method,
    ## but you can see how to make it work in more complicated cases...
    Y0tqteobj <- Y0tmethod(formla=formla, xformla=xformla,
                           t=t, tmin1=tmin1, tname=tname,
                           idname=idname, data=data, se=FALSE,...)
    dta1 <- setupTreatedData(data, Y0tqteobj, formla, xformla, idname, tname, t, tmin1, tmin2)
    if (is.null(y.seq)) y.seq <- sort(unique(c(dta1$Y1t,Y0tqteobj$y0t)))
    F.y1 <- F.Y1(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h)
    F.y0 <- F.Y0(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h)
    csa.res <- compute.csa.bounds(formla, t, tmin1, tmin2, tname,
                                  idname, data, delt.seq, y.seq, Y0tqteobj, F.y0=F.y0, F.y1=F.y1,
                                  h=h, xformla=xformla, firststep=firststep)


    if (se) {
        cat("boostrapping standard errors...\n")
        ########################################################
        ##
        ## empirical bootstrap
        ##
        ########################################################

        ## boot.out  <- pbapply::pblapply(1:bootiters, function(b) {
        ##     boot.data <- blockBootSample(data, idname)
        ##     boot.Y0tqteobj <- Y0tmethod(formla=formla, xformla=xformla,
        ##                        t=t, tmin1=tmin1, tname=tname,
        ##                        idname=idname, data=boot.data, se=FALSE,...)
        ##     boot.dta1 <- setupTreatedData(boot.data, boot.Y0tqteobj, formla, xformla, idname, tname, t, tmin1, tmin2)
        ##     boot.F.y0t <- F.Y0(firststep, xformla, boot.dta1, yvals=dta1$Y1t, tvals=dta1$Y0tmin1, h=h)
        ##     boot.F.y1t <- F.Y1(firststep, xformla, boot.dta1, yvals=dta1$Y1t, tvals=dta1$Y0tmin1, h=h)
        ##     csa.boot <- compute.csa.bounds(formla, t, tmin1, tmin2, tname,
        ##                                    idname, boot.data, delt.seq,
        ##                                    y.seq, boot.Y0tqteobj, F.y0=boot.F.y0t, F.y1=boot.F.y1t, h,
        ##                                    xformla, firststep)
        ##     csa.boot
        ## }, cl=cl)


        #########################################################
        ##
        ## numerical bootstrap (hong and li)
        ##
        #########################################################

        ##do some work out of loop...
        n <- nrow(dta1)
        en <- n^(-1/2.5) ## update this...
        dr1 <- F.Y1(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h, retF=FALSE)
        dr0 <- F.Y0(firststep, xformla, dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h, retF=FALSE)
        boot.out  <- pbapply::pblapply(1:bootiters, function(b) {
            boot.data <- BMisc::blockBootSample(data, idname)
            ##boot.data <- data
            boot.Y0tqteobj <- Y0tmethod(formla=formla, xformla=xformla,
                                        t=t, tmin1=tmin1, tname=tname,
                                        idname=idname, data=boot.data, se=FALSE,...)
            boot.dta1 <- setupTreatedData(boot.data, boot.Y0tqteobj, formla, xformla, idname, tname, t, tmin1, tmin2)
            boot.F.y0t.inner <- F.Y0(firststep, xformla, boot.dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h, retZ=TRUE)
            boot.F.y0t <- boot.F.y0t.inner$F.y0
            boot.F.y1t <- F.Y1(firststep, xformla, boot.dta1, yvals=y.seq, tvals=dta1$Y0tmin1, h=h)
            hat.F.y0t <- distreg::Fycondx(dr0$dr, dr0$yvals, xdf=boot.F.y0t.inner$xdf)
            hat.F.y1t <- distreg::Fycondx(dr1$dr, dr0$yvals, xdf=boot.dta1)

            hl.F.y0t <- lapply(1:nrow(boot.dta1), function(i) BMisc::makeDist(y.seq, hat.F.y0t[[i]](y.seq) - en*sqrt(n)*(hat.F.y0t[[i]](y.seq) - boot.F.y0t[[i]](y.seq)), rearrange=TRUE, force01=TRUE, method="linear"))
            hl.F.y1t <- lapply(1:nrow(boot.dta1), function(i) BMisc::makeDist(y.seq, hat.F.y1t[[i]](y.seq) - en*sqrt(n)*(hat.F.y1t[[i]](y.seq) - boot.F.y1t[[i]](y.seq)), rearrange=TRUE, force01=TRUE, method="linear"))


            ## there are slight differences due to interpolations...
            ##sapply(1:nrow(boot.dta1), function(i) all( hl.F.y0t[[i]](y.seq)== F.y0[[i]](y.seq+.5)))
            
            ## first term in hong and li
            csa.boot1 <- compute.csa.bounds(formla, t, tmin1, tmin2, tname,
                                            idname, boot.data, delt.seq,
                                            y.seq, boot.Y0tqteobj, F.y0=hl.F.y0t, F.y1=hl.F.y1t, h,
                                            xformla, firststep)

            hl <- list()
            
            hl$F.l <- (csa.boot1$F.l(delt.seq) - csa.res$F.l(delt.seq))/en
            hl$F.u <- (csa.boot1$F.u(delt.seq) - csa.res$F.u(delt.seq))/en
            hl$F.wd.l <- (csa.boot1$F.wd.l(delt.seq) - csa.res$F.wd.l(delt.seq))/en
            hl$F.wd.u <- (csa.boot1$F.wd.u(delt.seq) - csa.res$F.wd.u(delt.seq))/en
            
            hl
                
        }, cl=cl)

        
        ## get pointwise confidence intervals for DoTT and QoTT
        ## CSA lower
        boot.F.l <- t(sapply(1:bootiters, function(b) {
            boot.out[[b]]$F.l
        }))
        F.l.upper <- BMisc::makeDist(delt.seq, csa.res$F.l(delt.seq) - 
                                apply(boot.F.l, 2, quantile, alp)/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        F.l.lower <- BMisc::makeDist(delt.seq, csa.res$F.l(delt.seq) - 
                                apply(boot.F.l, 2, quantile, (1-alp))/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        ## CSA upper
        boot.F.u <- t(sapply(1:bootiters, function(b) {
          boot.out[[b]]$F.u
        }))
        F.u.upper <- BMisc::makeDist(delt.seq, csa.res$F.u(delt.seq) - 
                                apply(boot.F.u, 2, quantile, alp)/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        F.u.lower <- BMisc::makeDist(delt.seq, csa.res$F.u(delt.seq) - 
                                apply(boot.F.u, 2, quantile, (1-alp))/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        ## WD lower
        boot.F.wd.l <- t(sapply(1:bootiters, function(b) {
          boot.out[[b]]$F.wd.l
        }))
        F.wd.l.upper <- BMisc::makeDist(delt.seq, csa.res$F.wd.l(delt.seq) -
                                   apply(boot.F.wd.l, 2, quantile, alp)/sqrt(n),
                                 rearrange=TRUE, force01=TRUE)
        F.wd.l.lower <- BMisc::makeDist(delt.seq, csa.res$F.wd.l(delt.seq) - 
                                   apply(boot.F.wd.l, 2, quantile, (1-alp))/sqrt(n),
                                 rearrange=TRUE,force01=TRUE)
        ## WD upper
        boot.F.wd.u <- t(sapply(1:bootiters, function(b) {
          boot.out[[b]]$F.wd.u
        }))
        F.wd.u.upper <- BMisc::makeDist(delt.seq, csa.res$F.wd.u(delt.seq) - 
                                   apply(boot.F.wd.u, 2, quantile, alp)/sqrt(n),
                                 rearrange=TRUE, force01=TRUE)
        F.wd.u.lower <- BMisc::makeDist(delt.seq, csa.res$F.wd.u(delt.seq) - 
                                   apply(boot.F.wd.u, 2, quantile, (1-alp))/sqrt(n),
                                 rearrange=TRUE, force01=TRUE)
        
        
        ## add results to object to return
        csa.res$F.u.upper <- F.u.upper
        csa.res$F.u.lower <- F.u.lower
        csa.res$F.l.upper <- F.l.upper
        csa.res$F.l.lower <- F.l.lower
        csa.res$F.wd.u.upper <- F.wd.u.upper
        csa.res$F.wd.u.lower <- F.wd.u.lower
        csa.res$F.wd.l.upper <- F.wd.l.upper
        csa.res$F.wd.l.lower <- F.wd.l.lower
    }

    return(csa.res)


}


#' @title cia.bounds
#'
#' @description Bounds on the distribution of the treatment effect and on
#'  the quantile of the treatment effect under a conditional independence
#'  assumption (cia).  It takes in a data.frame, that should indicate whether
#'  or not an individual is treated; separates individuals into a treated
#'  and untreated group; runs distribution or quantile regression to
#'  estimate the conditional distributions; then computes bounds on the
#'  DoTT or QoTT.
#'
#' @inheritParams csa.bounds
#' @param formla y ~ d
#' @param xformla ~ x1 + x2
#' @param alp significance level for confidence intervals
#' @param link optional argument to pass to \code{distreg} for which link
#'  function to use when running distribution regressions
#' @param ... whatever extra arguments need to be passed to Y0tmethod
#'
#' @export
cia.bounds <- function(formla, xformla=~1, data, delt.seq, y.seq=NULL, 
                       firststep=c("dr","qr","ll"),
                       link="logit",
                       se=FALSE, bootiters=100,
                       cl=1,alp=.05,...) {

    ## setup data
    form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data)
    colnames(dta) <- c("y","treat") ## TODO: update when including covariates
    dta <- cbind.data.frame(dta, model.frame(xformla, data=data))
    dta1 <- subset(dta, treat==1)
    dta0 <- subset(dta, treat==0)                            
    if (is.null(y.seq)) y.seq <- sort(unique(dta1$y)) ## this is going to get "on the treated" results
    
    
    ## compute conditional distributions
    dformla <- BMisc::toformula("y", BMisc::rhs.vars(xformla))
    xdf <- dta1 ## (implies we get qott, etc. rather than qote)
    if (firststep=="dr") {
        dr1 <- distreg::distreg(dformla, data=dta1, yvals=y.seq,
                                      link=link)
        dr0 <- distreg::distreg(dformla, data=dta0, yvals=y.seq,
                                      link=link)
        
    } else {
        dr1 <- quantreg::rq(dformla, data=dta1, tau=seq(.01,.99,.01))
        dr0 <- quantreg::rq(dformla, data=dta0, tau=seq(.01,.99,.01))
    }
    F.y1 <- distreg::Fycondx(dr1, y.seq, xdf=xdf)
    F.y0 <- distreg::Fycondx(dr0, y.seq, xdf=xdf)

    ## compute estimates of bounds
    cia.res <- compute.DoTTBounds(delt.seq, y.seq, F.y1, F.y0)


    if (se) {
        cat("boostrapping standard errors...\n")
        #########################################################
        ##
        ## numerical bootstrap (hong and li)
        ##
        #########################################################

        ##do some work out of loop...
        n <- nrow(dta1)
        en <- n^(-1/2.5) ## update this...
        boot.out  <- pbapply::pblapply(1:bootiters, function(b) {
            ## setup bootstrap data
            boot.data <- simpleBoot(dta)
            boot.dta1 <- subset(boot.data, treat==1)
            boot.dta0 <- subset(boot.data, treat==0)                            

            ## roughly follow same steps as before, combining with Hong and Li
            ## compute conditional distributions
            boot.xdf <- boot.dta1 ## (implies we get qott, etc. rather than qote)
            if (firststep=="dr") {
                boot.dr1 <- distreg::distreg(dformla, data=boot.dta1, yvals=y.seq, link=link)
                boot.dr0 <- distreg::distreg(dformla, data=boot.dta0, yvals=y.seq, link=link)
                
            } else {
                boot.dr1 <- quantreg::rq(dformla, data=boot.dta1, tau=seq(.01,.99,.01))
                boot.dr0 <- quantreg::rq(dformla, data=boot.dta0, tau=seq(.01,.99,.01))
            }
            boot.F.y1 <- distreg::Fycondx(boot.dr1, y.seq, xdf=boot.xdf)
            boot.F.y0 <- distreg::Fycondx(boot.dr0, y.seq, xdf=boot.xdf)

            hat.F.y0 <- distreg::Fycondx(dr0, y.seq, xdf=boot.xdf)
            hat.F.y1 <- distreg::Fycondx(dr1, y.seq, xdf=boot.xdf)

            hl.F.y0 <- lapply(1:nrow(boot.dta1), function(i) BMisc::makeDist(y.seq, hat.F.y0[[i]](y.seq) - en*sqrt(n)*(hat.F.y0[[i]](y.seq) - boot.F.y0[[i]](y.seq)), rearrange=TRUE, force01=TRUE))
            hl.F.y1 <- lapply(1:nrow(boot.dta1), function(i) BMisc::makeDist(y.seq, hat.F.y1[[i]](y.seq) - en*sqrt(n)*(hat.F.y1[[i]](y.seq) - boot.F.y1[[i]](y.seq)), rearrange=TRUE, force01=TRUE))

            cia.boot <- compute.DoTTBounds(delt.seq, y.seq, F.y1=hl.F.y1, F.y0=hl.F.y0)
            
            hl <- list()
            ## note that these are not actually distributions; just use these to calculate confidence intervals later
            ## hl$F.l <- (csa.boot1$F.l(delt.seq) - csa.boot2$F.l(delt.seq))/en
            ## hl$F.u <- (csa.boot1$F.u(delt.seq) - csa.boot2$F.u(delt.seq))/en
            ## hl$F.wd.l <- (csa.boot1$F.wd.l(delt.seq) - csa.boot2$F.wd.l(delt.seq))/en
            ## hl$F.wd.u <- (csa.boot1$F.wd.u(delt.seq) - csa.boot2$F.wd.u(delt.seq))/en
            hl$F.l <- (cia.boot$F.l(delt.seq) - cia.res$F.l(delt.seq))/en
            hl$F.u <- (cia.boot$F.u(delt.seq) - cia.res$F.u(delt.seq))/en
            
            hl
            
        }, cl=cl)
        
        boot.F.l <- t(sapply(1:bootiters, function(b) {
          boot.out[[b]]$F.l
        }))
        F.l.upper <- BMisc::makeDist(delt.seq, cia.res$F.l(delt.seq) - 
                                apply(boot.F.l, 2, quantile, alp)/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        F.l.lower <- BMisc::makeDist(delt.seq, cia.res$F.l(delt.seq) - 
                                apply(boot.F.l, 2, quantile, (1-alp))/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        ## CSA upper
        boot.F.u <- t(sapply(1:bootiters, function(b) {
          boot.out[[b]]$F.u
        }))
        F.u.upper <- BMisc::makeDist(delt.seq, cia.res$F.u(delt.seq) - 
                                apply(boot.F.u, 2, quantile, alp)/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        F.u.lower <- BMisc::makeDist(delt.seq, cia.res$F.u(delt.seq) - 
                                apply(boot.F.u, 2, quantile, (1-alp))/sqrt(n),
                              rearrange = TRUE, force01 = TRUE)
        # ## WD lower
        # boot.F.wd.l <- t(sapply(1:bootiters, function(b) {
        #   boot.out[[b]]$F.wd.l
        # }))
        # F.wd.l.upper <- BMisc::makeDist(delt.seq, cia.res$F.wd.l(delt.seq) -
        #                            apply(boot.F.wd.l, 2, quantile, alp)/sqrt(n),
        #                          rearrange=TRUE, force01=TRUE)
        # F.wd.l.lower <- BMisc::makeDist(delt.seq, cia.res$F.wd.l(delt.seq) - 
        #                            apply(boot.F.wd.l, 2, quantile, (1-alp))/sqrt(n),
        #                          rearrange=TRUE,force01=TRUE)
        # ## WD upper
        # boot.F.wd.u <- t(sapply(1:bootiters, function(b) {
        #   boot.out[[b]]$F.wd.u
        # }))
        # F.wd.u.upper <- BMisc::makeDist(delt.seq, cia.res$F.wd.u(delt.seq) - 
        #                            apply(boot.F.wd.u, 2, quantile, alp)/sqrt(n),
        #                          rearrange=TRUE, force01=TRUE)
        # F.wd.u.lower <- BMisc::makeDist(delt.seq, cia.res$F.wd.u(delt.seq) - 
        #                            apply(boot.F.wd.u, 2, quantile, (1-alp))/sqrt(n),
        #                          rearrange=TRUE, force01=TRUE)
        
        ## add results to object to return
        cia.res$F.u.upper <- F.u.upper
        cia.res$F.u.lower <- F.u.lower
        cia.res$F.l.upper <- F.l.upper
        cia.res$F.l.lower <- F.l.lower
        # cia.res$F.wd.u.upper <- F.wd.u.upper
        # cia.res$F.wd.u.lower <- F.wd.u.lower
        # cia.res$F.wd.l.upper <- F.wd.l.upper
        # cia.res$F.wd.l.lower <- F.wd.l.lower
    }

    return(cia.res)

}



#' @title ggCSABounds
#'
#' @description plot bounds on the quantile of the treatment effect
#'  using ggplot2
#'
#' @param csaboundsobj an object returned from the csa.bounds method
#' @param tau vector of values between 0 and 1 to plot quantiles for
#' @param csalabel label for bounds under copula stability assumption
#' @param wdbounds boolean whether or not to also plot Williamson-Downs bounds
#' @param wdlabel label for worst case bounds
#' @param otherdist1 optional ecdf of the distribution of the treatment effect
#'  under cross sectional rank invariance
#' @param otherdist1label label for other distribution 1
#' @param otherdist2 optional ecdf of the distribution of the treatment effect
#'  under panel rank invariance
#' @param otherdist2label label for other distribution 2
#' @param plotSE whether or not to plot standard errors (default is \code{TRUE})
#' @param colors optional vector of colors to use in place of \code{ggplot}'s
#'  default scheme
#' @param legend whether or not to include a legend (default is \code{TRUE})
#'
#' @import ggplot2
#'
#' @export
ggCSABounds <- function(csaboundsobj, tau=seq(.05,.95,.05),
                        csalabel="CSA Bounds",
                        wdbounds=FALSE,
                        wdlabel="WD Bounds",
                        otherdist1=NULL,
                        otherdist1label="CS RI",
                        otherdist2=NULL,
                        otherdist2label="Panel RI", plotSE=TRUE,
                        colors=NULL, legend=TRUE) {
    tau <- seq(0.05, 0.95, .05)
    c <- csaboundsobj

    qu <- quantile(c$F.l, tau, type=1)
    ql <- quantile(c$F.u, tau, type=1)
    qwdu <- quantile(c$F.wd.l, tau, type=1)
    qwdl <- quantile(c$F.wd.u, tau, type=1)

    cmat <- data.frame(tau=tau, qu=qu, ql=ql, group=csalabel)
    cmat2 <- data.frame(tau=tau, qu=qwdu, ql=qwdl, group=wdlabel)
    if (wdbounds) {
        cmat <- rbind.data.frame(cmat, cmat2)
    }
    if (!is.null(otherdist1)) {
        cmat3 <- data.frame(tau=tau, qu=quantile(otherdist1, tau, type=1),
                            ql=ql, group=otherdist1label)
        cmat <- rbind.data.frame(cmat, cmat3)
    }
    if (!is.null(otherdist2)) {
        cmat4 <- data.frame(tau=tau, qu=quantile(otherdist2, tau, type=1),
                            ql=ql, group=otherdist2label)
        cmat <- rbind.data.frame(cmat, cmat4)
    }


    p <- ggplot(data=cmat) +
        geom_line(aes(x=tau, y=qu, color=factor(group)), size=1) +
        geom_line(aes(x=tau, y=ql, color=factor(group)), size=1) +
        scale_x_continuous(limits=c(0,1)) +
        theme_bw() +
        ylab("QoTT") + 
        theme(legend.title=element_blank())##%,
              ##legend.position="top",
    ##legend.direction="horizontal")

    if (!(is.null(colors))) {
        p <- p + scale_color_manual(values=colors)
    }

    if (!legend) {
        p <- p + theme(legend.position="none")
    }
    
    ## add confidence intervals, if they are provided
    if (plotSE) {
        if (!is.null(c$F.l.upper)) {
            dmat <- data.frame(tau=tau, uci=quantile(c$F.l.lower, probs=tau, type=1), lci=quantile(c$F.u.upper, probs=tau, type=1))
            p <- p + geom_line(data=dmat, mapping=aes(x=tau, y=uci), size=1, linetype=2)
            p <- p + geom_line(data=dmat, mapping=aes(x=tau, y=lci), size=1, linetype=2)
        }
    }

    p
}
