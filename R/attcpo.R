## kernel function
#' @title k
#'
#' @description kernel function
#'
#' @param z evaluate kernel at value z
#' @param h bandwidth
#'
#' @return gaussian kernel
#' @keywords internal
k <- function(z,h=1) {
    dnorm(z/h)
}

#' @title E.Y1
#'
#' @description a function for computing the conditional expectation
#'  of Y_1t given a value for Y_0tmin1
#'
#' @param ytmin1val scalar value to compute conditional expectation for
#' @param Y1t vector of treated potential outcomes for the treated group
#'  in period t
#' @param Y0tmin1 vector of untreated potential outcomes for the treated
#'  group in period t-1
#' @param h optional bandwidth paramater
#' @param method can be "level" or "rank", whether the conditional expectation
#'  is based on the level of Y0tmin1 or its rank
#'
#' @return the conditional expectation of y1 conditional on y0tmin1
#'
#' @examples
#' data(displacements)
#' ytmin1 <- 10
#' Y1t <- subset(displacements, year==2011 & treat==1)$learn
#' Y0tmin1 <- subset(displacements, year==2007 & treat==1)$learn
#' E.Y1(ytmin1, Y1t, Y0tmin1)
#'
#' @export
E.Y1 <- function(ytmin1val, Y1t, Y0tmin1, h=NULL, method="level") {
    n <- n.treated <- length(Y0tmin1)
    if (method=="rank") {
        Y0tmin1 <- order(Y0tmin1)/n
    }
    y <- ytmin1val
    X <- cbind(1, Y0tmin1-y)
    if (is.null(h)) {
        h <- 1.06*sd(Y0tmin1)*n.treated^(-1/4) ## check that this is right
    }
    K <- diag(k(Y0tmin1 - y,h), n.treated, n.treated)
    (solve(t(X)%*%K%*%X) %*% t(X)%*%K%*%Y1t)[1]
}

#' @title E.Y0
#'
#' @description a function for computing the conditional expectation of
#'   Y_0t given particular value of y_tmin1 under the Copula Stability
#'   Assumption
#'
#' @param ytmin1val the value to compute the conditional expectation for
#' @param Y0tmin1 a vector of untreated potential outcomes for the treated
#'  group in period t-1
#' @param Y0tmin2 a vector of untreated potential outcomes for the treated
#'  group in period t-2
#' @param Y0tqteobj a qte object which should have set F.treated.t.cf
#'  which is the counterfactual distribution of untreated potential outcomes
#'  for the treated group in period t
#' @inheritParams E.Y1
#'
#' @examples
#' data(displacements)
#' ytmin1 <- 10
#' Y0tmin1 <- subset(displacements, year==2007 & treat==1)$learn
#' Y0tmin2 <- subset(displacements, year==2003 & treat==1)$learn
#' cc <- qte::CiC(learn ~ treat,
#'                t=2011, tmin1=2007, tname="year",
#'                idname="id", panel=TRUE, data=displacements,
#'                probs=seq(.05,.95,.01),se=FALSE)
#' cc$F.treated.tmin2 <- ecdf(subset(displacements, year==2003 & treat==1)$learn)
#' cc$F.treated.tmin1 <- ecdf(subset(displacements, year==2007 & treat==1)$learn)
#' E.Y0(ytmin1, Y0tmin1, Y0tmin2, cc)
#'
#' @export
E.Y0 <- function(ytmin1val, Y0tmin1, Y0tmin2, Y0tqteobj,
                 h=NULL, method="level") {
    y <- ytmin1val
    ddid <- Y0tqteobj
    n <- n.treated <- length(Y0tmin1)
    if (method=="rank") {
        Y0tmin2 <- order(Y0tmin2)/n
        xtmin1 <- ytmin1val
    } else {
        xtmin1 <- quantile(ddid$F.treated.tmin2,
                           probs=ddid$F.treated.tmin1(ytmin1val), type=1)
    }

    Z <- quantile(ddid$F.treated.t.cf, probs=ddid$F.treated.tmin1(Y0tmin1),
                  type=1)
    X <- cbind(1, Y0tmin2-xtmin1)
    if (is.null(h)) {
        h <- 1.06*sd(X[,2])*n.treated^(-1/4)##sd(Y0tmin2)*n.treated^(-1/6) ## check that this is right
    }
    K <- diag(k(Y0tmin2-xtmin1,h), n.treated, n.treated)
    (solve(t(X)%*%K%*%X) %*% t(X)%*%K%*%Z)[1]
}

#' @title compute.attcpo
#'
#' @description does the heavy lifting for computing the attcpo
#'
#' @param Y1t treated group outcomes in period t
#' @param Y0tmin1 treated group outcomes in period t-1
#' @param Y0tmin2 treated group outcomes in period t-2
#' @param Y0tqteobj qte object to get counterfactual quantiles from
#' @param h optional bandwidth
#' @param ytmin1seq optional sequence of values to compute attcpo(y') for
#' @param yseqlength optional length of sequence of values to compute
#'  attcpo(y') for; default is 100
#' @param se boolean for whether or not to compute standard errors
#' @param iters number of bootstrap iterations to use to compute standard errors
#' @inheritParams E.Y1
#'
#' @import progress
#'
#' @keywords internal
#' @export
compute.attcpo <- function(Y1t, Y0tmin1, Y0tmin2,
                   Y0tqteobj, h=NULL, ytmin1seq=NULL,
                   yseqlength=NULL, se=TRUE, iters=100, method="level") {
    if (is.null(ytmin1seq)) {
        if ( is.null(yseqlength)) {
            yseqlength <- 100
        }
        y.seq <- seq(quantile(Y0tmin1,.03), quantile(Y0tmin1,.97),
                     length.out=yseqlength)
    } else {
        y.seq <- ytmin1seq
    }

    if (method=="rank") {
        y.seq <- order(y.seq)/length(y.seq)
    }

    print(h)
    
    e.y1 <- vapply(y.seq, E.Y1, 1.0, Y1t=Y1t, Y0tmin1=Y0tmin1, h=h, method=method)
    
    e.y0 <- vapply(y.seq, E.Y0, 1.0, Y0tmin1=Y0tmin1, Y0tmin2=Y0tmin2,
                   Y0tqteobj=Y0tqteobj, h=h, method=method)
    attcpo=e.y1-e.y0

    ##bootstrap
    if (se) {
        e.y1b <- matrix(0, nrow=iters, ncol=length(y.seq))
        e.y0b <- matrix(0, nrow=iters, ncol=length(y.seq))
        attcpob <- rep(0, iters)

        pb <- progress::progress_bar$new(total=iters) 
        
        for (i in 1:iters) {
            pb$tick()
            n <- length(Y1t)
            b <- sample(1:n, n, T)
            Y1tb <- Y1t[b]
            Y0tmin1b <- Y0tmin1[b]
            Y0tmin2b <- Y0tmin2[b]
            ## TODO:  need to estimate counterfactual distribution in each step too
            ##  but this is a problem in current setup, or maybe don't need to
            ##  do it...
            e.y1b[i,] <- vapply(y.seq, E.Y1, 1.0, Y1t=Y1tb, Y0tmin1=Y0tmin1b, h=h, method=method)
            e.y0b[i,] <- vapply(y.seq, E.Y0, 1.0, Y0tmin1=Y0tmin1b, Y0tmin2=Y0tmin2b,
                                Y0tqteobj=Y0tqteobj, h=h, method=method)
        }

        attcpob <- e.y1b - e.y0b

        attcpo.sd <- apply(attcpob, 2, sd)

        V <- var(attcpob)
    } else {
        attcpo.se <- NULL
        V <- NULL
    }
    
    return(ATTCPO.OBJ(y.seq=y.seq, attcpo=e.y1-e.y0, ey1=e.y1, ey0=e.y0,
                attcpo.se=attcpo.sd, V=V, Y0tqteobj=Y0tqteobj))
}


#' @title attcpo
#'
#' @description compute the Average Treatment Effect on the Treated
#'  Conditional on the previous outcome (ATT-CPO)
#'
#' @param formla e.g. y ~ treat
#' @param t the last time period
#' @param tmin1 the middle time period
#' @param tmin2 the first time period
#' @param tname the name of the column containing time periods in the data
#' @param data a data.frame
#' @param idname the name of the column containing an individual identifier over time
#' @param Y0tqteobj a qte object (from the qte package) containing the
#'  the counterfactual distribution of untreated potential outcomes for the
#'  treated group
#' @param h optional bandwidth
#' @param yseq optional sequence of y values, default is to use all unique
#'  yvalues in the data, though this can increase computation time
#' @param yseqlen optional length of y values to use, aids in automatically
#'  generating yseq if desired
#' @param se whether or not to compute standard errors
#' @param iters how many bootstrap iterations to use if computing standard errors; default is 100.
#' @param method should be either "levels" or "rank"; whether to compute the
#'  ATT-CPO using based on the levels of Y0tmin1 or the ranks of Y0tmin1;
#'  "levels" is the default.
#'
#' @examples
#' data(displacements)
#' cc <- qte::CiC(learn ~ treat,
#'                t=2011, tmin1=2007, tname="year",
#'                idname="id", panel=TRUE, data=displacements,
#'                probs=seq(.05,.95,.01),se=FALSE)
#' cc$F.treated.tmin1 <- ecdf(subset(displacements, year==2007 & treat==1)$learn)
#' cc$F.treated.tmin2 <- ecdf(subset(displacements, year==2003 & treat==1)$learn)
#' ac <- attcpo(learn ~ treat, 2011, 2007, 2003, "year", displacements,
#'         "id", cc, method="rank", yseqlen=10)
#' ac
#' ggattcpo(ac)
#'
#' @return att-cpo
#'
#' @export
attcpo <- function(formla, t, tmin1, tmin2, tname, data,
                   idname, Y0tqteobj,
                   h=NULL, yseq=NULL, yseqlen=100, se=TRUE, iters=100, method="level") {
    form <- as.formula(formla)
    dta <- model.frame(terms(form,data=data),data=data) 
    colnames(dta) <- c("y","treat")
    data <- cbind.data.frame(dta,data)
    data <- subset(data, treat==1) ## get treated group

    Y1t <- data[data[,tname]==t,]$y
    Y0tmin1 <- data[data[,tname]==tmin1,]$y
    Y0tmin2 <- data[data[,tname]==tmin2,]$y

    return(compute.attcpo(Y1t=Y1t, Y0tmin1=Y0tmin1, Y0tmin2=Y0tmin2,
                          Y0tqteobj=Y0tqteobj, h=h, ytmin1seq=yseq,
                          yseqlength=yseqlen, se=se, iters=iters, method=method))
}

#' @title ATTCPO.OBJ
#'
#' @description attcpo object
#'
#' @param y.seq possible values for y
#' @param attcpo value of attcpo
#' @param ey1 expectation of treated potential outcomes conditional on y_{t-1}
#' @param ey0 expectation of untreated potential outcomes conditional on y_{t-1}
#' @param attcpo.se standard errors for att.cpo
#' @param V variance matrix
#' @param Y0tqteobj a qte object to get the counterfactual distribution from
#'
#' @return attcpo object
#'
#' @keywords internal
#' @export
ATTCPO.OBJ <- function(y.seq, attcpo, ey1=NULL, ey0=NULL, attcpo.se=NULL, V=NULL,
                   Y0tqteobj=NULL) {
    out <- list(y.seq=y.seq, attcpo=ey1-ey0, ey1=ey1, ey0=ey0,
                attcpo.se=attcpo.se, V=V, Y0tqteobj=Y0tqteobj)
    class(out) <- "ATTCPO.OBJ"
    out
}

#' @title ggattcpo
#'
#' @description plot the ATT-CPO using ggplot2
#'
#' @param attcpoobj an attcpo object
#' @param ylim optional limits of the plot
#'
#' @import ggplot2
#'
#' @export
ggattcpo <- function(attcpoobj, ylim=NULL) {
    jd.attcpo <- attcpoobj
    cc <- attcpoobj$Y0tqteobj

    if (is.null(ylim)) {
        ymaxval <- max(jd.attcpo$attcpo) + max(1.96*jd.attcpo$attcpo.se)
        yminval <- min(jd.attcpo$attcpo) - max(1.96*jd.attcpo$attcpo.se)
    } else {
        ymaxval <- ylim[2]
        yminval <- ylim[1]
    }
    
    cmat <- data.frame(attcpo=jd.attcpo$attcpo, y.seq=jd.attcpo$y.seq, attcpo.se=jd.attcpo$attcpo.se, group="ATT-CPO")
    cutoff <- data.frame(x=c(min(cmat$y.seq),max(cmat$y.seq)),y=cc$ate,cutoff=factor("1"))
    p <- ggplot(data=cmat, aes(attcpo=attcpo, y.seq=y.seq,
                               ymin=attcpo-1.96*attcpo.se,
                               ymax=attcpo+1.96*attcpo.se,
                               group=group)) +
        geom_line(aes(x=y.seq, y=attcpo, color=factor(group)), size=1) +
        geom_line(aes(x=y.seq, y=attcpo+1.96*attcpo.se), data=cmat, linetype="dashed") +
        geom_line(aes(x=y.seq, y=attcpo-1.96*attcpo.se), data=cmat, linetype="dashed") +
        geom_point(aes(x=y.seq, y=attcpo, color=factor(group))) +
        ##geom_line(aes(x=x,y=y,linetype=cutoff), data=cutoff) +
        geom_hline(yintercept=cc$ate, color="blue", size=1, show.legend=TRUE) +
        scale_y_continuous("ATT-CPO", limits=c(yminval, ymaxval)) +
        ##scale_x_continuous(paste(t2,"Earnings",sep=" ")) + 
        theme_bw() +
        theme(legend.title=element_blank())
    ## legend.position="top",
    ## legend.direction="horizontal",
    ## panel.grid.minor=element_blank()
    ## )
    p
}
