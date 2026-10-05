##########
# Internal method
##########
.prepareInnerNested <- function(o,
                                n,
                                xlim,
                                ylim,
                                smooth,
                                a_in,
                                ...) {
  if (!exists("expsmooth") || !exists("mgks")) {
    expsmooth <- mgks <- function(x) {
      
    }
    stop("Please install the gamFactory package.")
  }
  
  gObj <- o$gObj
  sm <- gObj$smooth[[o$ism]]
  
  si <- sm$xt$si
  mk <- si$margin[[1]]   # the (only) margin of a 1D nested effect
  alpha <- si$alpha
  B <- mk$B
  
  da <- length(alpha)
  prange <- (sm$first.para:sm$last.para)[1:da]
  type <- .nestTypes(si)
  
  # Remove scale parameter of inner transformation
  if (type %in% c("mgks", "exp")) {
    prange <- prange[-1]
    da <- da - 1
    alpha <- alpha[-1]
    if (is.null(B)) {
      B <- diag(1, nrow = length(alpha))
    }
  }
  
  Va <- gObj$Vp[prange, prange, drop = FALSE]
  
  if (type == "si_nexp") {
    if (!smooth) {
      # coef plot (Only alpha_si)
      alpha_center <- mk$alpha_center
      n_si <- mk$n_si
      n_nexp <- mk$n_nexp
      alpha_si <- alpha[n_nexp + seq_len(n_si)]   # parameters are c(alpha_nexp, alpha_si)
      positive_si <- mk$positive_si
      
      if (is.null(alpha_center)) {
        alpha_center <- alpha_si * 0
      }
      
      alpha_si <- drop(mk$B_si %*% (alpha_si + alpha_center))
      
      Va_si <- mk$B_si %*% Va[(n_nexp + 1):(n_nexp + n_si), (n_nexp + 1):(n_nexp + n_si), drop = FALSE] %*% t(mk$B_si)
      se_si <- sqrt(pmax(0, diag(Va_si)))
      
      # Consistently use
      lower <- alpha_si - 2 * se_si
      upper <- alpha_si + 2 * se_si
      
      edf   <- sum(gObj$edf[prange])
      ylabel <- .subEDF(paste0("Inner_coef(", sm$term, ")"), edf)
      #number of edf = length(alpha_si + alpha_nexp)
      #but we only plot alpha_si here
      xlabel <- "Index (SI)"
      
      if (positive_si && !a_in) {
        exp_alpha <- exp(alpha_si)
        sum_exp <- sum(exp_alpha, na.rm = TRUE)
        
        alpha_si <- exp_alpha / sum_exp
        
        lower <- exp(lower) / sum_exp
        upper <- exp(upper) / sum_exp
      }
      
      out <- list(
        fit = unname(alpha_si),
        x = 1:length(alpha_si),
        se = rep(NA, length(alpha_si)), #use NA otherwise there is error in plot function
        lower = unname(lower), 
        upper = unname(upper),
        xlab = xlabel,
        ylab = ylabel,
        main = NULL,
        type = "si_nexpsm_coef"
      )
    } else {
      # smooth eff plot
      xa <- sm$xt$xa #values after smooth
      times <- 1:length(xa)
      edf   <- sum(gObj$edf[prange])
      ylabel <- .subEDF(paste0("expsm(", sm$term, ")"), edf)
      xlabel <- "Index"
      
      out <- list(
        fit = xa,
        x = times,
        xlab = xlabel,
        ylab = ylabel,
        main = NULL,
        type = "si_nexpsm_xa"
      )
    }
  }
  
  if (type == "exp") {
    inner <- expsmooth(
      y = mk$y,
      Xi = mk$W,
      beta = alpha,
      deriv = 1
    )
    fit <- inner$d0
    Jac <- inner$d1
    if (!is.null(mk$times)) {
      fit <- fit[1:max(mk$times)]
      Jac <- Jac[1:max(mk$times), ]
    }
    nobs <- length(fit)
    se <- sqrt(pmax(0, rowSums((Jac %*% Va) * Jac)))
    edf   <- sum(gObj$edf[prange])
    ylabel <- .subEDF(paste0("expsm(", sm$term, ")"), edf)
    xlabel <- "Index"
    
    if (!is.null(xlim)) {
      xlim <- sort(xlim)
      xlim[1] <- max(xlim[1], 1)
      xlim[2] <- min(xlim[2], nobs)
      ii <- which(1:nobs >= xlim[1] & 1:nobs <= xlim[2])
      nobs <- length(ii)
    } else {
      xlim <- c(1, nobs)
      ii <- 1:nobs
    }
    
    out <- list(
      "fit" = fit[ii],
      "x" = ii,
      "se" = se[ii],
      "p.resid" = mk$y[ii],
      "raw" = ii,
      "xlim" = xlim,
      xlab = xlabel,
      ylab = ylabel,
      main = NULL,
      type = "nexpsm"
    )
    
  }
  
  if (type == "si") {
    a0 <- mk$a0
    if (is.null(a0)) {
      a0 <- alpha * 0
    }
    alpha <- drop(B %*% (alpha + a0))
    Va <- B %*% Va %*% t(B)
    se <- sqrt(pmax(0, diag(Va)))
    edf   <- sum(gObj$edf[prange])
    ylabel <- .subEDF(paste0("Inner_coef(", sm$term, ")"), edf)
    xlabel <- "Index"
    
    out <- list(
      "fit" = alpha,
      "x" = 1:da,
      "se" = se,
      xlab = xlabel,
      ylab = ylabel,
      main = NULL,
      type = "si"
    )
    
  }
  
  # NOT CLEAR HOW TO DO THIS WITH DISTANCEs
  # if( type == "mgks" ){
  #   d <- ncol(si$X0)
  #   if( d != 2 ){ return( NULL ) }
  #   # ONLY 2D case handled at the moment!!
  #
  #   if( !is.null(xlim) ) xlim <- sort(xlim) else xlim <- range(si$X[ , 1])
  #   if( !is.null(ylim) ) ylim <- sort(ylim) else ylim <- range(si$X[ , 2])
  #
  #   xx <- rep(seq(xlim[1], xlim[2], length.out = n), n)
  #   yy <- rep(seq(ylim[1], ylim[2], length.out = n), rep(n, n))
  #   X <- cbind(xx, yy)
  #
  #   si$x <- as.matrix(si$x)
  #   if( ncol(si$x) > 1 ){ si$x <- colMeans(si$x) }
  #
  #   inner <- mgks(y = si$x, X = X, X0 = si$X0, beta = alpha[-1], deriv = 1)
  #   fit <- inner$d0
  #   Jac <- inner$d1
  #   se <- sqrt(pmax(0, rowSums((Jac %*% Va[-1, -1, drop = FALSE]) * Jac)))
  #   edf   <- sum(gObj$edf[prange[-1]])
  #
  #   mainlab <- .subEDF(paste0("mgks(", sm$term, ")"), edf)
  #   ylabel <- "X[ , 2]"
  #   xlabel <- "X[ , 1]"
  #   out <- list("fit" = fit, "X" = si$X, "se" = se, "x" = xx, "y" = yy,
  #               "p.resid" = si$x, "X0" = si$X0,
  #               "xlim" = xlim, "ylim" = ylim,
  #               "xlab" = xlabel, "ylab" = ylabel, "main" = mainlab, type = "mgks")
  # }
  
  return(out)
  
}