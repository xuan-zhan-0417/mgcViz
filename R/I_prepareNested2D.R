##########
# Internal methods that prepare plots of two dimensional nested effects
# (only "inter_le", i.e. s(si(x), exp(x)), at the moment).
#
# The effect is a smooth function of the two inner indices (z1, z2) = (si(x), exp(x)),
# with model matrix (before the columns of the inner parameters, which do not enter
# the effect as a function of z1, z2)
#
#    X0 = [ X1(z1) , X2(z2) , X1(z1) (x)_row X2(z2) ]
#
# and coefficients [beta_1, beta_2, beta_3]. Here we can plot the whole effect
# (part = "full"), the two margins X1 beta_1 and X2 beta_2 (part = "margin1", "margin2")
# or the interaction X12 beta_3 (part = "inter"). The coefficients of the inner 
# transformations are handled by .prepareNested2DCoef (part = "coef1", "coef2") and
# the series before/after the exponential smoothing by .prepareNested2DSmooth (part = "smooth").
#
.prepareNested2D <- function(o, part, n, n1, xlim, ylim, too.far, unconditional, ...) {

  gObj <- o$gObj
  sm <- gObj$smooth[[o$ism]]
  type <- class(o)[1]

  if (type != "inter_le") {
    stop("plot.nested2D is only available for effects of type \"inter_le\", not \"", type, "\".")
  }

  si <- sm$xt$si
  xa <- sm$xt$xa   # inner indices (z1, z2) at the observed covariates
  if (is.null(xa) || ncol(xa) != 2) {
    stop("The smooth effect does not contain the inner indices: was the model fitted with gam_nl()?")
  }
  raw <- data.frame(x = as.numeric(xa[, 1]), y = as.numeric(xa[, 2]))
  evalX <- function(z1, z2) sm$xt$basis$evalX(z1 = z1, z2 = z2, deriv = 0)$X0

  # Spline coefficients [beta_1, beta_2, beta_3]: the leading ones belong to the inner parameters
  prange <- (sm$first.para:sm$last.para)[-seq_len(length(si$alpha))]

  # Columns of X0 (and of [beta_1, beta_2, beta_3]) that we need
  d <- .nested2DBlockDims(evalX, raw)
  cols <- switch(part,
                 "full"    = seq_len(d$p1 + d$p2 + d$p1 * d$p2),
                 "margin1" = seq_len(d$p1),
                 "margin2" = d$p1 + seq_len(d$p2),
                 "inter"   = d$p1 + d$p2 + seq_len(d$p1 * d$p2))
  prange <- prange[cols]

  # Coefficients, covariance matrix and edf of the required part of the effect
  V <- gObj$Vp
  if (unconditional) {
    if (is.null(gObj$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      V <- gObj$Vc
    }
  }
  beta <- coef(gObj)[prange]
  Vb <- V[prange, prange, drop = FALSE]
  edf <- sum(gObj$edf[prange])

  t1 <- sm$term[1]
  t2 <- sm$term[2]

  # 1D effect of one of the margins: X_k(z_k) beta_k as a function of z_k
  if (part %in% c("margin1", "margin2")) {
    k <- if (part == "margin1") 1 else 2
    z <- raw[[k]]
    lim <- if (is.null(xlim)) range(z) else sort(xlim)
    zz <- seq(lim[1], lim[2], length.out = n1)
    zo <- rep(median(raw[[3 - k]]), n1)  # other index: it does not affect columns `cols`
    X <- (if (k == 1) evalX(zz, zo) else evalX(zo, zz))[, cols, drop = FALSE]

    se <- sqrt(pmax(0, rowSums((X %*% Vb) * X)))

    out <- list("fit" = drop(X %*% beta), "x" = zz, "se" = se, "raw" = z, "xlim" = lim,
                "xlab" = if (k == 1) paste0("si(", t1, ")") else paste0("exp(", t2, ")"),
                "ylab" = paste0("f", k, "(", if (k == 1) "si(" else "exp(",
                                if (k == 1) t1 else t2, "), ", round(edf, 2), ")"),
                "main" = NULL)
    return(out)
  }

  # 2D effect: whole effect or interaction, on a grid in (z1, z2) space. Here x (= z1)
  # runs fastest, as in .preparePlotSmooth2D
  n <- max(10, n)
  xlim <- if (is.null(xlim)) range(raw$x) else sort(xlim)
  ylim <- if (is.null(ylim)) range(raw$y) else sort(ylim)
  xm <- seq(xlim[1], xlim[2], length.out = n)
  ym <- seq(ylim[1], ylim[2], length.out = n)
  grid <- expand.grid(x = xm, y = ym)

  X <- evalX(grid$x, grid$y)[, cols, drop = FALSE]
  fit <- drop(X %*% beta)
  se <- sqrt(pmax(0, rowSums((X %*% Vb) * X)))

  # Do not plot locations that are too far from the data
  exclude <- if (too.far > 0) {
    exclude.too.far(grid$x, grid$y, raw$x, raw$y, dist = too.far)
  } else {
    rep(FALSE, n^2)
  }

  main <- if (part == "full") {
    .subEDF(paste0("s(si(", t1, "),exp(", t2, "))"), edf)
  } else {
    .subEDF(paste0("f12(si(", t1, "),exp(", t2, "))"), edf)
  }

  out <- list("fit" = fit, "x" = xm, "y" = ym, "se" = se, "raw" = raw,
              "xlim" = xlim, "ylim" = ylim, "exclude" = exclude,
              "xlab" = paste0("si(", t1, ")"),
              "ylab" = paste0("exp(", t2, ")"),
              "main" = main)
  return(out)

}

##########
# Number of columns of the marginal bases, p1 and p2, in X0 = [X1, X2, X1 (x)_row X2].
# They are not stored in the smooth object, so we work them out from the fact that
# ncol(X0) = (p1 + 1) * (p2 + 1) - 1, X1 depends on z1 only and X2 on z2 only. Hence
# the right split is the (unique) one such that the first p1 columns do not change with z2
# and the following p2 columns do not change with z1.
#
.nested2DBlockDims <- function(evalX, raw) {

  m <- 25
  z1 <- seq(min(raw$x), max(raw$x), length.out = m)
  z2 <- seq(min(raw$y), max(raw$y), length.out = m)

  XA <- evalX(z1, rep(median(raw$y), m))  # z1 varies, z2 fixed
  XB <- evalX(rep(median(raw$x), m), z2)  # z2 varies, z1 fixed

  ntot <- ncol(XA)
  isConst <- function(X) all(abs(sweep(X, 2, X[1, ])) < 1e-10)

  cand <- NULL
  for (p1 in seq_len(ntot)) {
    if ((ntot + 1) %% (p1 + 1) != 0) { next }
    p2 <- (ntot + 1) / (p1 + 1) - 1
    if (p2 < 1) { next }
    if (isConst(XB[, seq_len(p1), drop = FALSE]) &&
        isConst(XA[, p1 + seq_len(p2), drop = FALSE])) {
      cand <- rbind(cand, c(p1, p2))
    }
  }

  if (is.null(cand) || nrow(cand) != 1) {
    stop("Could not determine the dimensions of the marginal bases of the effect.")
  }

  return(list("p1" = cand[1, 1], "p2" = cand[1, 2]))

}

##########
# Coefficients of the inner transformation of one of the two margins: the weights of the single
# index si(x) (margin = 1) or the coefficients of the smoothing rate of exp(x) (margin = 2).
# The scaling parameter of exp(x), alpha_scale, is not included.
#
# The coefficients stored in the fit refer to the reparametrised model matrix X B, where B 
# diagonalises the penalty (and, for margin 1, a0 is a fixed shift), so we go back to the
# coefficients of the original covariates: B (alpha + a0), with covariance matrix B V B^T. 
#
.prepareNested2DCoef <- function(o, margin, unconditional) {
  
  gObj <- o$gObj
  sm <- gObj$smooth[[o$ism]]
  type <- class(o)[1]
  
  if (type != "inter_le") {
    stop("plot.nested2D is only available for effects of type \"inter_le\", not \"", type, "\".")
  }
  
  si <- sm$xt$si
  na1 <- si$na1
  na2 <- si$na2
  prange <- sm$first.para:sm$last.para
  
  if (margin == 1) {
    idx <- prange[seq_len(na1)]
    B <- si$B_1
    a0 <- si$a0_1
    term <- sm$term[1]
    main <- paste0("Weights of the single index si(", term, ")")
  } else {
    idx <- prange[na1 + 1 + seq_len(na2)]   # skip alpha_1 and alpha_scale
    B <- si$B_2
    a0 <- NULL
    term <- sm$term[2]
    main <- paste0("Smoothing rate coefficients of exp(", term, ")")
  }
  if (is.null(a0)) { a0 <- numeric(length(idx)) }
  
  V <- gObj$Vp
  if (unconditional) {
    if (is.null(gObj$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      V <- gObj$Vc
    }
  }
  
  alpha <- drop(B %*% (coef(gObj)[idx] + a0))
  Va <- B %*% V[idx, idx, drop = FALSE] %*% t(B)
  se <- sqrt(pmax(0, diag(Va)))
  edf <- sum(gObj$edf[idx])
  
  out <- list("fit" = unname(alpha), "x" = seq_along(alpha), "se" = unname(se),
              "xlab" = "Index", 
              "ylab" = .subEDF(paste0("Inner_coef(", term, ")"), edf), 
              "main" = main, "type" = "si")
  return(out)
  
}

##########
# The data of margin 2 before and after the exponential smoothing. The (possibly high frequency)
# series to be smoothed and the design matrix of the smoothing rate are stored in the smooth 
# object, so we just smooth the former using the fitted coefficients. The smoothed series is 
# returned on the same scale as the raw data (i.e. before the scaling by exp(alpha_scale) and the 
# centring that are applied to obtain the index z2). The time axis is the position in the series;
# if there are several observations per row of the response, the series is longer than the data.
#
.prepareNested2DSmooth <- function(o, xlim, ...) {
  
  gObj <- o$gObj
  sm <- gObj$smooth[[o$ism]]
  type <- class(o)[1]
  
  if (type != "inter_le") {
    stop("plot.nested2D is only available for effects of type \"inter_le\", not \"", type, "\".")
  }
  
  if (!exists("expsmooth")) {
    stop("Please install the gamFactory package.")
  }
  
  si <- sm$xt$si
  prange <- sm$first.para:sm$last.para
  alpha_2 <- coef(gObj)[prange[si$na1 + 1 + seq_len(si$na2)]]   # skip alpha_1 and alpha_scale
  
  raw <- si$y_raw
  smoothed <- drop(expsmooth(y = raw, Xi = si$W_2, beta = alpha_2)$d0)
  
  ii <- seq_along(raw)
  if (!is.null(xlim)) {
    xlim <- sort(xlim)
    ii <- ii[ii >= xlim[1] & ii <= xlim[2]]
  }
  
  out <- list("x" = ii, "raw" = raw[ii], "fit" = smoothed[ii], 
              "xlim" = range(ii),
              "xlab" = "Time index", 
              "ylab" = "Value", 
              "main" = paste0("Margin 2 (", sm$term[2], "): data before and after smoothing"))
  return(out)
  
}
