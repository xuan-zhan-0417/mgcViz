##########
# WHAT THIS FILE DOES
# --------------------
# Prepares the data behind every plot.nested2D() plot, for gamFactory's "inter_le"
# effect, s(si(x), exp(x)). That is the only effect type handled here.
#
# THE EFFECT, IN ONE PICTURE
# --------------------------
#   z1 = si(x)    -> margin 1 (a single index)
#   z2 = exp(x)   -> margin 2 (an adaptive exponential smooth)
#
#   model matrix:  X0 = [ X1(z1) , X2(z2) , X1(z1) (x)_row X2(z2) ]
#                          margin 1   margin 2      interaction block
#
#   effect      =  X1(z1) beta_1  +  X2(z2) beta_2  +  X12(z1,z2) beta_3
#
#   (columns for the inner si(x)/exp(x) parameters come before X0 in the real model
#   matrix, but the effect itself does not depend on them -- so we drop them below.)
#
# ONE FUNCTION PER "part" ARGUMENT OF plot.nested2D()
# ----------------------------------------------------
#   .prepareNested2D        "full" / "margin1" / "margin2" / "inter"
#   .prepareNested2DCoef     "coef1" / "coef2"   -- the si(x) / exp(x) weights
#   .prepareNested2DSmooth   "smooth"            -- data before/after exp smoothing
#
# p1, p2 below are the number of columns of X1(z1) and X2(z2) respectively (so that
# X0 has p1 + p2 + p1*p2 columns). gamFactory stores them on si$p1 / si$p2 (added
# specifically so mgcViz does not have to reverse-engineer them).
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

  # coef(gObj) for this smooth = [inner si(x)/exp(x) params, then beta_1, beta_2, beta_3].
  # We only want the beta_* part, so drop the inner-parameter columns:
  prange <- (sm$first.para:sm$last.para)[-seq_len(length(si$alpha))]

  # Which columns of X0 (and of beta_1/beta_2/beta_3) does this `part` need?
  if (is.null(si$p1) || is.null(si$p2)) {
    stop("si$p1 / si$p2 not found: please update gamFactory to a version whose ",
         ".build_n_inter_bspline_basis() stores the marginal basis dimensions on si.")
  }
  p1 <- si$p1; p2 <- si$p2
  cols <- switch(part,
                 "full"    = seq_len(p1 + p2 + p1 * p2),
                 "margin1" = seq_len(p1),
                 "margin2" = p1 + seq_len(p2),
                 "inter"   = p1 + p2 + seq_len(p1 * p2))
  prange <- prange[cols]

  # Grab the coefficients, covariance and edf for just those columns
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

  # ---- 1D case: margin1 or margin2 -> X_k(z_k) beta_k, drawn as a curve in z_k ----
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

  # ---- 2D case: full effect or interaction -> a surface over a (z1, z2) grid ----
  # x (= z1) runs fastest across the grid, matching .preparePlotSmooth2D's convention
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
# COEFFICIENTS OF THE INNER TRANSFORMATION (part = "coef1" or "coef2")
# ----------------------------------------------------------------------
#   margin = 1: the weights of the single index si(x)
#   margin = 2: the smoothing-rate coefficients of exp(x)
#               (NOT alpha_scale -- that is a separate scalar, not one of these)
#
# Why the B / a0 business: the model actually fits these on a REPARAMETRISED, rotated
# version of the covariates (B is the rotation matrix used to diagonalise the penalty).
# To report the coefficients on the ORIGINAL covariate scale, we rotate back:
#
#   original-scale coefficients  =  B %*% (alpha + a0)
#   their covariance matrix      =  B %*% V %*% t(B)
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
# MARGIN 2's DATA, BEFORE AND AFTER SMOOTHING (part = "smooth")
# -----------------------------------------------------------------
# The raw series and the smoothing-rate design matrix are already saved on the smooth
# object (si$y_raw, si$W_2). We just re-run the exponential smooth with the fitted
# coefficients to get the "after" series.
#
# Scale note: this is the RAW smoothed series, on the same scale as the data -- it is
# NOT z2 (z2 = exp(alpha_scale) * (this series - its mean)). We deliberately skip that
# last step, so the "before" and "after" lines can be compared on one common scale.
#
# x-axis note: it is just position-in-the-series (1, 2, 3, ...). If one response row
# corresponds to several observations of the series (a high-frequency case), this
# series is longer than the number of response rows.
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
