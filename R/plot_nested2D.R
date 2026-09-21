#'
#' Plotting two dimensional nested effects
#'
#' @description Plotting method for two dimensional nested smooth effects fitted with
#'              \code{gamFactory}. At the moment only effects of type \code{"inter_le"},
#'              that is \code{s(si(x), exp(x))}, an interaction between a single index and
#'              an adaptive exponential smooth, are supported. The effect is plotted as a
#'              function of the two inner indices, so the x axis is the single index
#'              \code{si(x)} and the y axis is the exponential smooth \code{exp(x)}.
#'              The output is of the same form as that of [plot.mgcv.smooth.2D] (or 
#'              [plot.mgcv.smooth.1D] for the margins), hence layers such as [l_fitRaster], 
#'              [l_fitContour], [l_rug], [l_points], [l_dens], [l_fitLine] and [l_ciLine] 
#'              can be added in the usual way.
#' @details The effect is \code{X1(z1) beta_1 + X2(z2) beta_2 + X12(z1, z2) beta_3}, where 
#'          \code{z1 = si(x)}, \code{z2 = exp(x)}, \code{X12} is the row-wise product of 
#'          \code{X1} and \code{X2}, and \code{beta_1, beta_2, beta_3} are the corresponding
#'          blocks of spline coefficients. The \code{part} argument selects what is plotted.
#'
#' @param x a smooth effect object, extracted using [mgcViz::sm].
#' @param part which part of the effect should be plotted:
#'             \itemize{
#'               \item{\code{"full"}}{ the whole effect, as a function of \code{(z1, z2)} (2D plot).}
#'               \item{\code{"margin1"}}{ the effect of the first margin, \code{X1 beta_1}, 
#'                     as a function of \code{z1 = si(x)} (1D plot).}
#'               \item{\code{"margin2"}}{ the effect of the second margin, \code{X2 beta_2}, 
#'                     as a function of \code{z2 = exp(x)} (1D plot).}
#'               \item{\code{"inter"}}{ the interaction, \code{X12 beta_3}, as a function of 
#'                     \code{(z1, z2)} (2D plot).}
#'               \item{\code{"coef1"}}{ the coefficients of the inner transformation of the first margin: 
#'                     the weights of the covariates in the single index \code{si(x)}, with 
#'                     confidence intervals (point and error bar plot).}
#'               \item{\code{"coef2"}}{ the coefficients of the inner transformation of the second margin: 
#'                     the coefficients of the covariates modelling the smoothing rate of \code{exp(x)},
#'                     with confidence intervals. The scaling parameter of \code{exp(x)} is not included.}
#'               \item{\code{"smooth"}}{ the data of the second margin before (grey line) and after 
#'                     (black line) the exponential smoothing, on the same plot. The x axis is the position 
#'                     in the series to be smoothed and both series are on the scale of the raw data.
#'                     Use \code{xlim} to zoom in on part of the series. No layers are available for 
#'                     this plot, but other \code{ggplot2} components can be added.}
#'             }
#'             The coefficients are on the scale of the original covariates. Their confidence intervals
#'             are based on the Bayesian covariance matrix of the fit, conditional on the 
#'             other parameters being estimated (the inner and outer parameters are estimated jointly, 
#'             but the smoothing parameters are treated as fixed unless \code{unconditional = TRUE}).
#' @param n sqrt of the number of grid points used to compute the 2D effect plots.
#' @param n1 number of grid points used to compute the 1D effect plots (margins).
#' @param xlim if supplied then this pair of numbers are used as the limits of the x axis 
#'             (\code{si(x)} for the 2D plots and for \code{"margin1"}, \code{exp(x)} for \code{"margin2"}).
#' @param ylim if supplied then this pair of numbers are used as the limits of the exponential
#'             smooth (y axis) for the 2D plots. Not used for the 1D plots.
#' @param maxpo maximum number of points that will be used by layers such as
#'              \code{l_rug()} and \code{l_points()}. If number of datapoints > \code{maxpo},
#'              then a subsample of \code{maxpo} points will be taken.
#' @param too.far if greater than 0 then this is used to determine when a location is too far
#'               from data to be plotted (2D plots only). The data are scaled into the unit square before
#'               deciding what to exclude, and too.far is a distance within the unit square.
#'               Setting to zero can make plotting faster for large datasets, but care
#'               is then needed with interpretation of plots.
#' @param trans monotonic function to apply to the effect before plotting.
#'              Monotonicity is not checked.
#' @param unconditional if \code{TRUE} then the smoothing parameter uncertainty corrected covariance
#'                      matrix is used to compute uncertainty bands, if available.
#'                      Otherwise the bands treat the smoothing parameters as fixed.
#' @param ... currently unused.
#' @return An object of class \code{c("plotSmooth", "gg")}.
#' @seealso [plot.mgcv.smooth.2D], [plot.mgcv.smooth.1D], [plot.nested1D].
#' @name plot.nested2D
#' @rdname plot.nested2D
#' @examples
#' \dontrun{
#' library(mgcViz)
#' library(gamFactory)
#'
#' # fit <- gam_nl(list(y ~ s_nest(X_l, X_e, trans = trans_inter_le()), ~ 1),
#' #               data = dat, family = fam_gaussian(), optimizer = "efs")
#' b <- getViz(fit)
#'
#' # Effect with contour and rug for the observed (si(x), exp(x)) pairs
#' plot(sm(b, 1)) + l_fitRaster() + l_fitContour() + l_rug()
#'
#' # Effect of the two margins, and of their interaction
#' plot(sm(b, 1), part = "margin1") + l_fitLine() + l_ciLine() + l_rug()
#' plot(sm(b, 1), part = "margin2") + l_fitLine() + l_ciLine() + l_rug()
#' plot(sm(b, 1), part = "inter") + l_fitRaster(pTrans = zto1(0.05, 2, 0.1), noiseup = TRUE) +
#'   l_rug() + l_fitContour()
#'
#' # Coefficients of the inner transformations, with confidence intervals
#' plot(sm(b, 1), part = "coef1") + l_ciBar() + l_fitPoints()
#' plot(sm(b, 1), part = "coef2") + l_ciBar() + l_fitPoints()
#'
#' # Data of the second margin before and after exponential smoothing
#' plot(sm(b, 1), part = "smooth")
#' plot(sm(b, 1), part = "smooth", xlim = c(1, 200))
#'
#' # Opacity proportional to the significance of the effect
#' plot(sm(b, 1)) + l_fitRaster(pTrans = zto1(0.05, 2, 0.1)) + l_fitContour() + l_points()
#'
#' # Joint density of the two indices, filling the plot
#' plot(sm(b, 1)) + l_dens(type = "joint") + l_points() + l_fitContour() +
#'   coord_cartesian(expand = FALSE)
#' }
#' @export plot.nested2D
#' @export
#'
plot.nested2D <- function(x, part = c("full", "margin1", "margin2", "inter", "coef1", "coef2", "smooth"), 
                          n = 40, n1 = 100, xlim = NULL, ylim = NULL, maxpo = 1e4,
                          too.far = 0.1, trans = identity, unconditional = FALSE, ...) {
  
  part <- match.arg(part)
  
  # Coefficients of the inner transformations: same plot as the weights of single index effects
  if (part %in% c("coef1", "coef2")) {
    P <- .prepareNested2DCoef(o = x, margin = as.numeric(substr(part, 5, 5)), 
                              unconditional = unconditional)
    out <- .plot.inner.nested.smooth.1D(P = P, trans = trans, maxpo = maxpo, ci = TRUE)
    class(out) <- c("plotSmooth", "gg")
    return(out)
  }
  
  # Data before and after exponential smoothing: both series on the same plot
  if (part == "smooth") {
    P <- .prepareNested2DSmooth(o = x, xlim = xlim, ...)
    out <- .plot.nested2D.smooth(P = P, trans = trans)
    class(out) <- c("plotSmooth", "gg")
    return(out)
  }
  
  # 1) Prepare data
  P <- .prepareNested2D(o = x, part = part, n = n, n1 = n1, xlim = xlim, ylim = ylim, 
                        too.far = too.far, unconditional = unconditional, ...)
  
  # 2) Produce output object: same structure as for standard 1D and 2D smooths
  if (part %in% c("margin1", "margin2")) {
    out <- .plot.outer.nested.1D(x = NULL, P = P, trans = trans, maxpo = maxpo)
  } else {
    out <- .plot.mgcv.smooth.2D(x = NULL, P = P, trans = trans, maxpo = maxpo)
  }
  
  class(out) <- c("plotSmooth", "gg")
  
  return(out)
  
}

########################
#' @noRd
.plot.nested2D.smooth <- function(P, trans) {
  
  lev <- c("Before smoothing", "After smoothing")
  
  .dat <- list()
  .dat$fit <- data.frame(x = P$x, y = P$fit, ty = trans(P$fit), se = NA)
  .dat$res <- data.frame(x = P$x, y = P$raw, sub = TRUE)
  .dat$misc <- list(trans = trans)
  
  # Raw series first, so that the smoothed one is drawn on top of it
  .long <- data.frame(x = rep(P$x, 2), 
                      y = c(trans(P$raw), trans(P$fit)),
                      series = factor(rep(lev, each = length(P$x)), levels = lev))
  
  .pl <- ggplot(data = .long, mapping = aes(x = x, y = y, colour = series)) +
    geom_line(na.rm = TRUE) +
    scale_colour_manual(values = c("grey60", "black"), name = NULL) +
    labs(title = P$main, x = P$xlab, y = P$ylab) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  return(list("ggObj" = .pl, "data" = .dat, "type" = c("nexp", "Series")))
  
}
