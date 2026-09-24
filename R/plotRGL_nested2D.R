#'
#' Visualizing 2D nested effects in 3D
#'
#' @description This method plots an interactive 3D representation of the whole effect
#'              (all three parts added together, see [plot.nested2D]) of a two dimensional
#'              nested smooth effect fitted with \code{gamFactory}. At the moment only
#'              effects of type \code{"inter_le"}, that is \code{s(si(x), exp(x))}, are
#'              supported. The output is built with the same internal renderer used by
#'              [plotRGL.mgcv.smooth.2D], so it has the same style: a red surface for the
#'              fitted effect and, if \code{se = TRUE}, two blue wireframe surfaces at
#'              \code{+/- se.mult} standard errors from it.
#' @param x a smooth effect object, extracted using [mgcViz::sm].
#' @param se when TRUE (default) upper and lower surfaces are added to the plot at \code{se.mult}
#'           (see below) standard deviations for the fitted surface.
#' @param n sqrt of the number of grid points used to compute the effect plot.
#' @param residuals if TRUE, then the (whole model) residuals will be added, as coloured spheres
#'                  on a plane below the fitted surface, at the observed \code{(si(x), exp(x))}
#'                  locations. See [plotRGL.mgcv.smooth.2D].
#' @param type the type of residuals that should be plotted. See [residuals.gamViz].
#' @param maxpo maximum number of residuals points that will be plotted.
#'              If number of datapoints > \code{maxpo}, then a subsample of \code{maxpo} points will be taken.
#' @param too.far if greater than 0 then this is used to determine when a location is too far
#'               from data to be plotted. The data are scaled into the unit square before deciding
#'               what to exclude, and too.far is a distance within the unit square.
#'               Unlike [plot.nested2D], the default here is 0 (no exclusion), as in
#'               [plotRGL.mgcv.smooth.2D].
#' @param xlab if supplied then this will be used as the x label of the plot.
#' @param ylab if supplied then this will be used as the y label of the plot.
#' @param main used as title for the plot if supplied.
#' @param xlim if supplied then this pair of numbers are used as the limits of \code{si(x)}
#'             (x axis) for the plot.
#' @param ylim if supplied then this pair of numbers are used as the limits of \code{exp(x)}
#'             (y axis) for the plot.
#' @param se.mult a positive number which will be the multiplier of the standard errors
#'                when calculating standard error surfaces.
#' @param trans currently unused (kept only for signature parity with
#'              [plotRGL.mgcv.smooth.2D], which accepts it but does not apply it either).
#' @param unconditional if \code{TRUE} then the smoothing parameter uncertainty corrected covariance
#'                      matrix is used to compute uncertainty bands, if available.
#'                      Otherwise the bands treat the smoothing parameters as fixed.
#' @param ... currently unused.
#' @return Returns \code{NULL} invisibly.
#' @seealso [plotRGL.mgcv.smooth.2D], [plot.nested2D].
#' @name plotRGL.nested2D
#' @rdname plotRGL.nested2D
#' @examples
#' \dontrun{
#' library(mgcViz)
#' library(gamFactory)
#' library(rgl)
#'
#' # fit <- gam_nl(list(y ~ s_nest(X_l, X_e, trans = trans_inter_le()), ~ 1),
#' #               data = dat, family = fam_gaussian(), optimizer = "efs")
#' b <- getViz(fit)
#'
#' plotRGL(sm(b, 1))
#'
#' # Also show the (whole model) residuals
#' rgl.close()
#' plotRGL(sm(b, 1), residuals = TRUE)
#' }
#' @export plotRGL.nested2D
#' @export
#'
plotRGL.nested2D <- function(x, se = TRUE, n = 40, residuals = FALSE, type = "auto",
                             maxpo = 1e3, too.far = 0, xlab = NULL, ylab = NULL,
                             main = NULL, xlim = NULL, ylim = NULL, se.mult = 1,
                             trans = identity, unconditional = FALSE, ...) {

  pack <- requireNamespace("rgl", quietly = TRUE)
  if (!pack) {
    message("Please install the rgl package to use this function.")
    return(NULL)
  }

  if (type == "auto") { type <- .getResTypeAndMethod(x$gObj$family)$type }

  # 1) Prepare data: the whole effect, on a grid in (si(x), exp(x)) space
  P <- .prepareNested2D(o = x, part = "full", n = n, n1 = 100, xlim = xlim, ylim = ylim,
                        too.far = too.far, unconditional = unconditional)

  if (!is.null(xlab)) { P$xlab <- xlab }
  if (!is.null(ylab)) { P$ylab <- ylab }
  if (!is.null(main)) { P$main <- main }

  # Locations too far from the data are not plotted, as in the standard 2D case (.createP)
  P$fit[P$exclude] <- NA
  P$se[P$exclude] <- NA
  P$se <- P$se * se.mult

  R <- list()
  if (residuals) {
    R <- .getResidualsPlotRGL(gamObj = x$gObj, X = P$raw, type = type, maxpo = maxpo,
                              xlimit = xlim, ylimit = ylim, exclude = rep(FALSE, nrow(P$raw)))
    P$raw <- R$raw
  }

  # 2) Actual plotting: reuse plotRGL.mgcv.smooth.2D's own internal renderer as-is, so this
  # has exactly the same look (and the same quirks -- e.g. `trans` is accepted here only for
  # signature parity with plotRGL.mgcv.smooth.2D; like that function, it is not actually
  # applied to the surfaces).
  P$plotCI <- se
  .plotRGL.mgcv.smooth.2D(P = P, trans = trans, res = R$res)

}
