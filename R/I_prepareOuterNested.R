
##########
# Internal method that prepare non-linear effects plots (only single index at the moment)
#
.prepareOuterNested <- function(o, n, xlim, ...){
  
  if(!exists("expsmooth") || !exists("mgks") ){
    expsmooth <- mgks <- function(x){}
    stop("Please install the gamFactory package.")
  }
  
  gObj <- o$gObj
  sm <- gObj$smooth[[ o$ism ]]
  si <- sm$xt$si
  mk <- si$margin[[1]]   # the (only) margin of a 1D nested effect
  dsi <- length( si$alpha )
  type <- .nestTypes(si)
  
  # x axis: the inner index at the data on the scale of the inner transformation, i.e. before the
  # exp(scale) * (. - xm) applied by exp and mgks margins; rescale() maps it back to the index
  sc <- if( is.null(mk$iscale) ) 1 else exp(si$alpha[mk$iscale])
  xm <- if( is.null(mk$iscale) ) 0 else mk$xm
  rescale <- function(x){ sc * (x - xm) }
  raw <- sm$xt$xa[ , 1] / sc + xm
  trnam <- switch(type, "si" = "proj", "exp" = "expsm", "mgks" = "mgks", "si_nexp" = "si_nexpsm")

  prange <- (sm$first.para:sm$last.para)[-(1:dsi)]
  beta <- coef( gObj )[ prange ]
  
  # Generate x sequence for prediction
  if (is.null(xlim)){ 
    xlim <- range(raw)
  } 
  xx <- seq(xlim[1], xlim[2], length = n) 
  
  # Compute outer model matrix
  X <- sm$xt$basis$evalX(z1 = rescale(xx), deriv = 0)$X0
  
  fit <- X %*% beta
  
  se <- sqrt(pmax(0, rowSums((X %*% gObj$Vp[prange, prange, drop = FALSE]) * X)))
  
  edf   <- sum(gObj$edf[prange])
  ylabel <- .subEDF(paste0("s(",trnam,"(", sm$term, "))"), edf)
  xlabel <- paste0(trnam,"(", sm$term, ")")
  out <- list("fit" = fit, "x" = xx, "se" = se, "raw" = raw, "xlim" = xlim,
              xlab = xlabel, ylab = ylabel, main = NULL)
  return(out)
  
}