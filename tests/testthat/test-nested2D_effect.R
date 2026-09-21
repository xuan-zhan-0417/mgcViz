context("plot.nested2D ")

test_that("plot.nested2D", {
  library(mgcViz)
  library(gamFactory)
  skip_if_not("trans_inter_le" %in% getNamespaceExports("gamFactory"),
              "installed gamFactory does not provide the inter_le effect")

  set.seed(2)
  n <- 300
  X_l <- matrix(rnorm(n * 3), n, 3)
  X_e <- cbind(rnorm(n), 1, runif(n))
  colnames(X_e) <- c("y", "x", "x")

  z1 <- as.vector(scale(X_l %*% c(1, -0.5, 0.25)))
  y <- sin(2 * z1) + 0.5 * X_e[, 1] + rnorm(n, sd = 0.1)

  dat <- data.frame(y = y)
  dat$X_l <- I(X_l)
  dat$X_e <- I(X_e)

  fit <- gam_nl(
    list(y ~ s_nest(X_l, X_e, trans = trans_inter_le()), ~ 1),
    family = fam_gaussian(),
    data = dat,
    optimizer = "efs"
  )

  viz <- getViz(fit)
  eff <- sm(viz, 1)

  expect_true(inherits(eff, "nested2D"))

  # Default plot has the same structure as standard 2D smooths
  pl <- plot(eff, n = 20)
  expect_s3_class(pl, "plotSmooth")
  expect_equal(pl$type, "2D")
  expect_equal(nrow(pl$data$fit), 20^2)
  expect_equal(nrow(pl$data$res), n)

  # Effect on the grid agrees with a direct computation from coefficients
  sm1 <- fit$smooth[[1]]
  X <- sm1$xt$basis$evalX(z1 = pl$data$fit$x, z2 = pl$data$fit$y)$X0
  X <- cbind(matrix(0, nrow(X), length(sm1$xt$si$alpha)), X)
  prange <- sm1$first.para:sm1$last.para
  ok <- !is.na(pl$data$fit$z)   # locations too far from the data are set to NA
  expect_true(any(ok) && !all(ok))
  expect_equal(pl$data$fit$z[ok], drop(X %*% coef(fit)[prange])[ok])

  # Separate effects: X1 beta_1 (margin 1), X2 beta_2 (margin 2), X12 beta_3 (interaction)
  # The blocks have p, p and p^2 columns, where p is the dimension of the marginal basis
  prange_sp <- prange[-seq_len(length(sm1$xt$si$alpha))]   # drop the inner parameters
  p <- sqrt(length(prange_sp) + 1) - 1
  expect_equal(p, round(p))
  blocks <- list(margin1 = seq_len(p), margin2 = p + seq_len(p), inter = 2 * p + seq_len(p^2))
  Vp <- fit$Vp[prange_sp, prange_sp]
  
  for (nm in c("margin1", "margin2")) {
    pl1 <- plot(eff, part = nm, n1 = 50)
    expect_equal(pl1$type, "1D")
    d <- pl1$data$fit
    expect_equal(nrow(d), 50)
    cc <- blocks[[nm]]
    X <- if (nm == "margin1") {
      sm1$xt$basis$evalX(z1 = d$x, z2 = rep(0, 50))$X0
    } else {
      sm1$xt$basis$evalX(z1 = rep(0, 50), z2 = d$x)$X0
    }
    X <- X[, cc]
    expect_equal(d$y, drop(X %*% coef(fit)[prange_sp][cc]))
    expect_equal(d$se, sqrt(rowSums((X %*% Vp[cc, cc]) * X)))
    expect_equal(nrow(pl1$data$res), n)
    expect_error(print(pl1 + l_fitLine() + l_ciLine() + l_rug()), NA)
  }
  
  pli <- plot(eff, part = "inter", n = 20)
  expect_equal(pli$type, "2D")
  X <- sm1$xt$basis$evalX(z1 = pli$data$fit$x, z2 = pli$data$fit$y)$X0[, blocks$inter]
  ok <- !is.na(pli$data$fit$z)
  expect_equal(pli$data$fit$z[ok], drop(X %*% coef(fit)[prange_sp][blocks$inter])[ok])
  expect_error(print(pli + l_fitRaster(pTrans = zto1(0.05, 2, 0.1), noiseup = TRUE) + 
                       l_rug() + l_fitContour()), NA)
  
  # The three parts add up to the whole effect
  plf <- plot(eff, n = 20, too.far = 0)
  Xf <- sm1$xt$basis$evalX(z1 = plf$data$fit$x, z2 = plf$data$fit$y)$X0
  expect_equal(plf$data$fit$z, drop(Xf %*% coef(fit)[prange_sp]))
  
  expect_error(plot(eff, part = "wrong"))
  
  # Coefficients of the inner transformations, on the scale of the original covariates
  si <- sm1$xt$si
  V <- fit$Vp
  ia1 <- sm1$first.para - 1 + seq_len(si$na1)
  ia2 <- sm1$first.para - 1 + si$na1 + 1 + seq_len(si$na2)
  
  pc1 <- plot(eff, part = "coef1")
  expect_s3_class(pc1, "plotSmooth")
  expect_equal(pc1$type, c("si", "Factor"))
  expect_equal(pc1$data$fit$y, drop(si$B_1 %*% (coef(fit)[ia1] + si$a0_1)), ignore_attr = TRUE)
  expect_equal(pc1$data$fit$se, sqrt(diag(si$B_1 %*% V[ia1, ia1] %*% t(si$B_1))), ignore_attr = TRUE)
  expect_equal(nrow(pc1$data$fit), si$na1)
  
  pc2 <- plot(eff, part = "coef2")
  expect_equal(pc2$data$fit$y, drop(si$B_2 %*% coef(fit)[ia2]), ignore_attr = TRUE)
  expect_equal(pc2$data$fit$se, sqrt(diag(si$B_2 %*% V[ia2, ia2] %*% t(si$B_2))), ignore_attr = TRUE)
  expect_equal(nrow(pc2$data$fit), si$na2)
  
  # With original-scale weights we recover the single index that the model stores
  z1 <- drop(scale(X_l, center = TRUE, scale = FALSE) %*% pc1$data$fit$y)
  expect_equal(z1 - mean(z1), drop(sm1$xt$xa[, 1] - mean(sm1$xt$xa[, 1])))
  
  expect_error(print(pc1 + l_ciBar() + l_fitPoints()), NA)
  expect_error(print(pc2), NA)
  
  # Layers and arguments
  expect_error(print(pl + l_fitRaster() + l_fitContour() + l_rug() + l_points()), NA)
  expect_error(print(plot(eff, n = 20) + l_dens(type = "joint") + l_fitContour()), NA)
  expect_error(print(plot(eff, n = 20, too.far = 0, trans = exp, unconditional = TRUE) + l_fitRaster()), NA)

  # Global plot
  expect_error(print(plot(viz), pages = 1), NA)
})
