## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(20221017)

## ----cran, eval=FALSE---------------------------------------------------------
#  install.packages("ktweedie")

## ----install_git, eval=FALSE--------------------------------------------------
#  devtools::install_github("ly129/ktweedie")

## ----setup--------------------------------------------------------------------
library(ktweedie)

## ----data, cache = FALSE------------------------------------------------------
data(dat)
x <- dat$x
y <- dat$y

## ----ktd_estimate1, cache = FALSE---------------------------------------------
fit.ktd <- ktd_estimate(x = x,
                        y = y,
                        kern = rbfdot(sigma = 0.1),
                        lam1 = c(0.01, 0.1, 1))
str(fit.ktd$estimates)

## ----sktd_est, cache = FALSE--------------------------------------------------
fit.sktd <- ktd_estimate(x = x,
                         y = y,
                         kern = rbfdot(sigma = 0.1),
                         lam1 = 5,
                         sparsity = TRUE,
                         lam2 = 1)

## ----sktd_wts, cache = FALSE--------------------------------------------------
fit.sktd$estimates[[1]]$weight

## ----laplace-kernel-----------------------------------------------------------
laplacedot(sigma = 1)

## ----one-d-cv, cache = FALSE--------------------------------------------------
ktd.cv1d <- ktd_cv(x = x,
                   y = y,
                   kern = laplacedot(sigma = 0.1),
                   lambda = c(0.0001, 0.001, 0.01, 0.1, 1),
                   nfolds = 5,
                   loss = "LL")
ktd.cv1d

## ----two-d-cv, cache = FALSE--------------------------------------------------
ktd.cv2d <- ktd_cv2d(x = x,
                     y = y,
                     kernfunc = laplacedot,
                     lambda = c(1e-5, 1e0),
                     sigma = c(1e-5, 1e0),
                     nfolds = 5,
                     ncoefs = 10,
                     loss = "MAD")
ktd.cv2d

## ----ktd_fit, cache = FALSE---------------------------------------------------
ktd.fit <- ktd_estimate(x = x,
                        y = y,
                        kern = laplacedot(sigma = ktd.cv2d$Best_sigma),
                        lam1 = ktd.cv2d$Best_lambda)
str(ktd.fit$estimates)

## ----sktd_fit, cache = FALSE--------------------------------------------------
sktd.cv2d <- ktd_cv2d(x = x,
                      y = y,
                      kernfunc = rbfdot,
                      lambda = c(1e-3, 1e0),
                      sigma = c(1e-3, 1e0),
                      nfolds = 5,
                      ncoefs = 10,
                      loss = "LL")

sktd.fit <- ktd_estimate(x = x,
                         y = y,
                         kern = rbfdot(sigma = sktd.cv2d$Best_sigma),
                         lam1 = sktd.cv2d$Best_lambda,
                         sparsity = TRUE,
                         lam2 = 1,
                         ftol = 1e-3,
                         partol = 1e-3,
                         innerpartol = 1e-5)

## ----fitting, cache = FALSE---------------------------------------------------
ktd.pred <- ktd_predict(ktd.fit, type = "response")
head(ktd.pred$prediction)

## ----fitting_new, cache = FALSE-----------------------------------------------
# Use a subset of the original x as newdata.
newdata <- x[1:6, ]
ktd.pred.new <- ktd_predict(ktd.fit,
                            newdata = newdata,
                            type = "response")
sktd.pred.new <- ktd_predict(sktd.fit,
                             newdata = newdata,
                             type = "response")
data.frame(ktweedie = ktd.pred.new$prediction,
           sktweedie = sktd.pred.new$prediction)

## ----solution-path, cache = FALSE, fig.height = 8, fig.width = 8--------------
nlam2 <- 10
lam2.seq <- 20 * 0.8^(1:nlam2 - 1)
wts <- matrix(NA, nrow = nlam2, ncol = ncol(x))
for (i in 1:nlam2) {
  sktd.tmp <- ktd_estimate(x = x,
                           y = y,
                           kern = rbfdot(sigma = sktd.cv2d$Best_sigma),
                           lam1 = sktd.cv2d$Best_lambda,
                           sparsity = TRUE,
                           lam2 = lam2.seq[i],
                           ftol = 1e-3,
                           partol = 1e-3,
                           innerpartol = 1e-5)
  wt.tmp <- sktd.tmp$estimates[[1]]$weight
  if (is.null(wt.tmp)) wts[i, ] <- 0 else wts[i, ] <- wt.tmp
}
# plot the solution path with graphics::matplot()
matplot(y = wts,
        x = lam2.seq,
        type = "l",
        log = "x",
        ylab = "Weights",
        xlab = expression(paste(lambda)),
        lwd = 2)
legend("topright",
       title = "w index",
       legend = 1:5,
       lty = 1:5,
       col = 1:6,
       lwd = 2)

