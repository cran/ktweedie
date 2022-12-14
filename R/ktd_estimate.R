#' Estimate kernel Tweedie model coefficients
#'
#' \code{ktd_estimate()} estimates the coefficients of the kernel Tweedie model \code{ktweedie} and the sparse kernel Tweedie model \code{sktweedie}. The log of the expected Tweedie mean is modeled by a function in the reproducing kernel Hilbert space. The \code{sktweedie} has an integrated feature selection component that induces sparsity by applying weights on the features and penalizing the weights.
#'
#' @param x Covariate matrix.
#' @param y Outcome vector (e.g. insurance cost).
#' @param kern Choice of kernel. See \code{\link{dots}} for details on supported kernel functions.
#' @param lam1 A vector of regularization coefficients.
#' @param rho The power parameter of the Tweedie model. Default is 1.5 and can take any real value between 1 and 2.
#' @param ftol Stopping criterion based on objective function value. Default is 1e-8. See Details.
#' @param partol Stopping criterion based on the coefficient values. Default is 1e-8. See Details.
#' @param abstol Stopping criterion based on absolute value of the objective function. Default is 0.
#' @param maxit Maximum number of iterations.
#' @param sparsity Logical If true, the \code{sktweedie} model with variable selection will be used. Default is false, for the \code{ktweedie} model.
#' @param lam2 Regularization coefficient for the sparsity-inducing penalty in the \code{sktweedie} model.
#' @param innerpartol Stopping criterion for the inner loops that update kernel parameters and weights based on the coefficient values. See Details.
#' @param innermaxit Maximum number of iterations for the inner loops that update kernel parameters and variable weights. See Details.
#' @param verbose Logical indicating whether to show details of each update.
#'
#' @returns A list of three items.
#' 1. \code{estimates}: a list containing the final objective function values and kernel Tweedie regression coefficients for each \code{lam1}.
#' 2. \code{data}: stores the inputs, including the predictor matrix, the kernel function used in the fitting and \code{lam1}.
#' 3. \code{sparsity}: a logical variable indicating whether the \code{ktweedie} or \code{sktweedie} is fitted.
#' @details
#' \code{ktd_estimate()} stops when the absolute difference between the objective function values of the last two updates is smaller than \code{ftol}, or the sum of absolute differences between the coefficients of the last two updates is smaller than \code{partol}, or the objective function values is below \code{abstol}, before \code{maxit} is reached. For the \code{sktweedie} model, there are inner loops for the update of kernel regression coefficients and regularization weights. The \code{innerpartol} and \code{innermaxit} arguments are the counterparts of \code{partol} and \code{maxit} for the inner loops.
#'
#' @seealso \code{\link{ktd_cv}}, \code{\link{ktd_cv2d}}, \code{\link{ktd_predict}}, \code{\link{rbfdot}}
#'
#' @examples
#' ###### ktweedie ######
#' # Provide a sequence of candidate values to the argument lam1.
#' # Provide a kernel object to the argument kern.
#' lam1.seq <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1)
#' fit.ktd <- ktd_estimate(x = dat$x, y = dat$y,
#'                         kern = rbfdot(sigma = 1e-8),
#'                         lam1 = lam1.seq)
#' ###### sktweedie ######
#' # Set sparsity to TRUE and a lam2 to control the level of sparsity
#' # Decrease lam2 if "WARNING: All weights are zero..."
#' fit.sktd <- ktd_estimate(x = dat$x,
#'                          y = dat$y,
#'                          kern = rbfdot(sigma = 0.1),
#'                          lam1 = 5,
#'                          sparsity = TRUE,
#'                          lam2 = 1)
#' # variables with fitted weight equal to 0 are not selected
#'
#' @export
#'
#'
ktd_estimate <- function(x, y, kern, lam1, rho = 1.5,
                        ftol = 1e-8, partol = 1e-8,
                        abstol = 0, maxit = 1e6,
                        sparsity = FALSE, lam2 = 0,
                        innerpartol = 1e-6, innermaxit = 1e6,
                        verbose = FALSE) {

  if (rho <= 1 | rho >= 2) stop("'rho' must be in (1, 2).")

  x <- as.matrix(x)
  dims <- dim(x)
  N <- dims[1]
  p <- dims[2]
  nz <- length(y)
  if (N != nz) stop("Sample sizes are different in x and y.")

  # number of lam1's to be fitted
  nhyper <- length(lam1)
  estimates <- vector(mode = "list", length = nhyper)

  # sparse kernel feature
  if (sparsity) {

    # only rbf kernel supported in sparse kernel
    if (!is(kern, "rbfkernel")) stop("Only RBF kernel is supported in sparse kernel feature. 'kern' has to be a 'rbfkernel' class.")

    # number of sigmas
    sigma <- kern@kpar$sigma
    nsig <- length(sigma)

    # check length conformity

    if (nsig != nhyper) {
      if (nsig == 1L) {
        sigma <- rep(sigma, nhyper)
      } else {
        stop("Number of sigmas to evaluate is different from that of lam1.")
      }
    }

    nlam2 <- length(lam2)
    if (nlam2 != nhyper) {
      if (nlam2 == 1L) {
        lam2 <- rep(lam2, nhyper)
      } else {
        stop("Number of lam2 to evaluate is different from that of lam1.")
      }
    }

    fit <- .Fortran("td_sk",
                    as.integer(N),
                    as.integer(p),
                    as.double(x),
                    as.double(y),
                    as.double(lam1),
                    as.double(lam2),
                    as.double(sigma),
                    as.integer(nhyper),
                    as.double(rho),
                    as.double(ftol),
                    as.double(partol),
                    as.double(abstol),
                    as.double(0),
                    as.double(innerpartol),
                    as.integer(maxit),
                    as.integer(innermaxit),
                    as.logical(verbose),
                    fn_final = double(nhyper),
                    gradf_final = double(N*nhyper),
                    param_final = double(N*nhyper),
                    wt_final = double(p*nhyper),
                    convergence = integer(nhyper)
    )

    # Prepare output
    estimates.name <- character(nhyper)
    for (i in 1:nhyper) {
      if (fit$convergence[i] > 3) {
        estimates[[i]] <- list(convergence = fit$convergence[i])
      } else {
        estimates[[i]] <- list(fn = fit$fn_final[i],
                               # grad = matrix(fit$gradf_final[((i-1)*N+1):(i*N)], ncol = 1),
                               coefficient = matrix(fit$param_final[((i-1)*N+1):(i*N)], ncol = 1),
                               weight = matrix(fit$wt_final[((i-1)*p+1):(i*p)], ncol = 1),
                               convergence = fit$convergence[i]

        )
      }
      estimates.name[i] <- paste("l1",lam1[i], "l2", lam2[i], class(kern), sigma[i])
    }
    names(estimates) <- estimates.name
  }

  if (!sparsity) {
    if (length(kern@kpar[[1]]) != 1) stop("Please supply one set of parameters in the kernel function.")

    lam1 <- sort(lam1, decreasing = TRUE)
    lam2 <- rep(0, nhyper)
    K <- as.matrix(kernelMatrix(kernel = kern, x))

    fit <- .Fortran("td_bfgs",
                    as.integer(N),
                    as.double(K),
                    as.double(y),
                    as.double(lam1),
                    as.integer(nhyper),
                    as.logical(FALSE),
                    as.double(rep(0,N)),
                    as.double(rho),
                    as.double(ftol),
                    as.double(partol),
                    as.double(abstol),
                    as.integer(maxit),
                    as.logical(verbose),
                    resfn = double(nhyper),
                    resgradf = double(N*nhyper),
                    resparam = double(N*nhyper),
                    resKa = double(N*nhyper),
                    conv = integer(nhyper)
    )

    # Prepare output
    estimates.name <- character(nhyper)
    for (i in 1:nhyper) {
      if (fit$conv[i]>3) {
        estimates[[i]] <- list(convergence = fit$conv[i])
      } else {
        estimates[[i]] <- list(fn = fit$resfn[i],
                               coefficient = matrix(fit$resparam[((i-1)*N+1):(i*N)], ncol = 1),
                               convergence = fit$conv[i]
        )
      }
      estimates.name[i] <- paste("lambda", lam1[i])
    }
    names(estimates) <- estimates.name
  }

  data <- list(x = x,
               kernel = kern,
               lambda1 = lam1)

  if (sparsity) {
    data$lambda2 = lam2
    data$sigma = sigma
  }
  result <- list(estimates = estimates,
                 data = data,
                 sparsity = sparsity)
  return(result)
}



