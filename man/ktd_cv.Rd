% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ktd_cv.R
\name{ktd_cv}
\alias{ktd_cv}
\title{Cross validation for tuning the regularization coefficient in the kernel Tweedie model}
\usage{
ktd_cv(x, y, kern, lambda, nfolds = 5, rho = 1.5, loss = "LL", ...)
}
\arguments{
\item{x}{Covariate matrix.}

\item{y}{Outcome vector (e.g. insurance cost).}

\item{kern}{Choice of kernel. See \code{\link{dots}} for details on supported kernel functions.}

\item{lambda}{A vector of candidate regularization coefficients used in cross-validation.}

\item{nfolds}{Number of folds in cross-validation. Default is 5.}

\item{rho}{The power parameter of the Tweedie model. Default is 1.5 and can take any real value between 1 and 2.}

\item{loss}{Criterion used in cross-validation. "LL" for log likelihood, "RMSE" for root mean squared error, "MAD" for mean absolute difference. Default is "LL".}

\item{...}{Optional arguments to be passed to \code{\link{ktd_estimate}()}.}
}
\value{
A list of two items.
\enumerate{
\item LL or RMSE or MAD: a vector of validation error based on the user-specified \code{loss}, named by the corresponding \code{lambda} values;
\item Best_lambda: the \code{lambda} value in the pair that generates the best loss;
}
}
\description{
\code{ktd_cv()} performs cross-validation to determine the optimal regularization coefficient of the \code{ktweedie} model.
}
\details{
\code{ktd_cv()} is a built-in wrapper for cross-validation for the choice of regularization coefficient.
}
\examples{
# Provide a sequence of candidate values to the argument lambda.
# ktd_cv() will perform cross-validation to determine which is the best.
( cv1d <- ktd_cv(x = dat$x, y = dat$y,
                 kern = rbfdot(sigma = 1e-8),
                 lambda = 10^(-8:-1),
                 nfolds = 5) )
}
\seealso{
\code{\link{ktd_cv2d}}, \code{\link{ktd_estimate}}, \code{\link{ktd_predict}}
}
