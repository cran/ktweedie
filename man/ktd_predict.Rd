% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ktd_predict.R
\name{ktd_predict}
\alias{ktd_predict}
\title{Predict outcome using fitted kernel Tweedie model}
\usage{
ktd_predict(model, newdata, which.lam1 = 1, type = "link")
}
\arguments{
\item{model}{Fitted model from \code{\link{ktd_estimate}}}

\item{newdata}{New x matrix for the prediction. If not provided, it will be the x matrix used to fit \code{model}.}

\item{which.lam1}{The index of the \code{lam1} in \code{model} used in the prediction. Default is 1.}

\item{type}{The type of prediction to be made - "\code{link}" for the linear predictor and "\code{response}" for the predicted outcome. Default is "\code{link}".}
}
\value{
A list named \code{prediction} containing the vector of predicted outcomes.
}
\description{
\code{ktd_predict()} predicts the outcome with fitted \code{ktweedie} or \code{sktweedie} model at the user supplied new data.
}
\details{
\code{ktd_predict()} uses the fitted model from \code{\link{ktd_estimate}} to estimate the mean outcome for new data points.
}
\examples{
# Fit a ktweedie model
fit <- ktd_estimate(x = dat$x, y = dat$y,
                    kern = rbfdot(sigma = 1e-6),
                    lam1 = 10^(-5:1))
# Generate newx at which predictions are to be made.
# The newdata should have the same dimension as the original trainig data.
newx <- matrix(rnorm(10 * ncol(dat$x)), nrow = 10)
pred <- ktd_predict(model = fit, newdata = newx,
                    which.lam1 = 3, type = "link")
}
\seealso{
\code{\link{ktd_estimate}}, \code{\link{ktd_cv}}, \code{\link{ktd_cv2d}}
}
