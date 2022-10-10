\name{sleep}
\alias{sleep}
\docType{data}
\title{Infant sleep status data
}
\description{
The sleep data  consist of sleep state measurements of a newborn infant together with his heart rate and temperature sampled every 30 seconds. The sleep states are classified as: (1) quiet sleep, (2) indeterminate sleep, (3) active sleep, (4) awake. The total number of observations is equal to 1024 and the objective is to predict the sleep state based on covariate information.
}
\usage{data(sleep)}
\format{
  A data frame with 1024 observations on the following 3 variables:
  \describe{
    \item{\code{heartrate}}{Heart rate.}
    \item{\code{sleep}}{An ordinal time series in the sense that the response increases from awake to active sleep, i.e., (4) < (1) < (2) < (3).}
    \item{\code{temperature}}{Temperature}
    } 
}
\source{
Fokianos, K. and Kedem, B. (2003). Regression theory for categorical time series. \emph{Statistical Science}, \bold{18}(3):357--376.
}
\examples{
data(sleep)
}
\keyword{datasets}
