% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_kernel.R
\name{kernel_rational}
\alias{kernel_rational}
\title{Generating A Single Point-wise Function Using Rational Quadratic}
\usage{
kernel_rational(Sigma, l, p)
}
\arguments{
\item{Sigma}{(matrix) The covariance matrix for neural network kernel.}

\item{l}{(numeric) A numeric number indicating the hyperparameter 
(flexibility) of a specific kernel.}

\item{p}{(integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
rational, alpha = p.}
}
\value{
\item{point_wise}{(function) A function calculating 
the relevance of two vectors.}
}
\description{
Generate point-wise functions for two vectors using rational kernel.
}
\details{
\bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
l^2}\Big)^{-\alpha}}
}
\references{
The MIT Press. Gaussian Processes for Machine Learning, 2006.
}
\author{
Wenying Deng
}
