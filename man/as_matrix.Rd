% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_matrix.R
\name{as_matrix}
\alias{as_matrix}
\title{Turn wide data to matrix}
\usage{
as_matrix(data, ...)
}
\arguments{
\item{data}{long data form}

\item{...}{extra args}
}
\value{
data in wide matrix form, with continuous seq of dates
}
\description{
Take long form input data and turn into wide form matrix with
continuous sequence of dates. This function assumes that users supply long
form count data where dates with 0 cases are explicit.
}
