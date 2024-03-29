% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustAR.R
\name{robustAR.conf}
\alias{robustAR.conf}
\title{Robust AR Test confidence set.}
\usage{
robustAR.conf(
  Y,
  X,
  Z,
  W,
  betagrid,
  alpha,
  type = c("Huber", "Tukey", "OLS"),
  weighting = c("hat", "mcd", "no"),
  homoskedasticity = "no"
)
}
\arguments{
\item{Y}{A vector of dependent variables.}

\item{X}{A vector of endogenous covariates.}

\item{Z}{A matrix/vector of instrumental variables.}

\item{W}{A matrix/vector of control variables.}

\item{betagrid}{A grid of beta values}

\item{alpha}{A scalar significance level between 0 and 1}

\item{type}{Specifies which downweighting function to use: Huber, Tukey or OLS.}

\item{weighting}{Specifies which weighting method to use: mcd, hat matrix or no weights.}

\item{homoskedasticity}{Assume homoskedasticity: yes or no.}
}
\value{
A list with a confidence set and corresponding p-values.
}
\description{
Automatically computes a confidence set for the robust AR test.
}
