% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_N.R
\name{coupled_update_N}
\alias{coupled_update_N}
\title{coupled_update_N}
\usage{
coupled_update_N(state1, state2, V, hyper, algotuning)
}
\value{
a list with 'state1' and 'state2', the two new states of the Markov chains
}
\description{
samples from a maximal coupling of the conditional distributions of N given the other variables,
for the two chains given as 'state1' and 'state2'; see \code{\link{update_N}} for the other parameters.
}
