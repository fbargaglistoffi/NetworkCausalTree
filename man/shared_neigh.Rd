% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{shared_neigh}
\alias{shared_neigh}
\title{Number of Shared Neighbors}
\usage{
shared_neigh(i, j, Ne_list)
}
\arguments{
\item{i}{Unit ID}

\item{j}{Unit ID}

\item{Ne_list}{List of N elements - where N is the sample size -
where each element i contains the IDs of the direct neighbors of unit i}
}
\value{
A numeric value representing the number of shared neighbors between
unit i and j.
}
\description{
Computes the number of shared neighbors between unit i and unit j.
}
