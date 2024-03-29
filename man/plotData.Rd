% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotData.R
\name{plotData}
\alias{plotData}
\alias{plotData,RprimerProfile-method}
\alias{plotData,RprimerOligo-method}
\alias{plotData,RprimerAssay-method}
\alias{plotData,RprimerMatchOligo-method}
\alias{plotData,RprimerMatchAssay-method}
\title{Plot an Rprimer object}
\usage{
plotData(x, ...)

\S4method{plotData}{RprimerProfile}(x, type = "overview", highlight = NULL, rc = FALSE)

\S4method{plotData}{RprimerOligo}(x)

\S4method{plotData}{RprimerAssay}(x)

\S4method{plotData}{RprimerMatchOligo}(x)

\S4method{plotData}{RprimerMatchAssay}(x)
}
\arguments{
\item{x}{An \code{RprimerProfile}, \code{RprimerOligo} \code{RprimerAssay},
\code{rprimerMatchOligo} or \code{RprimerMatchAssay} object.}

\item{...}{Optional arguments for \code{RprimerProfile} objects.}

\item{type}{For \code{Rprimeroligo} objects:
Type of plot: \code{"overview"}, or
\code{"nucleotide"}, defaults to \code{"overview"}.}

\item{highlight}{For \code{Rprimeroligo} objects:
If a specific region within an overview plot should be highlighted.
A numeric vector indicating the start and end position,
e.g. \code{c(100, 1000)}, defaults to \code{NULL}
(i.e., no highlight).}

\item{rc}{For \code{Rprimeroligo} objects, and \code{type = "nucleotide"}:
If the plotted sequence should be displayed
as reverse complement or not.
\code{TRUE} or {FALSE}, defaults to \code{FALSE}.

See examples below.}
}
\value{
A plot.
}
\description{
\code{plotData} visualizes objects from all different
Rprimer classes.
}
\section{Methods (by class)}{
\itemize{
\item \code{plotData(RprimerProfile)}: 

\item \code{plotData(RprimerOligo)}: 

\item \code{plotData(RprimerAssay)}: 

\item \code{plotData(RprimerMatchOligo)}: 

\item \code{plotData(RprimerMatchAssay)}: 

}}
\examples{
#### Plot an RprimerProfile object

data("exampleRprimerProfile")
prof <- exampleRprimerProfile

## Plot an overview
plotData(prof, highlight = c(500, 1000))

## Select a region of interest
roi <- prof[prof$position >= 500 & prof$position <= 550, ]

## Plot an overview of the roi
plotData(roi)

## Plot the nucleotide distribution of the roi
plotData(roi, type = "nucleotide")

#### Plot an RprimerOligo object

data("exampleRprimerOligo")
plotData(exampleRprimerOligo)

#### Plot an RprimerAssay object

data("exampleRprimerAssay")
plotData(exampleRprimerAssay)

#### Plot an RprimerMatchOligo object

data("exampleRprimerMatchOligo")
plotData(exampleRprimerMatchOligo)

#### Plot an RprimerMatchAssay object

data("exampleRprimerMatchAssay")
plotData(exampleRprimerMatchAssay)
}
