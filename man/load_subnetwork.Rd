% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_subnetwork.R
\name{load_subnetwork}
\alias{load_subnetwork}
\title{Loads the subnetwork information from pairs of drugs and diseases}
\usage{
load_subnetwork(DISEASE1, DISEASE2, DRUG1, DRUG2)
}
\arguments{
\item{DISEASE1}{name of the first disease. Should be compatible with
DisNor (https://disnor.uniroma2.it/) disease names.}

\item{DISEASE2}{name of the second disease. Should be compatible with
DisNor (https://disnor.uniroma2.it/) disease names.}

\item{DRUG1}{name of the first drug. Should be compatible with Drugbank
(https://go.drugbank.com/) drug names.}

\item{DRUG2}{name of the second drug. Should be compatible with Drugbank
(https://go.drugbank.com/) drug names.}
}
\value{
list
}
\description{
Loads the subnetwork information from pairs of drugs and diseases
}
\examples{
\dontrun{
library(drugdealR)
load_subnetwork('Myeloid Leukemia, Chronic', 'Multiple Sclerosis', 'Imatinib', 'Tandutinib')
}
}
