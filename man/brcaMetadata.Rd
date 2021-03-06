% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{brcaMetadata}
\alias{brcaMetadata}
\title{A data.frame with patient metadata for breast 
cancer patients.}
\format{
A data.frame object
}
\usage{
data(brcaMetadata)
}
\description{
Has metadata for patients for which DNA methylation
or chromatin accessibility data was included as package data (329 patients).
Rows are patients,
with TCGA patient identifiers as row names and the column "subject_ID". 
Also includes columns: ER_status, ER_percent, age_diagnosis, days_to_death,
and days_to_last_follow_up.
Metadata is from The Cancer Genome Atlas
(TCGA-BRCA, https://portal.gdc.cancer.gov/).
}
\keyword{datasets}
