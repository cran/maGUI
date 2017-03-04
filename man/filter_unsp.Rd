\name{filter_unsp}
\alias{filter_unsp}
\title{
UnSpecific Filtering on Normalized Microarray Data
}
\description{
Performs UnSpecific filtering on normalized microarray data. It includes Standard Deviation and Expression filtering.
}
\details{
Function can perform UnSpecific filtering on normalized microarray data if exists object of Affymetrix data as dat2Affy.m, Agilent_OneColor as datAgOne2.m, Agilent_TwoColor as datAgTwo2.m, Illumina_Beadarray as datIllBA2.m2, Illumina_Lumi as lumi_NQ.m, Nimblegen as data.matrix_Nimblegen2.m, Series_Matrix as data.matrixNorm.m and Online_Data as data.matrix_onlineNorm.m. The outputs are as .f to the normalized data. It also generates data with p-value less than 0.01 as statistically significant genes.
}
\author{
Dhammapal Bharne
}
