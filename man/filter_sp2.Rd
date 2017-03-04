\name{filter_sp2}
\alias{filter_sp2}
\title{
Specific Filtering on Normalized Microarray Data
}
\description{
Performs specific filtering on normalized microarray data. It takes control and test sample names of one or mulitple groups as input.
}
\details{
Function can perform specific filtering on normalized microarray data if exists object new of Affymetrix data as use.dat2Affy.m, Agilent_OneColor as use.datAgOne2.m, Agilent_TwoColor as use.datAgTwo2.m, Illumina_Beadarray as use.datIllBA2.m2, Illumina_Lumi as use.lumi_NQ.m, Nimblegen as use.data.matrix_Nimblegen2.m, Series_Matrix as use.data.matrixNorm.m and Online_Data as use.data.matrix_onlineNorm.m. The outputs are as .f object to the normalized data. It also generates data with p-value less than 0.01 as statistically significant genes.
}
\author{
Dhammapal Bharne
}
