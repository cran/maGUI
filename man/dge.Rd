\name{dge}
\alias{dge}
\title{
Differntial Gene Expression of Filtered Data
}
\description{
Performs differential gene expression of filtered microarray data. Parameters such as number of genes to list, the amount of log fold change, p-value cut off, adjustment method and sort order can be set from the dialog box.
}
\details{
Function can perform differential gene expression on filtered microarray data if exists object of Affymetrix data as dat2Affy.f, Agilent_OneColor as datAgOne2.f, Agilent_TwoColor as datAgTwo2.f, Illumina_Beadarray as datIllBA2.f, Illumina_Lumi as lumi_NQ.f, Nimblegen as data.matrix_Nimblegen2.f, Series_Matrix as data.matrixNorm.f and Online_Data as data.matrix_onlineNorm.f. The outputs are DE_Affy for Affymetrix, DE_Ag1 for Agilent_OneColor, DE_Ag2 for Agilent_TwoColor, DE_Il_B for Illumina_Beadarray, DE_Il_L for Illumina_Lumi, DE_N for Nimblegen, DE_S for Series_Matrix and DE_O for Online procedures.
}
\author{
Dhammapal Bharne
}
