\name{maGUI}
\alias{maGUI}
\title{
A Graphical User Interface for Microarray Data Analysis and Annotation
}
\description{
Provides a comprehensive graphical user interface for analysis of Affymetrix, Agilent, Illumina, Nimblegen and other microarray data. It can perform miscellaneous tasks such gene set enrichment and test analysis, identifying gene symbols and building co-expression network. It can also estimate sample size for atleast two-fold expression change. Further, it can identify protein-protein associations from two normalized microarray data. The current version is its slenderized form for easy launching and implementation.
}
\details{
Loading of Affymetrix data requires raw files with .CEL extension, Agilent-OneColor requires raw files with foreground mean signal and background median signal values, Agilent-TwoColor requires raw files with Agilent source, Nimblegen require raw files with .xys or .pair extensions, Illumina Beadarray and Lumi data requires non-normalized files, Series Matrix requires platform file and series matrix file and Online data requires a GSE number. The loaded microarray data can pre-processed and analyzed. Miscellaneous applications such as gene set enrichment analysis, gene set test analysis, annotations and gene symbol can also be performed online or through GEO metadata SQLite database. Other miscellaneous applications include sample size estimation and co-regulation network. Protein-protein associations can also be predicted from two normalized microarray data.
}
\author{
Dhammapal Bharne
}
