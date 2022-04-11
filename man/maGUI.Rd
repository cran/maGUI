\name{maGUI}
\alias{maGUI}
\title{
A Graphical User Interface for Microarray Data Analysis and Annotation
}
\description{
Provides a graphical user interface for pre-processing and analyses of DNA microarray data of various platforms. It performs functional annotation and pathway analysis of differentially expressed genes. It identifies gene symbols, performs sample size estimation and also builds co-expression network.
}
\details{
Loading of Affymetrix data requires selecting the folder with raw files of .CEL extension, Agilent-OneColor with raw files containing foreground mean signal and background median signal values, Agilent-TwoColor with Agilent source and Nimblegen with .xys or .pair extensions. Loading of Illumina Beadarray and Lumi data require selecting the files with non-normalized data while Series Matrix requires selecting platform and series matrix files. Online data can be loaded by providing a GSE number and a live connection. If the loaded microarray data is two group, it can be pre-processed and analyzed automatically using the default parameters and unspecific filtering. Independent steps with different parameters are also available for pre-processing and analyzing any group microarray data. A user can perform miscellaneous applications such as gene set enrichment analysis, test analysis, pathway, gene symbol and co-expression network. For functional annotations and pathways, the package scans for the local GEOmetadb.sqlite, otherwise, it downloads and uses the latest version from the online source. If it fails in either case, it uses the internal data frame for functional annotations and pathways. Additionally, it uses the internal data frame of Bioconductor software, KEGG organism names and codes for pathway analyses. The annotated data can be visualized in the form of graphs which highlights enriched GO terms of differentially expressed genes in yellow color. The pathways can be visualized in the form of graphs which highlights differentially expressed genes in red and green color palette. The package has inbuilt functionalities to visualize tables and graphs in the graphical user interface. Further, a user can add, remove or delete child items from tree widget to observe distinguishable results. At any point, the required namespaces are automatically loaded if the live connection is available. 
}
\author{
Dhammapal Bharne
}
