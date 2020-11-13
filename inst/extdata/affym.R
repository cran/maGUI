affym<-function(h,...){

#	choose_folder()
#	folderchoose=NULL;
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
#	folderchoose<<-folder
#	folder_Affy<<-folderchoose
	folder_Affy<<-folder
#	folderchoose=NULL
	new_datAffy<<-NULL
	if(folder_Affy==""){
		folder_Affy<<-NULL
		} else { try(({
	setwd(folder_Affy)
	svalue(sb)<-"				Please wait while loading.."
	new_datAffy<<-ReadAffy()
	if(length(new_datAffy)!=0){datAffy<<-new_datAffy}
	ann_Affy<<-annotation(datAffy)
	 }),silent=TRUE)
	if(length(new_datAffy)==0)gmessage("Could not load this type of data..!\nTry Series_matrix or Online_Data methods",title="Loading error",icon="error")
	 }
	if(length(datAffy)!=0){
		delete(g1,g1_1)
		l$Affymetrix<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_datAffy)!=0)
		{
			if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
			{
						
	dat2Affy.m=NULL;
	aqc=NULL;
	use.dat2Affy.m=NULL;
	xf=NULL;
	dat2Affy.f=NULL;
	design_Affy=NULL;
	ttx=NULL;
	dat2Affy.s=NULL;
	DE_Affy=NULL;
	DE_Affy2=NULL;
	pca_Affy=NULL;
	sample.dist_Affy=NULL;
	sample.clust_Affy=NULL;
	use_DE_Affy=NULL;
	Clas_Affy=NULL;
	design=NULL;

	try(({folder_Affy<-folder_Affy;datAffy<-datAffy;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)
	setwd(folder_Affy)
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	dat2Affy<-justRMA()
	dat2Affy.m<-exprs(dat2Affy)
	dat2Affy.m<<-as.matrix(dat2Affy.m)
	rownames(dat2Affy.m)<<-as.character(rownames(dat2Affy.m))
	try(if(length(rownames(dat2Affy.m))==0)rownames(dat2Affy.m)<<-1:length(dat2Affy.m[,1]),silent=TRUE)

#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	rm(aqc)
	try(aqc<<-qc(datAffy),silent=TRUE)
		err<-try(plot(aqc),silent=TRUE)
		if(length(grep("Error in",err))!=0){
		plot(dat2Affy.m)
	}
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.dat2Affy.m)
	use.dat2Affy.m<<-dat2Affy.m
	try(({
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.dat2Affy.m,fff)
		dat.fo<-use.dat2Affy.m[i,]
		i<-genefilter(-use.dat2Affy.m,fff)
		dat.fu<-use.dat2Affy.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(dat2Affy.f)
		dat2Affy.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(dat2Affy.f)==0)
	{
		try(({
			rsd<-rowSds(use.dat2Affy.m)
			i<-rsd>=2
			dat.f<-use.dat2Affy.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(dat2Affy.f)
			dat2Affy.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(dat2Affy.f)!=0){
		rm(ttx)
		err<-try(ttx<<-topTable(dat2Affy.f,number=nrow(use.dat2Affy.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(dat2Affy.s)
		err=NULL;
		err<-try(dat2Affy.s<<-use.dat2Affy.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			dat2Affy.s<<-use.dat2Affy.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_Affy)
	err<-try(DE_Affy<<-topTable(dat2Affy.f),silent=TRUE)	
	err<-try(DE_Affy2<<-topTable(dat2Affy.f,number=nrow(dat2Affy.f)),silent=TRUE)	
	
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_Affy)
	pca_Affy<<-prcomp(t(dat2Affy.m))
	plot(pca_Affy,main="Affymetrix PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(dat2Affy.m)) 
	rm(sample.dist_Affy)
	sample.dist_Affy<<-as.dist(1-cor(dat2Affy.m,method="pearson"))
	rm(sample.clust_Affy)
	sample.clust_Affy<<-hclust(sample.dist_Affy,method="complete")
	plot(sample.clust_Affy)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_Affy)
	er_x<-try(use_DE_Affy<<-as.matrix(DE_Affy),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_Affy)<<-as.character(DE_Affy)
			use_DE_Affy<<-as.matrix(DE_Affy)
	}
	rm(Clas_Affy)	
	err2<-try(Clas_Affy<<-use.dat2Affy.m[rownames(use_DE_Affy),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_Affy,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Affymetrix$Normalization<<-list()
		l$Affymetrix$QC_Plot<<-list()
		l$Affymetrix$Filtered<<-list()
		l$Affymetrix$Stat_Significant<<-list()
		l$Affymetrix$DGE<<-list()
		l$Affymetrix$PCA_Plot<<-list()
		l$Affymetrix$Cluster_Plot<<-list()
		l$Affymetrix$Classification<<-list()
#		tr<-gtree(offspring=tree,container=g1_1)
#		size(tr)<-c(300,480)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
}
		} 
	svalue(sb)<-"Done"
	on.exit(setwd(cur_dir))
}

if (!requireNamespace("affy", quietly = TRUE))
{
	install.packages("affy")
} else {
	library(affy)
}
if (!requireNamespace("genefilter", quietly = TRUE))
{
	install.packages("genefilter")
} else {
	library(genefilter)
}
if (!requireNamespace("limma", quietly = TRUE))
{
	install.packages("limma")
} else {
	library(limma)
}
if (!requireNamespace("simpleaffy", quietly = TRUE))
{
	install.packages("simpleaffy")
} else {
	library(simpleaffy)
}
if (!requireNamespace("stats", quietly = TRUE))
{
	install.packages("stats")
} else {
	library(stats)
}
if (!requireNamespace("tcltk", quietly = TRUE))
{
	install.packages("tcltk")
} else {
	library(tcltk)
}
affym()
