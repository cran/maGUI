agone<-function(h,...){

#	choose_folder()
#	folderchoose=NULL;
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
#	folderchoose<<-folder
#	folder_Ag1<<-folderchoose
	folder_Ag1<<-folder
#	folderchoose=NULL
	new_datAgOne<<-NULL
	if(folder_Ag1==""){
		folder_Ag1<<-NULL
		} else { try(({
	setwd(folder_Ag1)
	svalue(sb)<-"				Please wait while loading.."
	new_datAgOne<<-read.maimages(files=dir(folder_Ag1),columns=list(G="gMeanSignal",Gb="gBGMedianSignal",R="gMeanSignal",Rb="gBGMedianSignal"))
	if(length(new_datAgOne)!=0){datAgOne<<-new_datAgOne}
	z<-as(datAgOne,"NChannelSet")
	ann_Ag1<<-annotation(z)
	 }),silent=TRUE)
	if(length(new_datAgOne)==0)gmessage("Could not load this type of data..!\nTry Series_matrix or Online_Data methods",title="Loading error",icon="error")
	}
	if(length(datAgOne)!=0){
		delete(g1,g1_1)
		l$Agilent_OneColor<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_datAgOne)!=0)
		{
			if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
			{
			
	datAgOne2=NULL;
	datAgOne2.m=NULL;
	use.datAgOne2.m=NULL;
	xf=NULL;
	datAgOne2.f=NULL;
	design_Ag1=NULL;
	ttx=NULL;
	datAgOne2.s=NULL;
	DE_Ag1=NULL;
	DE_Ag1_2=NULL;
	pca_Ag1=NULL;
	sample.dist_Ag1=NULL;
	sample.clust_Ag1=NULL;
	use_DE_Ag1=NULL;
	Clas_Ag1=NULL;
	design=NULL;

	try(({folder_Ag1<-folder_Ag1;datAgOne<-datAgOne;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Ag1)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	try(datAgOne2<<-backgroundCorrect(datAgOne,"normexp"),silent=TRUE)
	if(length(datAgOne2)==0)datAgOne2<-datAgOne
	datAgOne2<-normalizeBetweenArrays(datAgOne2$R,method="quantile")
	datAgOne2<-log(datAgOne2)
	rm(datAgOne2.m)
	datAgOne2.m<<-datAgOne2
	datAgOne2.m<<-as.matrix(datAgOne2.m)
	rownames(datAgOne2.m)<<-as.character(rownames(datAgOne2.m))
	try(if(length(rownames(datAgOne2.m))==0)rownames(datAgOne2.m)<<-1:length(datAgOne2.m[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(datAgOne2.m)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.datAgOne2.m)
	use.datAgOne2.m<<-datAgOne2.m
	try(({
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.datAgOne2.m,fff)
		dat.fo<-use.datAgOne2.m[i,]
		i<-genefilter(-use.datAgOne2.m,fff)
		dat.fu<-use.datAgOne2.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(datAgOne2.f)
		datAgOne2.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(datAgOne2.f)==0)
	{
		try(({
			rsd<-rowSds(use.datAgOne2.m)
			i<-rsd>=2
			dat.f<-use.datAgOne2.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(datAgOne2.f)
			datAgOne2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(datAgOne2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-topTable(datAgOne2.f,number=nrow(use.datAgOne2.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(datAgOne2.s)
		err=NULL;
		err<-try(datAgOne2.s<<-use.datAgOne2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			datAgOne2.s<<-use.datAgOne2.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_Ag1)
	err<-try(DE_Ag1<<-topTable(datAgOne2.f),silent=TRUE)	
	err<-try(DE_Ag1_2<<-topTable(datAgOne2.f,number=nrow(datAgOne2.f)),silent=TRUE)	
	
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_Ag1)
	pca_Ag1<<-prcomp(t(datAgOne2.m))
	plot(pca_Ag1,main="Agilent_OneColor PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(datAgOne2.m)) 
	rm(sample.dist_Ag1)
	sample.dist_Ag1<<-as.dist(1-cor(datAgOne2.m,method="pearson"))
	rm(sample.clust_Ag1)
	sample.clust_Ag1<<-hclust(sample.dist_Ag1,method="complete")
	plot(sample.clust_Ag1)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_Ag1)
	er_x<-try(use_DE_Ag1<<-as.matrix(DE_Ag1),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_Ag1)<<-as.character(DE_Ag1)
			use_DE_Ag1<<-as.matrix(DE_Ag1)
	}
	rm(Clas_Ag1)	
	err2<-try(Clas_Ag1<<-use.datAgOne2.m[rownames(use_DE_Ag1),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_Ag1,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Agilent_OneColor$Normalization<<-list()
		l$Agilent_OneColor$QC_Plot<<-list()
		l$Agilent_OneColor$Filtered<<-list()
		l$Agilent_OneColor$Stat_Significant<<-list()
		l$Agilent_OneColor$DGE<<-list()
		l$Agilent_OneColor$PCA_Plot<<-list()
		l$Agilent_OneColor$Cluster_Plot<<-list()
		l$Agilent_OneColor$Classification<<-list()
#		tr<<-gtree(offspring=tree,container=g1_1)
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
agone()
