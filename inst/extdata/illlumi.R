illlumi<-function(h,...){

#	choose_file()
#	filechoose=NULL;
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
#	filechoose<<-file
#	file<-filechoose
#	filechoose=NULL
	new_lumi_data<<-NULL
	if(file==""){
		file<<-NULL
		} else { try(({
	folder_Il_L<<-dirname(file)
	setwd(folder_Il_L)
	svalue(sb)<-"				Please wait while loading.."
	new_lumi_data<<-lumiR(file)
	summary(new_lumi_data,'QC')
	if(length(new_lumi_data)!=0){lumi_data<<-new_lumi_data}
	ann_Il_L<<-annotation(lumi_data)
	 }),silent=TRUE)
	if(length(new_lumi_data)==0)gmessage("Could not load this type of file..!\nTry Series_matrix or Online_Data methods",title="Loading error",icon="error")
	 }
	if(length(lumi_data)!=0){
		delete(g1,g1_1)
		l$Illumina_Lumi<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_lumi_data)!=0)
		{
			if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
			{
			
	lumi_NQ=NULL;
	lumi_NQ.m=NULL;
	use.lumi_NQ.m=NULL;
	xf=NULL;
	lumi_NQ.f=NULL;
	design_Il_L=NULL;
	ttx=NULL;
	lumi_NQ.s=NULL;
	DE_Il_L=NULL;
	DE_Il_L2=NULL;
	pca_Il_L=NULL;
	sample.dist_Il_L=NULL;
	sample.clust_Il_L=NULL;
	use_DE_Il_L=NULL;
	Clas_Il_L=NULL;
	design=NULL;

	try(({folder_Il_L<-folder_Il_L;lumi_data<-lumi_data;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Il_L)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
#	rm(lumi_NQ)
#	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
#	summary(lumi_NQ,'QC')
#	rm(lumi_NQ.m)
#	lumi_NQ.m<<-lumi_NQ
#	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	lumi_NQ<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	summary(lumi_NQ,'QC')
	rm(lumi_NQ.m)
	lumi_NQ.m<<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
	lumi_NQ.m<<-as.matrix(lumi_NQ.m)
	rownames(lumi_NQ.m)<<-as.character(rownames(lumi_NQ.m))
	try(if(length(rownames(lumi_NQ.m))==0)rownames(lumi_NQ.m)<<-1:length(lumi_NQ.m[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	plot(lumi_NQ.m,what="boxplot")
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.lumi_NQ.m)
	use.lumi_NQ.m<<-lumi_NQ.m
	try(({
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.lumi_NQ.m,fff)
		dat.fo<-use.lumi_NQ.m[i,]
		i<-genefilter(-use.lumi_NQ.m,fff)
		dat.fu<-use.lumi_NQ.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(lumi_NQ.f)
		lumi_NQ.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(lumi_NQ.f)==0)
	{
		try(({
			rsd<-rowSds(use.lumi_NQ.m)
			i<-rsd>=2
			dat.f<-use.lumi_NQ.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(lumi_NQ.f)
			lumi_NQ.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(lumi_NQ.f)!=0){
		rm(ttx)
		err<-try(ttx<<-topTable(lumi_NQ.f,number=nrow(use.lumi_NQ.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(lumi_NQ.s)
		err=NULL;
		err<-try(lumi_NQ.s<<-use.lumi_NQ.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			lumi_NQ.s<<-use.lumi_NQ.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_Il_L)
	err<-try(DE_Il_L<<-topTable(lumi_NQ.f),silent=TRUE)	
	err<-try(DE_Il_L2<<-topTable(lumi_NQ.f,number=nrow(lumi_NQ.f)),silent=TRUE)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_Il_L)
	pca_Il_L<<-prcomp(t(lumi_NQ.m))
	plot(pca_Il_L,main="Illumina_Lumi PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(lumi_NQ.m)) 
	rm(sample.dist_Il_L)
	sample.dist_Il_L<<-as.dist(1-cor(lumi_NQ.m,method="pearson"))
	rm(sample.clust_Il_L)
	sample.clust_Il_L<<-hclust(sample.dist_Il_L,method="complete")
	plot(sample.clust_Il_L)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_Il_L)
	er_x<-try(use_DE_Il_L<<-as.matrix(DE_Il_L),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_Il_L)<<-as.character(DE_Il_L)
			use_DE_Il_L<<-as.matrix(DE_Il_L)
	}
	rm(Clas_Il_L)	
	err2<-try(Clas_Il_L<<-use.lumi_NQ.m[rownames(use_DE_Il_L),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_Il_L,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Illumina_Lumi$Normalization<<-list()
		l$Illumina_Lumi$QC_Plot<<-list()
		l$Illumina_Lumi$Filtered<<-list()
		l$Illumina_Lumi$Stat_Significant<<-list()
		l$Illumina_Lumi$DGE<<-list()
		l$Illumina_Lumi$PCA_Plot<<-list()
		l$Illumina_Lumi$Cluster_Plot<<-list()
		l$Illumina_Lumi$Classification<<-list()
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
if (!requireNamespace("lumi", quietly = TRUE))
{
	install.packages("lumi")
} else {
	library(lumi)
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
illlumi()
