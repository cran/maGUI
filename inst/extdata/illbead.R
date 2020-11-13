illbead<-function(h,...){

#	choose_file()
#	filechoose=NULL;
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
#	filechoose<<-file
#	file<-filechoose
#	filechoose=NULL
	new_datIllBA<<-NULL
	if(file==""){
		file<<-NULL
		} else { try(({
	folder_Il_B<<-dirname(file)
	setwd(folder_Il_B)
	svalue(sb)<-"				Please wait while loading.."
	probes<-c("ProbeID","TargetID","ProbeId","TargetId","Probe_ID","Probe_Id","Target_ID","Target_Id","ID_REF","Id_Ref")
	x<-read.table(file,sep="\t",fill=TRUE)
	y<-which(x[,4]!="")
	ToSkip<-y[1]-1
	x2<-read.table(file,sep="\t",header=TRUE,skip=ToSkip,fill=TRUE)
	ProbeID=probes[which(probes %in% colnames(x2))]
	new_datIllBA<<-readBeadSummaryData(file,sep="\t",skip=ToSkip,ProbeID=ProbeID)
	if(length(new_datIllBA)!=0){datIllBA<<-new_datIllBA}
	ann_Il_B<<-annotation(datIllBA)
	 }),silent=TRUE)
	if(length(new_datIllBA)==0)gmessage("Could not load this type of file..!\nTry Series_matrix or Online_Data methods",title="Loading error",icon="error")
	 }
	if(length(datIllBA)!=0){
		delete(g1,g1_1)
		l$Illumina_Beadarray<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_datIllBA)!=0)
		{
			if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
			{
			
	datIllBA2=NULL;
	datIllBA2.m2=NULL;
	use.datIllBA2.m2=NULL;
	xf=NULL;
	datIllBA2.f=NULL;
	design_Il_B=NULL;
	ttx=NULL;
	datIllBA2.s=NULL;
	DE_Il_B=NULL;
	DE_Il_B2=NULL;
	pca_Il_B=NULL;
	sample.dist_Il_B=NULL;
	sample.clust_Il_B=NULL;
	use_DE_Il_B=NULL;
	Clas_Il_B=NULL;
	design=NULL;

	try(({folder_Il_B<-folder_Il_B;datIllBA<-datIllBA;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Il_B)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
#	rm(datIllBA2)
	datIllBA2<-normaliseIllumina(datIllBA,method="quantile")
	datIllBA2.m<-exprs(datIllBA2)
	rm(datIllBA2.m2)
	datIllBA2.m2<<-log2(datIllBA2.m)
	datIllBA2.m2<<-as.matrix(datIllBA2.m2)
	rownames(datIllBA2.m2)<<-as.character(rownames(datIllBA2.m2))
	try(if(length(rownames(datIllBA2.m2))==0)rownames(datIllBA2.m2)<<-1:length(datIllBA2.m2[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(datIllBA2.m2)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.datIllBA2.m2)
	use.datIllBA2.m2<<-datIllBA2.m2
	try(({
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.datIllBA2.m2,fff)
		dat.fo<-use.datIllBA2.m2[i,]
		i<-genefilter(-use.datIllBA2.m2,fff)
		dat.fu<-use.datIllBA2.m2[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(datIllBA2.f)
		datIllBA2.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(datIllBA2.f)==0)
	{
		try(({
			rsd<-rowSds(use.datIllBA2.m2)
			i<-rsd>=2
			dat.f<-use.datIllBA2.m2[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(datIllBA2.f)
			datIllBA2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(datIllBA2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-topTable(datIllBA2.f,number=nrow(use.datIllBA2.m2)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(datIllBA2.s)
		err=NULL;
		err<-try(datIllBA2.s<<-use.datIllBA2.m2[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			datIllBA2.s<<-use.datIllBA2.m2[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_Il_B)
	err<-try(DE_Il_B<<-topTable(datIllBA2.f),silent=TRUE)
	err<-try(DE_Il_B2<<-topTable(datIllBA2.f,number=nrow(datIllBA2.f)),silent=TRUE)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_Il_B)
	pca_Il_B<<-prcomp(t(datIllBA2.m2))
	plot(pca_Il_B,main="Illumina_Beadarray PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(datIllBA2.m2)) 
	rm(sample.dist_Il_B)
	sample.dist_Il_B<<-as.dist(1-cor(datIllBA2.m2,method="pearson"))
	rm(sample.clust_Il_B)
	sample.clust_Il_B<<-hclust(sample.dist_Il_B,method="complete")
	plot(sample.clust_Il_B)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_Il_B)
	er_x<-try(use_DE_Il_B<<-as.matrix(DE_Il_B),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_Il_B)<<-as.character(DE_Il_B)
			use_DE_Il_B<<-as.matrix(DE_Il_B)
	}
	rm(Clas_Il_B)	
	err2<-try(Clas_Il_B<<-use.datIllBA2.m2[rownames(use_DE_Il_B),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_Il_B,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Illumina_Beadarray$Normalization<<-list()
		l$Illumina_Beadarray$QC_Plot<<-list()
		l$Illumina_Beadarray$Filtered<<-list()
		l$Illumina_Beadarray$Stat_Significant<<-list()
		l$Illumina_Beadarray$DGE<<-list()
		l$Illumina_Beadarray$PCA_Plot<<-list()
		l$Illumina_Beadarray$Cluster_Plot<<-list()
		l$Illumina_Beadarray$Classification<<-list()
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

if (!requireNamespace("beadarray", quietly = TRUE))
{
	install.packages("beadarray")
} else {
	library(beadarray)
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
illbead()
