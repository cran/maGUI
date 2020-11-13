dmo<-function(h,...){

	gse_no<<-ginput("Enter GSE Name",icon="question",handler=function(h,...)h$input)
	if(length(grep("GSE",gse_no))!=0)
	{
		svalue(sb)<-"				Please wait while loading.."
		new_data.matrix_online<<-NULL
		try(({
		err<-try(new_gse<-getGEO(gse_no,GSEMatrix=TRUE,destdir=cur_dir)[[1]],silent=TRUE)
		if(length(new_gse)!=0){
		gse<-new_gse
		ann_O<<-annotation(gse)
		}
		ORF<-featureData(gse)@data$ORF
		use_probe<-which(is.na(ORF)==F & match(ORF,"",nomatch=0)==0)
		new_data.matrix_online<<-exprs(gse)[use_probe,]
		data.matrix_online<<-exprs(gse)[use_probe,]
		rownames(new_data.matrix_online)<<-ORF[use_probe]
		if(length(new_data.matrix_online)!=0){data.matrix_online<<-new_data.matrix_online}
		rm(ORF,use_probe)
	 	}),silent=TRUE)
		if(length(new_data.matrix_online)==0)gmessage("Could not load this type of data..!",title="Loading error",icon="error")
		svalue(sb)<-"Done"
	 }
	if(length(data.matrix_online)!=0){
		delete(g1,g1_1)
		l$Online_Data<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_data.matrix_online)!=0)
		{
		if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
		{
			svalue(sb)<-"				Please wait while log transforming.."
			Med<-median(data.matrix_online,na.rm=TRUE)
			if(Med>16)
			{
				data.matrix_onlineLog<<-log2(data.matrix_online)
				} else {
					data.matrix_onlineLog<<-data.matrix_online
					}
			rm(Med)
			svalue(sb)<-"Done"
			
			svalue(sb)<-"				Please wait while imputing.."
			na.length<-length(which(is.na(data.matrix_onlineLog)==TRUE))
			if(na.length > 0)data.matrix_onlineImp<<-impute.knn(data.matrix_onlineLog,k=10,rowmax=0.5,colmax=0.3)$data
			if(na.length < 0)data.matrix_onlineImp<<-data.matrix_onlineLog
			svalue(sb)<-"Done"

	data.matrix_onlineNorm=NULL;
	data.matrix_onlineNorm.m=NULL;
	use.data.matrix_onlineNorm.m=NULL;
	xf=NULL;
	data.matrix_onlineNorm.f=NULL;
	design_O=NULL;
	ttx=NULL;
	data.matrix_onlineNorm.s=NULL;
	DE_O=NULL;
	DE_O2=NULL;
	pca_O=NULL;
	sample.dist_O=NULL;
	sample.clust_O=NULL;
	use_DE_O=NULL;
	Clas_O=NULL;
	design=NULL;

	try(({data.matrix_onlineImp<-data.matrix_onlineImp;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	data.matrix_onlineNorm<-normalizeBetweenArrays(data.matrix_onlineImp,method="quantile")
	tmp<-aggregate(data.matrix_onlineNorm,list(rownames(data.matrix_onlineNorm)),median)
	rm(data.matrix_onlineNorm.m)
	data.matrix_onlineNorm.m<<-as.matrix(tmp[,-1])
	rownames(data.matrix_onlineNorm.m)<<-tmp[,1]
	rownames(data.matrix_onlineNorm.m)<<-as.character(rownames(data.matrix_onlineNorm.m))
	try(if(length(rownames(data.matrix_onlineNorm.m))==0)rownames(data.matrix_onlineNorm.m)<<-1:length(data.matrix_onlineNorm.m[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(data.matrix_onlineNorm.m)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.data.matrix_onlineNorm.m)
	use.data.matrix_onlineNorm.m<<-data.matrix_onlineNorm.m
	try(({
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.data.matrix_onlineNorm.m,fff)
		dat.fo<-use.data.matrix_onlineNorm.m[i,]
		i<-genefilter(-use.data.matrix_onlineNorm.m,fff)
		dat.fu<-use.data.matrix_onlineNorm.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(data.matrix_onlineNorm.f)
		data.matrix_onlineNorm.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(data.matrix_onlineNorm.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrix_onlineNorm.m)
			i<-rsd>=2
			dat.f<-use.data.matrix_onlineNorm.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(data.matrix_onlineNorm.f)
			data.matrix_onlineNorm.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(data.matrix_onlineNorm.f)!=0){
		rm(ttx)
		err<-try(ttx<<-topTable(data.matrix_onlineNorm.f,number=nrow(use.data.matrix_onlineNorm.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(data.matrix_onlineNorm.s)
		err=NULL;
		err<-try(data.matrix_onlineNorm.s<<-use.data.matrix_onlineNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			data.matrix_onlineNorm.s<<-use.data.matrix_onlineNorm.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_O)
	err<-try(DE_O<<-topTable(data.matrix_onlineNorm.f),silent=TRUE)	
	rm(DE_O2)
	err<-try(DE_O2<<-topTable(data.matrix_onlineNorm.f,number=nrow(data.matrix_onlineNorm.f)),silent=TRUE)
	rownames(DE_O)<<-DE_O[,1]
	DE_O<<-DE_O[,-1]
	try(rownames(DE_O2)<<-DE_O2[,1],silent=TRUE)
	DE_O2<<-DE_O2[,-1]
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_O)
	pca_O<<-prcomp(t(data.matrix_onlineNorm.m))
	plot(pca_O,main="Online_Data PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrix_onlineNorm.m)) 
	rm(sample.dist_O)
	sample.dist_O<<-as.dist(1-cor(data.matrix_onlineNorm.m,method="pearson"))
	rm(sample.clust_O)
	sample.clust_O<<-hclust(sample.dist_O,method="complete")
	plot(sample.clust_O)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_O)
	er_x<-try(use_DE_O<<-as.matrix(DE_O),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_O)<<-as.character(DE_O)
			use_DE_O<<-as.matrix(DE_O)
	}
	rm(Clas_O)	
	err2<-try(Clas_O<<-use.data.matrix_onlineNorm.m[rownames(use_DE_O),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_O,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Online_Data$Normalization<<-list()
		l$Online_Data$QC_Plot<<-list()
		l$Online_Data$Filtered<<-list()
		l$Online_Data$Stat_Significant<<-list()
		l$Online_Data$DGE<<-list()
		l$Online_Data$PCA_Plot<<-list()
		l$Online_Data$Cluster_Plot<<-list()
		l$Online_Data	$Classification<<-list()
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
if (!requireNamespace("GEOquery", quietly = TRUE))
{
	install.packages("GEOquery")
} else {
	library(GEOquery)
}
if (!requireNamespace("impute", quietly = TRUE))
{
	install.packages("impute")
} else {
	library(impute)
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
dmo()
