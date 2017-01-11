dmsm<-function(h,...){
	data.matrixNorm=NULL;
	data.matrixNorm.m=NULL;
	use.data.matrixNorm.m=NULL;
	xf=NULL;
	data.matrixNorm.f=NULL;
	design_S=NULL;
	ttx=NULL;
	data.matrixNorm.s=NULL;
	DE_S=NULL;
	DE_S2=NULL;
	pca_S=NULL;
	sample.dist_S=NULL;
	sample.clust_S=NULL;
	use_DE_S=NULL;
	Clas_S=NULL;
	design=NULL;

	try(({folder_S<-folder_S;data.matrixImp<-data.matrixImp;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_S)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	if(length(data.matrixImp)!=0){
	data.matrixNorm<-normalizeBetweenArrays(data.matrixImp,method="quantile")
	tmp<-aggregate(data.matrixNorm,list(rownames(data.matrixNorm)),median)
	rm(data.matrixNorm.m)
	data.matrixNorm.m<<-as.matrix(tmp[,-1])
	rownames(data.matrixNorm.m)<<-tmp[,1]
	try(if(length(rownames(data.matrixNorm.m))==0)rownames(data.matrixNorm.m)<<-1:length(data.matrixNorm.m[,1]),silent=TRUE)
	}
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(data.matrixNorm.m)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.data.matrixNorm.m)
	use.data.matrixNorm.m<<-data.matrixNorm.m
	try(({
		ff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.data.matrixNorm.m,ff)
		dat.fo<-use.data.matrixNorm.m[i,]
		i<-genefilter(-use.data.matrixNorm.m,ff)
		dat.fu<-use.data.matrixNorm.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(data.matrixNorm.f)
		data.matrixNorm.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(data.matrixNorm.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrixNorm.m)
			i<-rsd>=2
			dat.f<-use.data.matrixNorm.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(data.matrixNorm.f)
			data.matrixNorm.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(data.matrixNorm.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(data.matrixNorm.s)
		err=NULL;
		err<-try(data.matrixNorm.s<<-use.data.matrixNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			data.matrixNorm.s<<-use.data.matrixNorm.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_S)
	err<-try(DE_S<<-toptable(data.matrixNorm.f),silent=TRUE)	
	err<-try(DE_S2<<-toptable(data.matrixNorm.f,number=nrow(data.matrixNorm.f)),silent=TRUE)
	rownames(DE_S)<<-DE_S[,1]
	DE_S<<-DE_S[,-1]
	rownames(DE_S2)<<-DE_S2[,1]
	DE_S2<<-DE_S2[,-1]
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_S)
	pca_S<<-prcomp(t(data.matrixNorm.m))
	plot(pca_S,main="Series_Matrix PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrixNorm.m)) 
	rm(sample.dist_S)
	sample.dist_S<<-as.dist(1-cor(data.matrixNorm.m,method="pearson"))
	rm(sample.clust_S)
	sample.clust_S<<-hclust(sample.dist_S,method="complete")
	plot(sample.clust_S)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_S)
	er_x<-try(use_DE_S<<-as.matrix(DE_S),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_S)<<-as.character(DE_S)
			use_DE_S<<-as.matrix(DE_S)
	}
	rm(Clas_S)	
	err2<-try(Clas_S<<-use.data.matrixNorm.m[rownames(use_DE_S),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_S,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Series_Matrix$Normalization<<-list()
		l$Series_Matrix$QC_Plot<<-list()
		l$Series_Matrix$Filtered<<-list()
		l$Series_Matrix$Stat_Significant<<-list()
		l$Series_Matrix$DGE<<-list()
		l$Series_Matrix$PCA_Plot<<-list()
		l$Series_Matrix$Cluster_Plot<<-list()
		l$Series_Matrix	$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}

