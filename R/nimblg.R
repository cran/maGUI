nimblg<-function(h,...){
	data.matrix_Nimblegen2=NULL;
	data.matrix_Nimblegen2.m=NULL;
	use.data.matrix_Nimblegen2.m=NULL;
	xf=NULL;
	data.matrix_Nimblegen2.f=NULL;
	design_N=NULL;
	ttx=NULL;
	data.matrix_Nimblegen2.s=NULL;
	DE_N=NULL;
	DE_N2=NULL;
	pca_N=NULL;
	sample.dist_N=NULL;
	sample.clust_N=NULL;
	use_DE_N=NULL;
	Clas_N=NULL;
	design=NULL;

	try(({folder_N<-folder_N;data.matrix_Nimblegen<-data.matrix_Nimblegen;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_N)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	data.matrix_Nimblegen2<<-normalizeBetweenArrays(data.matrix_Nimblegen,method="quantile")
	rm(data.matrix_Nimblegen2.m)
	data.matrix_Nimblegen2.m<<-normalizeBetweenArrays(data.matrix_Nimblegen,method="quantile")
	data.matrix_Nimblegen2.m<<-as.matrix(data.matrix_Nimblegen2.m)
	rownames(data.matrix_Nimblegen2.m)<<-as.character(rownames(data.matrix_Nimblegen2.m))
	try(if(length(rownames(data.matrix_Nimblegen2.m))==0)rownames(data.matrix_Nimblegen2.m)<<-1:length(data.matrix_Nimblegen2.m[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(data.matrix_Nimblegen2.m)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.data.matrix_Nimblegen2.m)
	use.data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m
	try(({
		ff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.data.matrix_Nimblegen2.m,ff)
		dat.fo<-use.data.matrix_Nimblegen2.m[i,]
		i<-genefilter(-use.data.matrix_Nimblegen2.m,ff)
		dat.fu<-use.data.matrix_Nimblegen2.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(data.matrix_Nimblegen2.f)
		data.matrix_Nimblegen2.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(data.matrix_Nimblegen2.f)==0)
	{
		try(({
			rsd<-rowSds(use.data.matrix_Nimblegen2.m)
			i<-rsd>=2
			dat.f<-use.data.matrix_Nimblegen2.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(data.matrix_Nimblegen2.f)
			data.matrix_Nimblegen2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(data.matrix_Nimblegen2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(data.matrix_Nimblegen2.s)
		err=NULL;
		err<-try(data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_N)
	err<-try(DE_N<<-toptable(data.matrix_Nimblegen2.f),silent=TRUE)	
	err<-try(DE_N2<<-toptable(data.matrix_Nimblegen2.f,number=nrow(data.matrix_Nimblegen2.f)),silent=TRUE)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_N)
	pca_N<<-prcomp(t(data.matrix_Nimblegen2.m))
	plot(pca_N,main="Nimblegen PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(data.matrix_Nimblegen2.m)) 
	rm(sample.dist_N)
	sample.dist_N<<-as.dist(1-cor(data.matrix_Nimblegen2.m,method="pearson"))
	rm(sample.clust_N)
	sample.clust_N<<-hclust(sample.dist_N,method="complete")
	plot(sample.clust_N)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_N)
	er_x<-try(use_DE_N<<-as.matrix(DE_N),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_N)<<-as.character(DE_N)
			use_DE_N<<-as.matrix(DE_N)
	}
	rm(Clas_N)	
	err2<-try(Clas_N<<-use.data.matrix_Nimblegen2.m[rownames(use_DE_N),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_N,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Nimblegen$Normalization<<-list()
		l$Nimblegen$QC_Plot<<-list()
		l$Nimblegen$Filtered<<-list()
		l$Nimblegen$Stat_Significant<<-list()
		l$Nimblegen$DGE<<-list()
		l$Nimblegen$PCA_Plot<<-list()
		l$Nimblegen$Cluster_Plot<<-list()
		l$Nimblegen$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}

