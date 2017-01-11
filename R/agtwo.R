agtwo<-function(h,...){
	datAgTwo2=NULL;
	datAgTwo2.m=NULL;
	use.datAgTwo2.m=NULL;
	xf=NULL;
	datAgTwo2.f=NULL;
	design_Ag2=NULL;
	ttx=NULL;
	datAgTwo2.s=NULL;
	DE_Ag2=NULL;
	DE_Ag2_2=NULL;
	pca_Ag2=NULL;
	sample.dist_Ag2=NULL;
	sample.clust_Ag2=NULL;
	use_DE_Ag2=NULL;
	Clas_Ag2=NULL;
	design=NULL;

	try(({folder_Ag2<-folder_Ag2;datAgTwo<-datAgTwo;heatcol<-heatcol;l<-l;tree<-tree;
	}),silent=TRUE)

	setwd(folder_Ag2)	
#	galert("				Please wait while normalizing",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while normalizing.."
#	Sys.sleep(1)
	try(datAgTwo2<<-backgroundCorrect(datAgTwo,"normexp"),silent=TRUE)
	if(length(datAgTwo2)==0)datAgTwo2<-datAgTwo
	datAgTwo2<-normalizeBetweenArrays(datAgTwo2$R,method="quantile")
	datAgTwo2<-log(datAgTwo2)
	rm(datAgTwo2.m)
	datAgTwo2.m<<-datAgTwo2
	datAgTwo2.m<<-as.matrix(datAgTwo2.m)
	rownames(datAgTwo2.m)<<-as.character(rownames(datAgTwo2.m))
	try(if(length(rownames(datAgTwo2.m))==0)rownames(datAgTwo2.m)<<-1:length(datAgTwo2.m[,1]),silent=TRUE)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while QC check",title="Pre-processing",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while QC check.."
#	Sys.sleep(1)
	boxplot(datAgTwo2.m)
#	galert("Done",title="Pre-processing",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Filtering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Filtering.."
#	Sys.sleep(1)
	rm(use.datAgTwo2.m)
	use.datAgTwo2.m<<-datAgTwo2.m
	try(({
		ff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.datAgTwo2.m,ff)
		dat.fo<-use.datAgTwo2.m[i,]
		i<-genefilter(-use.datAgTwo2.m,ff)
		dat.fu<-use.datAgTwo2.m[i,]
		dat.f<-rbind(dat.fo,dat.fu)
		fit<-lmFit(dat.f)
		fit2<-eBayes(fit)
		rm(datAgTwo2.f)
		datAgTwo2.f<<-fit2
#		print("Expression Filtering")
	}),silent=TRUE)
	if(length(datAgTwo2.f)==0)
	{
		try(({
			rsd<-rowSds(use.datAgTwo2.m)
			i<-rsd>=2
			dat.f<-use.datAgTwo2.m[i,]
			fit<-lmFit(dat.f)
			fit2<-eBayes(fit)
			rm(datAgTwo2.f)
			datAgTwo2.f<<-fit2
#			print("Standard Deviation Filtering")
		}),silent=TRUE)
	}

	if(length(datAgTwo2.f)!=0){
		rm(ttx)
		err<-try(ttx<<-toptable(datAgTwo2.f,number=nrow(use.datAgTwo2.m)),silent=TRUE)
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		rm(datAgTwo2.s)
		err=NULL;
		err<-try(datAgTwo2.s<<-use.datAgTwo2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rn<-as.numeric(rn)
			datAgTwo2.s<<-use.datAgTwo2.m[rn,]
		}
	}
	rm(err)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"
	
#	galert("				Please wait while DGE",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while DGE.."
#	Sys.sleep(1)
	rm(DE_Ag2)
	err<-try(DE_Ag2<<-toptable(datAgTwo2.f),silent=TRUE)	
	err<-try(DE_Ag2_2<<-toptable(datAgTwo2.f,number=nrow(datAgTwo2.f)),silent=TRUE)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while PCA",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while PCA.."
#	Sys.sleep(1)
	rm(pca_Ag2)
	pca_Ag2<<-prcomp(t(datAgTwo2.m))
	plot(pca_Ag2,main="Agilent_TwoColor PCA")
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Clustering",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Clustering.."
#	Sys.sleep(1)

	sample.dist_c1<-dist(t(datAgTwo2.m)) 
	rm(sample.dist_Ag2)
	sample.dist_Ag2<<-as.dist(1-cor(datAgTwo2.m,method="pearson"))
	rm(sample.clust_Ag2)
	sample.clust_Ag2<<-hclust(sample.dist_Ag2,method="complete")
	plot(sample.clust_Ag2)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

#	galert("				Please wait while Classification",title="Analysis",delay=5,parent=c(600,400))
	svalue(sb)<-"				Please wait while Classification.."
#	Sys.sleep(1)

	rm(use_DE_Ag2)
	er_x<-try(use_DE_Ag2<<-as.matrix(DE_Ag2),silent=TRUE)
	if(length(grep("Error",er_x))!=0)
	{
			rownames(DE_Ag2)<<-as.character(DE_Ag2)
			use_DE_Ag2<<-as.matrix(DE_Ag2)
	}
	rm(Clas_Ag2)	
	err2<-try(Clas_Ag2<<-use.datAgTwo2.m[rownames(use_DE_Ag2),],silent=TRUE)
	#heatcol<<-colorRampPalette(c("Green","Red"))(32) ##############
	heatmap(Clas_Ag2,margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
#	galert("Done",title="Analysis",delay=3,parent=c(600,400))
	svalue(sb)<-"Done"

	try(
	({
		visible(g1_1)<-FALSE
		l$Agilent_TwoColor$Normalization<<-list()
		l$Agilent_TwoColor$QC_Plot<<-list()
		l$Agilent_TwoColor$Filtered<<-list()
		l$Agilent_TwoColor$Stat_Significant<<-list()
		l$Agilent_TwoColor$DGE<<-list()
		l$Agilent_TwoColor$PCA_Plot<<-list()
		l$Agilent_TwoColor$Cluster_Plot<<-list()
		l$Agilent_TwoColor$Classification<<-list()
		tr<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,400)
		visible(g1_1)<-TRUE
		display()
		}),silent=TRUE)
}
