dmsm<-function(h,...){

	gconfirm("Select Platform Soft file",title="Loading",icon="info")
#	choose_file()
#	filechoose=NULL;
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
#	filechoose<<-file
#	file_soft<<-filechoose
	file_soft<<-file
#	filechoose=NULL
	if(file_soft==""){
		file_soft<<-NULL
		} else { try(({
	folder_S<<-dirname(file_soft)
	setwd(folder_S)	
	svalue(sb)<-"				Please wait while loading.."
	a2=readLines(file_soft)
	b2=length(a2)
	begin2<-which(gsub("\t","",a2)=="!platform_table_begin")
	end2<-which(gsub("\t","",a2)=="!platform_table_end")-1
	size2=end2-begin2-1
	new_y<<-NULL
	new_y<<-read.table(file_soft,skip=begin2,nrows=size2,header=TRUE,sep="\t",row.names=NULL)
	if(length(new_y)!=0){y<<-new_y}
	if(length(new_y)==0)gmessage("Could not load this type of file..!",title="Loading error",icon="error")
	}),silent=TRUE)
	}
	
	gconfirm("Select Series Matrix file",title="Loading",icon="info")
#	choose_file()
#	filechoose=NULL;
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
#	filechoose<<-file
#	file_series_mat<<-filechoose
	file_series_mat<<-file
#	filechoose=NULL
	if(file_series_mat==""){
		file_series_mat<<-NULL
		} else { try(({
	directory1<-dirname(file_series_mat)
	setwd(directory1)
	new_data.matrix<<-NULL;
	svalue(sb)<-"				Please wait while loading.."
	a=readLines(file_series_mat)
	b=length(a)
	begin<-which(gsub("\t","",a)=="!series_matrix_table_begin")
	end<-which(gsub("\t","",a)=="series_matrix_ends!")-1
	size=end-begin-1
	x=read.table(file_series_mat,skip=begin,nrows=size,header=TRUE,sep="\t",row.names=NULL)
	row.names(x)=x[,1]
	x=x[,-1]
	new_gse<-NULL
	new_gse<-new("ExpressionSet",exprs=as.matrix(x))
	if(length(new_gse)!=0){
	gse<-new_gse
	ann_S<<-annotation(gse)
	}
	table=exprs(gse)
	orf_ids<-cbind(y$ID,as.character(y$ORF))
	table2=table
	for(i in 1:length(orf_ids[,1]))
	{
		for(j in 1:length(orf_ids[,1]))
		{
			if(rownames(table2)[i]==orf_ids[j,1])
			{
				rownames(table2)[i]<-orf_ids[j,2]
				break
			}
		}
	}	
	table3<-table2[-(which(row.names(table2)=="")),]
	new_data.matrix<<-table3
	data.matrix_series<<-table3
	if(length(new_data.matrix)!=0){data.matrix_series<<-new_data.matrix}
	if(length(new_data.matrix)==0)gmessage("Could not load this type of file..!",title="Loading error",icon="error")
	  }),silent=TRUE)
	}
	if(length(data.matrix_series)!=0){
		delete(g1,g1_1)
		l$Series_Matrix<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(data.matrix_series)!=0)
		{
		if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
		{
			svalue(sb)<-"				Please wait while log transforming.."
			Med<-median(data.matrix_series,na.rm=TRUE)
			if(Med>16)
			{
				data.matrixLog<<-log2(data.matrix_series)
				} else {
				data.matrixLog<<-data.matrix_series
				}
			rm(Med)
			svalue(sb)<-"Done"
	
			svalue(sb)<-"				Please wait while imputing.."
			na.length<<-length(which(is.na(data.matrixLog)==TRUE))
			if(na.length > 0)data.matrixImp<<-impute.knn(data.matrixLog,k=10,rowmax=0.5,colmax=0.3)$data
			if(na.length <= 0)data.matrixImp<<-data.matrixLog
			svalue(sb)<-"Done"

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
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.data.matrixNorm.m,fff)
		dat.fo<-use.data.matrixNorm.m[i,]
		i<-genefilter(-use.data.matrixNorm.m,fff)
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
		err<-try(ttx<<-topTable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m)),silent=TRUE)
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
	err<-try(DE_S<<-topTable(data.matrixNorm.f),silent=TRUE)
	rm(DE_S2)	
	err<-try(DE_S2<<-topTable(data.matrixNorm.f,number=nrow(data.matrixNorm.f)),silent=TRUE)
	rownames(DE_S)<<-DE_S[,1]
	DE_S<<-DE_S[,-1]
	try(rownames(DE_S2)<<-DE_S2[,1],silent=TRUE)
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
if (!requireNamespace("Biobase", quietly = TRUE))
{
	install.packages("Biocbase")
} else {
	library(Biobase)
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
dmsm()
