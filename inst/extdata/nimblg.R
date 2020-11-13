nimblg<-function(h,...){

#	choose_folder()
#	folderchoose=NULL;
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
#	folderchoose<<-folder
#	folder_N<<-folderchoose
	folder_N<<-folder
#	folderchoose=NULL
	new_data.matrix_Nimblegen<<-NULL
	if(folder_N==""){
		folder_Ag2<<-NULL
		} else { try(({
	setwd(folder_N)
	svalue(sb)<-"				Please wait while loading.."
	pair2xys(list.files(pattern='\\.pair$'))
	ndf<-list.files(pattern=".ndf")
	seed<-new("NgsExpressionPDInfoPkgSeed",ndfFile=ndf,xysFile=list.files(pattern='xys$')[1])
	ndf2=gsub(".ndf","",ndf)
	ndf3=gsub("_",".",ndf2)
	pkg=paste("pd",ndf3,sep=".")
	pkg2<-tolower(pkg)
	if(length(dir(pkg2))==0)makePdInfoPackage(seed,destDir=".")	
	install.packages(pkg2,type='source',repos=NULL)
	pkg3<-pkg2
	xys=list.xysfiles()
	new_rawData<<-read.xysfiles(xys,pkgname=pkg3)
	if(length(new_rawData)!=0)rawData<<-new_rawData
	z<-annotation(rawData)
	z1<-strsplit(z,"\\.")
	zz<-z1[[1]][2]
	ann_N<<-toupper(zz)
	try(
	({
		res<-rma(rawData)
		new_data.matrix_Nimblegen<<-exprs(res)
	}),silent=TRUE)
	if(length(new_data.matrix_Nimblegen)==0)
	{
		new_data.matrix_Nimblegen<<-exprs(rawData)
	}
	if(length(new_data.matrix_Nimblegen)!=0){data.matrix_Nimblegen<<-new_data.matrix_Nimblegen}
	 }),silent=TRUE)
	if(length(new_data.matrix_Nimblegen)==0)gmessage("Could not load this type of data..!\nTry Series_matrix or Online_Data methods",title="Loading error",icon="error")
	 }
	if(length(data.matrix_Nimblegen)!=0){
		delete(g1,g1_1)
		l$Nimblegen<<-list()
		g1_1<<-ggroup(container=g1)
		tr<<-gtree(offspring=tree,container=g1_1)
		size(tr)<-c(300,480)
		if(length(new_data.matrix_Nimblegen)!=0)
		{
			if(gconfirm("Do you want to pre-process and analyze automatically",icon="question")==TRUE)
			{
			
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
		fff<-pOverA(A=1,p=0.05)
		i<-genefilter(use.data.matrix_Nimblegen2.m,fff)
		dat.fo<-use.data.matrix_Nimblegen2.m[i,]
		i<-genefilter(-use.data.matrix_Nimblegen2.m,fff)
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
		err<-try(ttx<<-topTable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m)),silent=TRUE)
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
	err<-try(DE_N<<-topTable(data.matrix_Nimblegen2.f),silent=TRUE)	
	err<-try(DE_N2<<-topTable(data.matrix_Nimblegen2.f,number=nrow(data.matrix_Nimblegen2.f)),silent=TRUE)
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

pair2xys<-function(pairFiles,outdir=getwd(),verbose=TRUE){
	if(verbose)message('Output directory:',outdir)
	for(pairFile in pairFiles)
	{
		if(verbose)message('Processing',basename(pairFile))
		header<-readLines(pairFile,n=1)
		pair<-read.delim(pairFile,header=TRUE,sep='\t',stringsAsFactors=FALSE,comment.char='#')
		maxX<-max(pair$X)
		maxY<-max(pair$Y)
		xys<-expand.grid(X=1:maxX,Y=1:maxY)
		xys<-merge(xys,pair[,c('X','Y','PM')],all.x=TRUE)
		xys<-pair[,c('X','Y','PM')]
		names(xys)<-c('X','Y','SIGNAL')
		xys$COUNT<-ifelse(is.na(xys$SIGNAL),NA_integer_,1L)
		xys<-xys[with(xys,order(Y,X)),]
		rownames(xys)<-NULL
		xysFile<-file.path(outdir,gsub('\\.pair$','\\.xys',basename(pairFile)))
		if(verbose)message('Writing',basename(xysFile))
		writeLines(header,con=xysFile)
		suppressWarnings(write.table(xys,file=xysFile,sep='\t',row.names=FALSE,quote=FALSE,append=TRUE))
	}
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
if (!requireNamespace("oligo", quietly = TRUE))
{
	install.packages("oligo")
} else {
	library(lumi)
}
if (!requireNamespace("pdInfoBuilder", quietly = TRUE))
{
	install.packages("pdInfoBuilder")
} else {
	library(pdInfoBuilder)
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
nimblg()
