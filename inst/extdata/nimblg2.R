nimblg<-function(h,...){
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
	tnbl_i<<-tnbl_i+1
	folder_N[[tnbl_i]]<<-folder
	new_data.matrix_Nimblegen<<-NULL
	if(folder_N[[tnbl_i]]==""){
		folder_N[[tnbl_i]]<-NULL
	} else {
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
			setwd(folder_N[[tnbl_i]])
			pair2xys(list.files(pattern='\\.pair$'))
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			ndf<-list.files(pattern=".ndf")
			seed<-new("NgsExpressionPDInfoPkgSeed",ndfFile=ndf,xysFile=list.files(pattern='xys$')[1])
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			ndf2=gsub(".ndf","",ndf)
			ndf3=gsub("_",".",ndf2)
			pkg=paste("pd",ndf3,sep=".")
			pkg2<-tolower(pkg)
			if(length(dir(pkg2))==0)makePdInfoPackage(seed,destDir=".")	
			install.packages(pkg2,type='source',repos=NULL)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
			pkg3<-pkg2
			xys=list.xysfiles()
			new_rawData<-read.xysfiles(xys,pkgname=pkg3)
			if(length(new_rawData)!=0){rawData=NULL;rawData<-new_rawData}
			z<-annotation(rawData)
			z1<-strsplit(z,"\\.")
			zz<-z1[[1]][2]
			ann_N[[tnbl_i]]<<-toupper(zz)
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			try(
			({
				res<-rma(rawData)
				new_data.matrix_Nimblegen<-exprs(res)
			}),silent=TRUE)
			if(length(new_data.matrix_Nimblegen)==0)
			{
				new_data.matrix_Nimblegen<-exprs(rawData)
			}
			if(length(new_data.matrix_Nimblegen)!=0){data.matrix_Nimblegen=NULL;data.matrix_Nimblegen<-new_data.matrix_Nimblegen}
		 }),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_data.matrix_Nimblegen)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	}
	if(length(data.matrix_Nimblegen)!=0){
		if(length(new_data.matrix_Nimblegen)!=0){

		prnt<-paste("Nimblegen",tnbl_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_data.matrix_Nimblegen[[tnbl_i]]<<-data.matrix_Nimblegen
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				data.matrix_Nimblegen2.m=NULL;use.data.matrix_Nimblegen2.m=NULL;xf=NULL;data.matrix_Nimblegen2.f=NULL;data.matrix_Nimblegen2.s=NULL;
				design=NULL;design_N=NULL;ttx=NULL;DE_N=NULL;DE_N_2=NULL;use_DE_N=NULL;
				sample.dist_N=NULL;
#				pca_N=NULL;
#				sample.clust_N=NULL;
#				clas_N=NULL;
				try(({folder_N[[tnbl_i]]<-folder_N[[tnbl_i]];data.matrix_Nimblegen<-data.matrix_Nimblegen;heatcol<-heatcol; }),silent=TRUE)
				setwd(folder_N[[tnbl_i]])	
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				data.matrix_Nimblegen2<-normalizeBetweenArrays(data.matrix_Nimblegen,method="quantile")
				data.matrix_Nimblegen2.m<-normalizeBetweenArrays(data.matrix_Nimblegen2,method="quantile")
				data.matrix_Nimblegen2.m<-as.matrix(data.matrix_Nimblegen2.m)
				try({
					if(length(rownames(data.matrix_Nimblegen2.m))==0){
						rownames(data.matrix_Nimblegen2.m)<-1:length(data.matrix_Nimblegen2.m[,1])
					} else {
						rownames(data.matrix_Nimblegen2.m)<-as.character(rownames(data.matrix_Nimblegen2.m))
					}
				},silent=TRUE)
				ta_nbl_norm<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.m)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.m)[2]){ 
					if(i==0){ ta_nbl_norm[[i,j]]<-colnames(data.matrix_Nimblegen2.m)[j] } else { 
					if(j==0){ ta_nbl_norm[[i,j]]<-rownames(data.matrix_Nimblegen2.m)[i] } else {
					tem<-data.matrix_Nimblegen2.m[i,j]
					ta_nbl_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }

#				data.matrix_Nimblegen2.m<<-data.matrix_Nimblegen2.m
				data.matrix_Nimblegen2.m2[[tnbl_i]]<<-data.matrix_Nimblegen2.m
				
				ta_nbl_norm2[[tnbl_i]]<<-ta_nbl_norm
#				

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(data.matrix_Nimblegen2.m)},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_Nimblegen2.m)},hscale=1.7,vscale=1.07)
				}
p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.data.matrix_Nimblegen2.m,fff)
					dat.fo<-use.data.matrix_Nimblegen2.m[i,]
					i<-genefilter(-use.data.matrix_Nimblegen2.m,fff)
					dat.fu<-use.data.matrix_Nimblegen2.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_N[[tnbl_i]]<<-fit2
					data.matrix_Nimblegen2.f<-fit2
				}),silent=TRUE)
				if(length(data.matrix_Nimblegen2.f)==0)
				{
					try(({
						rsd<-rowSds(use.data.matrix_Nimblegen2.m)
						i<-rsd>=2
						dat.f<-use.data.matrix_Nimblegen2.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_N[[tnbl_i]]<<-fit2
						data.matrix_Nimblegen2.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(data.matrix_Nimblegen2.f)	
				rownames(ps1)<-rownames(data.matrix_Nimblegen2.f)
				data.matrix_Nimblegen2.f2=NULL
				data.matrix_Nimblegen2.f2<-ps1
				ta_nbl_filt=NULL
				ta_nbl_filt<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.f2)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.f2)[2]){ 
					if(i==0){ ta_nbl_filt[[i,j]]<-colnames(data.matrix_Nimblegen2.f2)[j] } else { 
					if(j==0){ ta_nbl_filt[[i,j]]<-rownames(data.matrix_Nimblegen2.f2)[i] } else {
					tem<-data.matrix_Nimblegen2.f2[i,j]
					ta_nbl_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				data.matrix_Nimblegen2.f3[[tnbl_i]]<<-data.matrix_Nimblegen2.f2
				
				ta_nbl_filt2[[tnbl_i]]<<-ta_nbl_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(data.matrix_Nimblegen2.f)!=0){
					err<-try({
							ttx<-topTable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_Nimblegen2.m)==TRUE)]
							data.matrix_Nimblegen2.s<-use.data.matrix_Nimblegen2.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_Nimblegen2.m)==TRUE)])
						data.matrix_Nimblegen2.s<-data.matrix_Nimblegen2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				data.matrix_Nimblegen2.s2=NULL
				data.matrix_Nimblegen2.s2<-data.matrix_Nimblegen2.s
				ta_nbl_stat=NULL
				ta_nbl_stat<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.s2)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.s2)[2]){ 
					if(i==0){ ta_nbl_stat[[i,j]]<-colnames(data.matrix_Nimblegen2.s2)[j] } else { 
					if(j==0){ ta_nbl_stat[[i,j]]<-rownames(data.matrix_Nimblegen2.s2)[i] } else {
					tem<-data.matrix_Nimblegen2.s2[i,j]
					ta_nbl_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				data.matrix_Nimblegen2.s3[[tnbl_i]]<<-data.matrix_Nimblegen2.s2
				
				ta_nbl_stat2[[tnbl_i]]<<-ta_nbl_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_N<-topTable(data.matrix_Nimblegen2.f),silent=TRUE)	
				err<-try(DE_NAll[[tnbl_i]]<<-topTable(data.matrix_Nimblegen2.f,number=nrow(data.matrix_Nimblegen2.f)),silent=TRUE)
				ta_nbl_dge=NULL
				ta_nbl_dge<-tclArray()
				for(i in 0:dim(DE_N)[1]){ for(j in 0:dim(DE_N)[2]){ 
					if(i==0){ ta_nbl_dge[[i,j]]<-colnames(DE_N)[j] } else { 
					if(j==0){ ta_nbl_dge[[i,j]]<-rownames(DE_N)[i] } else {
					tem<-DE_N[i,j]
					ta_nbl_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_N_2[[tnbl_i]]<<-DE_N
				
				ta_nbl_dge2[[tnbl_i]]<<-ta_nbl_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_N[[tnbl_i]]<<-prcomp(t(data.matrix_Nimblegen2.m))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_N[[tnbl_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_N[[tnbl_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(data.matrix_Nimblegen2.m)) 
				sample.dist_N<-as.dist(1-cor(data.matrix_Nimblegen2.m,method="pearson"))
				sample.clust_N[[tnbl_i]]<<-hclust(sample.dist_N,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_N[[tnbl_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_N[[tnbl_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_N<-as.matrix(DE_N),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_N)<-as.character(DE_N)
					use_DE_N<-as.matrix(DE_N)
				}
				rnames<-rownames(use_DE_N)
				clas_N[[tnbl_i]]<<-as.matrix(data.matrix_Nimblegen2.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_N[[tnbl_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_N[[tnbl_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
				}
p_enter_class<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Classification"))
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
			}
		}
	}
	display()
	
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
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("oligo",quietly=TRUE)){BiocManager::install("oligo",ask=FALSE,update=FALSE);library(oligo)} else {library(oligo)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("pdInfoBuilder",quietly=TRUE)){BiocManager::install("pdInfoBuilder",ask=FALSE,update=FALSE);library(pdInfoBuilder)} else {library(pdInfoBuilder)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("genefilter",quietly=TRUE)){BiocManager::install("genefilter",ask=FALSE,update=FALSE);library(genefilter)} else {library(genefilter)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("limma",quietly=TRUE)){BiocManager::install("limma",ask=FALSE,update=FALSE);library(limma)} else {library(limma)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("stats",quietly=TRUE)){install.packages("stats",ask=FALSE,update=FALSE);library(stats)} else {library(stats)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk",ask=FALSE,update=FALSE);library(tcltk)} else {library(tcltk)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
nimblg()
