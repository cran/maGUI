dmsm<-function(h,...){
	tkmessageBox(title="Loading",message="Select Platform Soft file",icon="info",type="ok")
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
	tsmt_i<<-tsmt_i+1
	file_soft<<-file
	if(file_soft==""){
		file_soft<-NULL
	} else { 
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			folder_S[[tsmt_i]]<<-dirname(file_soft)
			setwd(folder_S[[tsmt_i]])	
			a2=readLines(file_soft)
			b2=length(a2)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			begin2<-which(gsub("\t","",a2)=="!platform_table_begin")
			end2<-which(gsub("\t","",a2)=="!platform_table_end")-1
			size2=end2-begin2-1
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			new_y=NULL
			new_y<-read.table(file_soft,skip=begin2,nrows=size2,header=TRUE,sep="\t",row.names=NULL)
			if(length(new_y)!=0){y=NULL;y<-new_y}
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
			if(length(new_y)==0)tkmessageBox(title="Error",message="Could not load this type of data!",icon="error",type="ok")
		}),silent=TRUE)
	}
	tkmessageBox(title="Loading",message="Select Series Matrix file",icon="info",type="ok")
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
	file_series_mat<-file
	if(file_series_mat==""){
		file_series_mat<-NULL
	} else { 
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
			directory1<-dirname(file_series_mat)
			setwd(directory1)
			new_data.matrix=NULL;
			data.matrix_series=NULL
			a=readLines(file_series_mat)
			b=length(a)
			begin<-which(gsub("\t","",a)=='!series_matrix_table_begin')
			end<-which(gsub("\t","",a)=='series_matrix_ends!')-1
			size=end-begin-1
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			x=read.table(file_series_mat,skip=begin,nrows=size,header=TRUE,sep="\t",row.names=NULL)
			row.names(x)=x[,1]
			x=x[,-1]
			new_gse<-NULL
			new_gse<-new("ExpressionSet",exprs=as.matrix(x))
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
			if(length(new_gse)!=0){
				gse<-new_gse
				ann_S[[tsmt_i]]<<-annotation(gse)
			}
			table=exprs(gse)
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			orf_ids<-cbind(y$ID,as.character(y$ORF))
			orf_ids2<-orf_ids[which(orf_ids[,2]!=""),]
			table2<-table[orf_ids2[,1],]
			rownames(table2)<-orf_ids2[,2]
			table3<-table2
			new_data.matrix<-table3
			
			data.matrix_series<-table3
			if(length(new_data.matrix)!=0){data.matrix_series=NULL;data.matrix_series<-new_data.matrix}
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
			if(length(new_data.matrix)==0)tkmessageBox(title="Error",message="Could not load this type of data!",icon="error",type="ok")
		}),silent=TRUE)
	}
	if(length(data.matrix_series)!=0){
		if(length(data.matrix_series)!=0){

		prnt<-paste("Series-Matrix",tsmt_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_data.matrix_series[[tsmt_i]]<<-data.matrix_series
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

		autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
		if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Log Transformation...",label="",min=0,max=1,initial=0,width=500)
			Med<-median(data.matrix_series,na.rm=TRUE)
			if(Med>16)
			{
				data.matrixLog<-log2(data.matrix_series)
			} else {
				data.matrixLog<-data.matrix_series
			}
setTkProgressBar(pb_ma,0.1,title="KNN Imputation...",label="")
			na.length<-length(which(is.na(data.matrixLog)==TRUE))
			if(na.length > 0)data.matrixImp<-impute.knn(data.matrixLog,k=10,rowmax=0.5,colmax=0.3)$data
			if(na.length <= 0)data.matrixImp<-data.matrixLog
			data.matrixNorm=NULL;
			data.matrixNorm.m=NULL;use.data.matrixNorm.m=NULL;xf=NULL;data.matrixNorm.f=NULL;data.matrixNorm.s=NULL;
			design=NULL;design_S=NULL;ttx=NULL;DE_S=NULL;DE_S_2=NULL;use_DE_S=NULL;
			sample.dist_S=NULL;
#			pca_S=NULL;
#			sample.clust_S=NULL;
#			clas_S=NULL;
			try(({folder_S[[tsmt_i]]<-folder_S[[tsmt_i]];data.matrixImp<-data.matrixImp;heatcol<-heatcol;}),silent=TRUE)
			setwd(folder_S[[tsmt_i]])	
			if(length(data.matrixImp)!=0){
setTkProgressBar(pb_ma,0.2,title="Normalization...",label="")
				data.matrixNorm<-normalizeBetweenArrays(data.matrixImp,method="quantile")
				tmp<-aggregate(data.matrixNorm,list(rownames(data.matrixNorm)),median)
				rownames(tmp)<-tmp[,1]
				data.matrixNorm.m<-as.matrix(tmp[,-1])
				try({
					if(length(rownames(data.matrixNorm.m))==0){
						rownames(data.matrixNorm.m)<-1:length(data.matrixNorm.m[,1])
					} else {
						rownames(data.matrixNorm.m)<-as.character(rownames(data.matrixNorm.m))
					}
				},silent=TRUE)
			}
			ta_smt_norm<-tclArray()
			for(i in 0:dim(data.matrixNorm.m)[1]){ for(j in 0:dim(data.matrixNorm.m)[2]){ 
				if(i==0){ ta_smt_norm[[i,j]]<-colnames(data.matrixNorm.m)[j] } else { 
				if(j==0){ ta_smt_norm[[i,j]]<-rownames(data.matrixNorm.m)[i] } else {
					tem<-data.matrixNorm.m[i,j]
					ta_smt_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
				} }
			} }

#			data.matrixNorm.m<<-data.matrixNorm.m
			data.matrixNorm.m2[[tsmt_i]]<<-data.matrixNorm.m

			ta_smt_norm2[[tsmt_i]]<<-ta_smt_norm
#				

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(data.matrixNorm.m2[[tsmt_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrixNorm.m2[[tsmt_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.data.matrixNorm.m<-data.matrixNorm.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.data.matrixNorm.m,fff)
					dat.fo<-use.data.matrixNorm.m[i,]
					i<-genefilter(-use.data.matrixNorm.m,fff)
					dat.fu<-use.data.matrixNorm.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_S[[tsmt_i]]<<-fit2
					data.matrixNorm.f<-fit2
				}),silent=TRUE)
				if(length(data.matrixNorm.f)==0)
				{
					try(({
						rsd<-rowSds(use.data.matrixNorm.m)
						i<-rsd>=2
						dat.f<-use.data.matrixNorm.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_S[[tsmt_i]]<<-fit2
						data.matrixNorm.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(data.matrixNorm.f)	
				rownames(ps1)<-rownames(data.matrixNorm.f)
				data.matrixNorm.f2=NULL
				data.matrixNorm.f2<-ps1
				ta_smt_filt=NULL
				ta_smt_filt<-tclArray()
				for(i in 0:dim(data.matrixNorm.f2)[1]){ for(j in 0:dim(data.matrixNorm.f2)[2]){ 
					if(i==0){ ta_smt_filt[[i,j]]<-colnames(data.matrixNorm.f2)[j] } else { 
					if(j==0){ ta_smt_filt[[i,j]]<-rownames(data.matrixNorm.f2)[i] } else {
					tem<-data.matrixNorm.f2[i,j]
					ta_smt_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				data.matrixNorm.f3[[tsmt_i]]<<-data.matrixNorm.f2
				
				ta_smt_filt2[[tsmt_i]]<<-ta_smt_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(data.matrixNorm.f)!=0){
					err<-try({
						ttx<-topTable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrixNorm.m)==TRUE)]
						data.matrixNorm.s<-use.data.matrixNorm.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrixNorm.m)==TRUE)])
						data.matrixNorm.s<-data.matrixNorm.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				data.matrixNorm.s2=NULL
				data.matrixNorm.s2<-data.matrixNorm.s
				ta_smt_stat=NULL
				ta_smt_stat<-tclArray()
				for(i in 0:dim(data.matrixNorm.s2)[1]){ for(j in 0:dim(data.matrixNorm.s2)[2]){ 
					if(i==0){ ta_smt_stat[[i,j]]<-colnames(data.matrixNorm.s2)[j] } else { 
					if(j==0){ ta_smt_stat[[i,j]]<-rownames(data.matrixNorm.s2)[i] } else {
					tem<-data.matrixNorm.s2[i,j]
					ta_smt_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				data.matrixNorm.s3[[tsmt_i]]<<-data.matrixNorm.s2
				
				ta_smt_stat2[[tsmt_i]]<<-ta_smt_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_S<-topTable(data.matrixNorm.f),silent=TRUE)	
				err<-try(DE_SAll[[tsmt_i]]<<-topTable(data.matrixNorm.f,number=nrow(data.matrixNorm.f)),silent=TRUE)
				ta_smt_dge=NULL
				ta_smt_dge<-tclArray()
				for(i in 0:dim(DE_S)[1]){ for(j in 0:dim(DE_S)[2]){ 
					if(i==0){ ta_smt_dge[[i,j]]<-colnames(DE_S)[j] } else { 
					if(j==0){ ta_smt_dge[[i,j]]<-rownames(DE_S)[i] } else {
					tem<-DE_S[i,j]
					ta_smt_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_S_2[[tsmt_i]]<<-DE_S
				
				ta_smt_dge2[[tsmt_i]]<<-ta_smt_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_S[[tsmt_i]]<<-prcomp(t(data.matrixNorm.m[[tsmt_i]]))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_S[[tsmt_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_S[[tsmt_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(data.matrixNorm.m)) 
				sample.dist_S<-as.dist(1-cor(data.matrixNorm.m,method="pearson"))
				sample.clust_S[[tsmt_i]]<<-hclust(sample.dist_S,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_S[[tsmt_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_S[[tsmt_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_S<-as.matrix(DE_S),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_S)<-as.character(DE_S)
					use_DE_S<-as.matrix(DE_S)
				}
				rnames<-rownames(use_DE_S)
				clas_S[[tsmt_i]]<<-as.matrix(data.matrixNorm.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_S[[tsmt_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_S[[tsmt_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("genefilter",quietly=TRUE)){BiocManager::install("genefilter",ask=FALSE,update=FALSE);library(genefilter)} else {library(genefilter)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("impute",quietly=TRUE)){BiocManager::install("impute",ask=FALSE,update=FALSE);library(impute)} else {library(impute)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("limma",quietly=TRUE)){BiocManager::install("limma",ask=FALSE,update=FALSE);library(limma)} else {library(limma)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("Biobase",quietly=TRUE)){BiocManager::install("Biobase",ask=FALSE,update=FALSE);library(Biobase)} else {library(Biobase)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("stats",quietly=TRUE)){install.packages("stats",ask=FALSE,update=FALSE);library(stats)} else {library(stats)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk",ask=FALSE,update=FALSE);library(tcltk)} else {library(tcltk)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
dmsm()
