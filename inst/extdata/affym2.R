affym<-function(h,...){
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
	taf_i<<-taf_i+1
	folder_Affy[[taf_i]]<<-folder
	if(folder_Affy[[taf_i]]==""){
		folder_Affy<-NULL
	} else { 
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			setwd(folder_Affy[[taf_i]])
			datAffy<-ReadAffy()
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			new_datAffy=NULL
			new_datAffy<-datAffy
			if(length(new_datAffy)!=0){datAffy=NULL;datAffy<-new_datAffy}
			
			ann_Affy[[taf_i]]<<-annotation(datAffy)
		 }),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_datAffy)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	}
	if(length(datAffy)!=0){
		if(length(new_datAffy)!=0){

		prnt<-paste("Affymetrix",taf_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_datAffy[[taf_i]]<<-datAffy
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))
#p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(paste("Affymetrix",taf_i,sep="_")))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				dat2Affy.m=NULL;use.dat2Affy.m=NULL;xf=NULL;dat2Affy.f=NULL;dat2Affy.s=NULL;
				design=NULL;design_Affy=NULL;ttx=NULL;DE_Affy=NULL;DE_Affy_2=NULL;use_DE_Affy=NULL;				
				sample.dist_Affy=NULL;
				try(({folder_Affy[[taf_i]]<-folder_Affy[[taf_i]];datAffy<-datAffy;heatcol<-heatcol;}),silent=TRUE)
				setwd(folder_Affy[[taf_i]])
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				dat2Affy<-justRMA()
				dat2Affy.m<-exprs(dat2Affy)
				dat2Affy.m<-as.matrix(dat2Affy.m)
				try({
					if(length(rownames(dat2Affy.m))==0){
						rownames(dat2Affy.m)<-1:length(dat2Affy.m[,1])
					} else {
						rownames(dat2Affy.m)<-as.character(rownames(dat2Affy.m))
					}
				},silent=TRUE)
				ta_affy_norm=NULL
				ta_affy_norm<-tclArray()
				for(i in 0:dim(dat2Affy.m)[1]){ for(j in 0:dim(dat2Affy.m)[2]){ 
					if(i==0){ ta_affy_norm[[i,j]]<-colnames(dat2Affy.m)[j] } else { 
					if(j==0){ ta_affy_norm[[i,j]]<-rownames(dat2Affy.m)[i] } else {
					tem<-dat2Affy.m[i,j]
					ta_affy_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
#				dat2Affy.m<<-dat2Affy.m
				dat2Affy.m2[[taf_i]]<<-dat2Affy.m
				
				ta_affy_norm2[[taf_i]]<<-ta_affy_norm
#				
				
p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(dat2Affy.m2[[taf_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(dat2Affy.m2[[taf_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.dat2Affy.m<-dat2Affy.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.dat2Affy.m,fff)
					dat.fo<-use.dat2Affy.m[i,]
					i<-genefilter(-use.dat2Affy.m,fff)
					dat.fu<-use.dat2Affy.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_Affy[[taf_i]]<<-fit2
					dat2Affy.f<-fit2
				}),silent=TRUE)
				if(length(dat2Affy.f)==0)
				{
					try(({
						rsd<-rowSds(use.dat2Affy.m)
						i<-rsd>=2
						dat.f<-use.dat2Affy.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_Affy[[taf_i]]<<-fit2
						dat2Affy.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(dat2Affy.f)
				rownames(ps1)<-rownames(dat2Affy.f)
				dat2Affy.f2=NULL
				dat2Affy.f2<-ps1
				ta_affy_filt=NULL
				ta_affy_filt<-tclArray()
				for(i in 0:dim(dat2Affy.f2)[1]){ for(j in 0:dim(dat2Affy.f2)[2]){ 
					if(i==0){ ta_affy_filt[[i,j]]<-colnames(dat2Affy.f2)[j] } else { 
					if(j==0){ ta_affy_filt[[i,j]]<-rownames(dat2Affy.f2)[i] } else {
					tem<-dat2Affy.f2[i,j]
					ta_affy_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }

				
				dat2Affy.f3[[taf_i]]<<-dat2Affy.f2
				
				ta_affy_filt2[[taf_i]]<<-ta_affy_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(dat2Affy.f)!=0){
					err<-try({
						ttx<-topTable(dat2Affy.f,number=nrow(use.dat2Affy.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(dat2Affy.m)==TRUE)]
						dat2Affy.s<-use.dat2Affy.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(dat2Affy.f,number=nrow(use.dat2Affy.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(dat2Affy.m)==TRUE)])
						dat2Affy.s<-dat2Affy.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				dat2Affy.s2=NULL
				dat2Affy.s2<-dat2Affy.s
				ta_affy_stat=NULL
				ta_affy_stat<-tclArray()
				for(i in 0:dim(dat2Affy.s2)[1]){ for(j in 0:dim(dat2Affy.s2)[2]){ 
					if(i==0){ ta_affy_stat[[i,j]]<-colnames(dat2Affy.s2)[j] } else { 
					if(j==0){ ta_affy_stat[[i,j]]<-rownames(dat2Affy.s2)[i] } else {
					tem<-dat2Affy.s2[i,j]
					ta_affy_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }

				
				dat2Affy.s3[[taf_i]]<<-dat2Affy.s2
				
				ta_affy_stat2[[taf_i]]<<-ta_affy_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_Affy<-topTable(dat2Affy.f),silent=TRUE)	
				err<-try(DE_AffyAll[[taf_i]]<<-topTable(dat2Affy.f,number=nrow(dat2Affy.f)),silent=TRUE)	
				ta_affy_dge=NULL
				ta_affy_dge<-tclArray()
				for(i in 0:dim(DE_Affy)[1]){ for(j in 0:dim(DE_Affy)[2]){ 
					if(i==0){ ta_affy_dge[[i,j]]<-colnames(DE_Affy)[j] } else { 
					if(j==0){ ta_affy_dge[[i,j]]<-rownames(DE_Affy)[i] } else {
					tem<-DE_Affy[i,j]
					ta_affy_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_Affy_2[[taf_i]]<<-DE_Affy
				
				ta_affy_dge2[[taf_i]]<<-ta_affy_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_Affy[[taf_i]]<<-prcomp(t(dat2Affy.m))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_Affy[[taf_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Affy[[taf_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(dat2Affy.m))
				sample.dist_Affy<-as.dist(1-cor(dat2Affy.m,method="pearson"))
				sample.clust_Affy[[taf_i]]<<-hclust(sample.dist_Affy,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_Affy[[taf_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Affy[[taf_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_Affy<-as.matrix(DE_Affy),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_Affy)<-as.character(rownames(DE_Affy))
					use_DE_Affy<-as.matrix(DE_Affy)
				}
				rnames<-rownames(use_DE_Affy)
				clas_Affy[[taf_i]]<<-as.matrix(use.dat2Affy.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_Affy[[taf_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_Affy[[taf_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("affy",quietly=TRUE)){BiocManager::install("affy",ask=FALSE,update=FALSE);library(affy)} else {library(affy)}
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
affym()
