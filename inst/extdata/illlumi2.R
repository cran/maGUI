illlumi<-function(h,...){
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
	till_i<<-till_i+1
	folder_Il_L[[till_i]]<<-dirname(file)
	if(file==""){
		folder_Il_L<-NULL
	} else { 
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
			new_lumi_data=NULL
			setwd(folder_Il_L[[till_i]])
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			new_lumi_data<-lumiR(file)
			summary(new_lumi_data,'QC')
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			if(length(new_lumi_data)!=0){lumi_data=NULL;lumi_data<-new_lumi_data}
			ann_Il_L[[till_i]]<<-annotation(lumi_data)
		 }),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_lumi_data)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	}
	if(length(lumi_data)!=0){
			if(length(new_lumi_data)!=0){

		prnt<-paste("Illumina-Lumi",till_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_lumi_data[[till_i]]<<-lumi_data
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				lumi_NQ.m=NULL;use.lumi_NQ.m=NULL;xf=NULL;lumi_NQ.f=NULL;lumi_NQ.s=NULL;
				design=NULL;design_Il_L=NULL;ttx=NULL;DE_Il_L=NULL;DE_Il_L_2=NULL;use_DE_Il_L=NULL;
				sample.dist_Il_L=NULL;
#				pca_Il_L=NULL;
#				sample.clust_Il_L=NULL;
#				clas_Il_L=NULL;
				try(({folder_Il_L[[till_i]]<-folder_Il_L[[till_i]];lumi_data<-lumi_data;heatcol<-heatcol;}),silent=TRUE)
				setwd(folder_Il_L[[till_i]])	
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				lumi_NQ<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
				summary(lumi_NQ,'QC')
				lumi_NQ.m<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
				lumi_NQ.m<-as.matrix(lumi_NQ.m)
				try({
					if(length(rownames(lumi_NQ.m))==0){
						rownames(lumi_NQ.m)<-1:length(lumi_NQ.m[,1])
					} else {
						rownames(lumi_NQ.m)<-as.character(rownames(lumi_NQ.m))
					}
				},silent=TRUE)

				ta_ill_norm=NULL
				ta_ill_norm<-tclArray()
				for(i in 0:dim(lumi_NQ.m)[1]){ for(j in 0:dim(lumi_NQ.m)[2]){ 
					if(i==0){ ta_ill_norm[[i,j]]<-colnames(lumi_NQ.m)[j] } else { 
					if(j==0){ ta_ill_norm[[i,j]]<-rownames(lumi_NQ.m)[i] } else {
					tem<-lumi_NQ.m[i,j]
					ta_ill_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				lumi_NQ.m2[[till_i]]<<-lumi_NQ.m
				
				ta_ill_norm2[[till_i]]<<-ta_ill_norm

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(lumi_NQ.m2[[till_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(lumi_NQ.m2[[till_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.lumi_NQ.m<-lumi_NQ.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.lumi_NQ.m,fff)
					dat.fo<-use.lumi_NQ.m[i,]
					i<-genefilter(-use.lumi_NQ.m,fff)
					dat.fu<-use.lumi_NQ.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_Il_L[[till_i]]<<-fit2
					lumi_NQ.f<-fit2
				}),silent=TRUE)
				if(length(lumi_NQ.f)==0)
				{
					try(({
						rsd<-rowSds(use.lumi_NQ.m)
						i<-rsd>=2
						dat.f<-use.lumi_NQ.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_Il_L[[till_i]]<<-fit2
						lumi_NQ.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(lumi_NQ.f)	
				rownames(ps1)<-rownames(lumi_NQ.f)
				lumi_NQ.f2=NULL
				lumi_NQ.f2<-ps1
				ta_ill_filt=NULL
				ta_ill_filt<-tclArray()
				for(i in 0:dim(lumi_NQ.f2)[1]){ for(j in 0:dim(lumi_NQ.f2)[2]){ 
					if(i==0){ ta_ill_filt[[i,j]]<-colnames(lumi_NQ.f2)[j] } else { 
					if(j==0){ ta_ill_filt[[i,j]]<-rownames(lumi_NQ.f2)[i] } else {
					tem<-lumi_NQ.f2[i,j]
					ta_ill_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				lumi_NQ.f3[[till_i]]<<-lumi_NQ.f2
				
				ta_ill_filt2[[till_i]]<<-ta_ill_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(lumi_NQ.f)!=0){
					err<-try({
							ttx<-topTable(lumi_NQ.f,number=nrow(use.lumi_NQ.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(lumi_NQ.m)==TRUE)]
							lumi_NQ.s<-use.lumi_NQ.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(lumi_NQ.f,number=nrow(use.lumi_NQ.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(lumi_NQ.m)==TRUE)])
						lumi_NQ.s<-lumi_NQ.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				lumi_NQ.s2=NULL
				lumi_NQ.s2<-lumi_NQ.s
				ta_ill_stat=NULL
				ta_ill_stat<-tclArray()
				for(i in 0:dim(lumi_NQ.s2)[1]){ for(j in 0:dim(lumi_NQ.s2)[2]){ 
					if(i==0){ ta_ill_stat[[i,j]]<-colnames(lumi_NQ.s2)[j] } else { 
					if(j==0){ ta_ill_stat[[i,j]]<-rownames(lumi_NQ.s2)[i] } else {
					tem<-lumi_NQ.s2[i,j]
					ta_ill_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				lumi_NQ.s3[[till_i]]<<-lumi_NQ.s2
				
				ta_ill_stat2[[till_i]]<<-ta_ill_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_Il_L<-topTable(lumi_NQ.f),silent=TRUE)	
				err<-try(DE_Il_LAll[[till_i]]<<-topTable(lumi_NQ.f,number=nrow(lumi_NQ.f)),silent=TRUE)
				ta_ill_dge=NULL
				ta_ill_dge<-tclArray()
				for(i in 0:dim(DE_Il_L)[1]){ for(j in 0:dim(DE_Il_L)[2]){ 
					if(i==0){ ta_ill_dge[[i,j]]<-colnames(DE_Il_L)[j] } else { 
					if(j==0){ ta_ill_dge[[i,j]]<-rownames(DE_Il_L)[i] } else {
					tem<-DE_Il_L[i,j]
					ta_ill_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_Il_L_2[[till_i]]<<-DE_Il_L
				
				ta_ill_dge2[[till_i]]<<-ta_ill_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_Il_L[[till_i]]<<-prcomp(t(lumi_NQ.m))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_Il_L[[till_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Il_L[[till_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(lumi_NQ.m)) 
				sample.dist_Il_L<-as.dist(1-cor(lumi_NQ.m,method="pearson"))
				sample.clust_Il_L[[till_i]]<<-hclust(sample.dist_Il_L,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_Il_L[[till_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Il_L[[till_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_Il_L<-as.matrix(DE_Il_L),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_Il_L)<-as.character(DE_Il_L)
					use_DE_Il_L<-as.matrix(DE_Il_L)
				}
				rnames<-rownames(use_DE_Il_L)
				clas_Il_L[[till_i]]<<-as.matrix(lumi_NQ.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_Il_L[[till_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_Il_L[[till_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
if(!requireNamespace("lumi",quietly=TRUE)){BiocManager::install("lumi",ask=FALSE,update=FALSE);library(lumi)} else {library(lumi)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("genefilter",quietly=TRUE)){BiocManager::install("genefilter",ask=FALSE,update=FALSE);library(genefilter)} else {library(genefilter)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("limma",quietly=TRUE)){BiocManager::install("limma",ask=FALSE,update=FALSE);library(limma)} else {library(limma)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("stats",quietly=TRUE)){install.packages("stats",ask=FALSE,update=FALSE);library(stats)} else {library(stats)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk",ask=FALSE,update=FALSE);library(tcltk)} else {library(tcltk)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
illlumi()
