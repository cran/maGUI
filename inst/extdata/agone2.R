agone<-function(h,...){
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
	tag1_i<<-tag1_i+1
	folder_Ag1[[tag1_i]]<<-folder
	if(folder_Ag1[[tag1_i]]==""){
		folder_Ag1<-NULL
	} else {
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			new_datAgOne=NULL
			setwd(folder_Ag1[[tag1_i]])
			new_datAgOne<-read.maimages(files=dir(folder_Ag1[[tag1_i]]),columns=list(G="gMeanSignal",Gb="gBGMedianSignal",R="gMeanSignal",Rb="gBGMedianSignal"))
			if(length(new_datAgOne)!=0){datAgOne=NULL;datAgOne<-new_datAgOne}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
			z<-as(datAgOne,"NChannelSet")
			ann_Ag1[[tag1_i]]<<-annotation(z)
		 }),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_datAgOne)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	}
	if(length(datAgOne)!=0){
		if(length(new_datAgOne)!=0){

		prnt<-paste("Agilent-OneColor",tag1_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_datAgOne[[tag1_i]]<<-datAgOne
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				datAgOne2.m=NULL;use.datAgOne2.m=NULL;xf=NULL;datAgOne2.f=NULL;datAgOne2.s=NULL;
				design=NULL;design_Ag1=NULL;ttx=NULL;DE_Ag1=NULL;DE_Ag1_2=NULL;use_DE_Ag1=NULL;
				sample.dist_Ag1=NULL;
#				pca_Ag1=NULL;
#				sample.clust_Ag1=NULL;
#				clas_Ag1=NULL;
				try(({folder_Ag1[[tag1_i]]<-folder_Ag1[[tag1_i]];datAgOne<-datAgOne;heatcol<-heatcol; }),silent=TRUE)
				setwd(folder_Ag1[[tag1_i]])	
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				try(datAgOne2<-backgroundCorrect(datAgOne,"normexp"),silent=TRUE)
				if(length(datAgOne2)==0)datAgOne2<-datAgOne
				datAgOne2<-normalizeBetweenArrays(datAgOne2$R,method="quantile")
				datAgOne2<-log(datAgOne2)
				datAgOne2.m<-datAgOne2
				datAgOne2.m<-as.matrix(datAgOne2.m)
				try({
					if(length(rownames(datAgOne2.m))==0){
						rownames(datAgOne2.m)<-1:length(datAgOne2.m[,1])
					} else {
						rownames(datAgOne2.m)<-as.character(rownames(datAgOne2.m))
					}
				},silent=TRUE)
				ta_ag1_norm=NULL
				ta_ag1_norm<-tclArray()
				for(i in 0:dim(datAgOne2.m)[1]){ for(j in 0:dim(datAgOne2.m)[2]){ 
					if(i==0){ ta_ag1_norm[[i,j]]<-colnames(datAgOne2.m)[j] } else { 
					if(j==0){ ta_ag1_norm[[i,j]]<-rownames(datAgOne2.m)[i] } else {
					tem<-datAgOne2.m[i,j]
					ta_ag1_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
#				datAgOne2.m[[tag1_i]]<<-datAgOne2.m[[tag1_i]]
				datAgOne2.m2[[tag1_i]]<<-datAgOne2.m
				
				ta_ag1_norm2[[tag1_i]]<<-ta_ag1_norm
#				

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datAgOne2.m2[[tag1_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgOne2.m2[[tag1_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.datAgOne2.m<-datAgOne2.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.datAgOne2.m,fff)
					dat.fo<-use.datAgOne2.m[i,]
					i<-genefilter(-use.datAgOne2.m,fff)
					dat.fu<-use.datAgOne2.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_Ag1[[tag1_i]]<<-fit2
					datAgOne2.f<-fit2
				}),silent=TRUE)
				if(length(datAgOne2.f)==0)
				{
					try(({
						rsd<-rowSds(use.datAgOne2.m)
						i<-rsd>=2
						dat.f<-use.datAgOne2.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_Ag1[[tag1_i]]<<-fit2
						datAgOne2.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(datAgOne2.f)	
				rownames(ps1)<-rownames(datAgOne2.f)
				datAgOne2.f2=NULL
				datAgOne2.f2<-ps1
				ta_ag1_filt=NULL
				ta_ag1_filt<-tclArray()
				for(i in 0:dim(datAgOne2.f2)[1]){ for(j in 0:dim(datAgOne2.f2)[2]){ 
					if(i==0){ ta_ag1_filt[[i,j]]<-colnames(datAgOne2.f2)[j] } else { 
					if(j==0){ ta_ag1_filt[[i,j]]<-rownames(datAgOne2.f2)[i] } else {
					tem<-datAgOne2.f2[i,j]
					ta_ag1_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
			
				datAgOne2.f3[[tag1_i]]<<-datAgOne2.f2
				
				ta_ag1_filt2[[tag1_i]]<<-ta_ag1_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(datAgOne2.f)!=0){
					err<-try({
						ttx<-topTable(datAgOne2.f,number=nrow(use.datAgOne2.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datAgOne2.m)==TRUE)]
						datAgOne2.s<-use.datAgOne2.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datAgOne2.f,number=nrow(use.datAgOne2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datAgOne2.m)==TRUE)])
						datAgOne2.s<-datAgOne2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datAgOne2.s2=NULL
				datAgOne2.s2<-datAgOne2.s
				ta_ag1_stat=NULL
				ta_ag1_stat<-tclArray()
				for(i in 0:dim(datAgOne2.s2)[1]){ for(j in 0:dim(datAgOne2.s2)[2]){ 
					if(i==0){ ta_ag1_stat[[i,j]]<-colnames(datAgOne2.s2)[j] } else { 
					if(j==0){ ta_ag1_stat[[i,j]]<-rownames(datAgOne2.s2)[i] } else {
					tem<-datAgOne2.s2[i,j]
					ta_ag1_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
			
				datAgOne2.s3[[tag1_i]]<<-datAgOne2.s2
				
				ta_ag1_stat2[[tag1_i]]<<-ta_ag1_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_Ag1<-topTable(datAgOne2.f),silent=TRUE)	
				err<-try(DE_Ag1All[[tag1_i]]<<-topTable(datAgOne2.f,number=nrow(datAgOne2.f)),silent=TRUE)	
				ta_ag1_dge=NULL
				ta_ag1_dge<-tclArray()
				for(i in 0:dim(DE_Ag1)[1]){ for(j in 0:dim(DE_Ag1)[2]){ 
					if(i==0){ ta_ag1_dge[[i,j]]<-colnames(DE_Ag1)[j] } else { 
					if(j==0){ ta_ag1_dge[[i,j]]<-rownames(DE_Ag1)[i] } else {
					tem<-DE_Ag1[i,j]
					ta_ag1_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_Ag1_2[[tag1_i]]<<-DE_Ag1
				
				ta_ag1_dge2[[tag1_i]]<<-ta_ag1_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_Ag1[[tag1_i]]<<-prcomp(t(datAgOne2.m[[tag1_i]]))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_Ag1[[tag1_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Ag1[[tag1_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(datAgOne2.m)) 
				sample.dist_Ag1<-as.dist(1-cor(datAgOne2.m,method="pearson"))
				sample.clust_Ag1[[tag1_i]]<<-hclust(sample.dist_Ag1,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_Ag1[[tag1_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Ag1[[tag1_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_Ag1<-as.matrix(DE_Ag1),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_Ag1)<-as.character(rownames(DE_Ag1))
					use_DE_Ag1<-as.matrix(DE_Ag1)
				}
				rnames<-rownames(use_DE_Ag1)
				clas_Ag1[[tag1_i]]<<-as.matrix(use.datAgOne2.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_Ag1[[tag1_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_Ag1[[tag1_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
agone()
