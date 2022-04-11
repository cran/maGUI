agtwo<-function(h,...){
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
	tag2_i<<-tag2_i+1
	folder_Ag2[[tag2_i]]<<-folder
	if(folder_Ag2[[tag2_i]]==""){
		folder_Ag2<-NULL
	} else {
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			new_datAgTwo=NULL
			setwd(folder_Ag2[[tag2_i]])
			new_datAgTwo<-read.maimages(files=dir(folder_Ag2[[tag2_i]]),source="agilent")
			if(length(new_datAgTwo)!=0){datAgTwo=NULL;datAgTwo<-new_datAgTwo}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
			z<-as(datAgTwo,"NChannelSet")
			ann_Ag2[[tag2_i]]<<-annotation(z)
		 }),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_datAgTwo)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	}
	if(length(datAgTwo)!=0){
		if(length(new_datAgTwo)!=0){

		prnt<-paste("Agilent-TwoColor",tag2_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_datAgTwo[[tag2_i]]<<-datAgTwo
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				datAgTwo2.m=NULL;use.datAgTwo2.m=NULL;xf=NULL;datAgTwo2.f=NULL;datAgTwo2.s=NULL;
				design=NULL;design_Ag2=NULL;ttx=NULL;DE_Ag2=NULL;DE_Ag2_2=NULL;use_DE_Ag2=NULL;
				sample.dist_Ag2=NULL;
#				pca_Ag2=NULL;
#				sample.clust_Ag2=NULL;
#				clas_Ag2=NULL;
				try(({folder_Ag2[[tag2_i]]<-folder_Ag2[[tag2_i]];datAgTwo<-datAgTwo;heatcol<-heatcol; }),silent=TRUE)
				setwd(folder_Ag2[[tag2_i]])
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				try(datAgTwo2<-backgroundCorrect(datAgTwo,"normexp"),silent=TRUE)
				if(length(datAgTwo2)==0)datAgTwo2<-datAgTwo
				datAgTwo2<-normalizeBetweenArrays(datAgTwo2$R,method="quantile")
				datAgTwo2<-log(datAgTwo2)
				datAgTwo2.m<-datAgTwo2
				datAgTwo2.m<-as.matrix(datAgTwo2.m)
				try({
					if(length(rownames(datAgTwo2.m))==0){
						rownames(datAgTwo2.m)<-1:length(datAgTwo2.m[,1])
					} else {
						rownames(datAgTwo2.m)<-as.character(rownames(datAgTwo2.m))
					}
				},silent=TRUE)
				ta_ag2_norm=NULL
				ta_ag2_norm<-tclArray()
				for(i in 0:dim(datAgTwo2.m)[1]){ for(j in 0:dim(datAgTwo2.m)[2]){ 
					if(i==0){ ta_ag2_norm[[i,j]]<-colnames(datAgTwo2.m)[j] } else { 
					if(j==0){ ta_ag2_norm[[i,j]]<-rownames(datAgTwo2.m)[i] } else {
					tem<-datAgTwo2.m[i,j]
					ta_ag2_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }

#				datAgTwo2.m[[tag2_i]]<<-datAgTwo2.m[[tag2_i]]
				datAgTwo2.m2[[tag2_i]]<<-datAgTwo2.m
				
				ta_ag2_norm2[[tag2_i]]<<-ta_ag2_norm
#				

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datAgTwo2.m2[[tag2_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgTwo2.m2[[tag2_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.datAgTwo2.m<-datAgTwo2.m
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.datAgTwo2.m,fff)
					dat.fo<-use.datAgTwo2.m[i,]
					i<-genefilter(-use.datAgTwo2.m,fff)
					dat.fu<-use.datAgTwo2.m[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_Ag2[[tag2_i]]<<-fit2
					datAgTwo2.f<-fit2
				}),silent=TRUE)
				if(length(datAgTwo2.f)==0)
				{
					try(({
						rsd<-rowSds(use.datAgTwo2.m)
						i<-rsd>=2
						dat.f<-use.datAgTwo2.m[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_Ag2[[tag2_i]]<<-fit2
						datAgTwo2.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(datAgTwo2.f)
				rownames(ps1)<-rownames(datAgTwo2.f)
				datAgTwo2.f2=NULL
				datAgTwo2.f2<-ps1
				ta_ag2_filt=NULL
				ta_ag2_filt<-tclArray()
				for(i in 0:dim(datAgTwo2.f2)[1]){ for(j in 0:dim(datAgTwo2.f2)[2]){ 
					if(i==0){ ta_ag2_filt[[i,j]]<-colnames(datAgTwo2.f2)[j] } else { 
					if(j==0){ ta_ag2_filt[[i,j]]<-rownames(datAgTwo2.f2)[i] } else {
					tem<-datAgTwo2.f2[i,j]
					ta_ag2_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				datAgTwo2.f3[[tag2_i]]<<-datAgTwo2.f2
				
				ta_ag2_filt2[[tag2_i]]<<-ta_ag2_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(datAgTwo2.f)!=0){
					err<-try({
							ttx<-topTable(datAgTwo2.f,number=nrow(use.datAgTwo2.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datAgTwo2.m)==TRUE)]
							datAgTwo2.s<-use.datAgTwo2.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datAgTwo2.f,number=nrow(use.datAgTwo2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datAgTwo2.m)==TRUE)])
						datAgTwo2.s<-datAgTwo2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datAgTwo2.s2=NULL
				datAgTwo2.s2<-datAgTwo2.s
				ta_ag2_stat=NULL
				ta_ag2_stat<-tclArray()
				for(i in 0:dim(datAgTwo2.s2)[1]){ for(j in 0:dim(datAgTwo2.s2)[2]){ 
					if(i==0){ ta_ag2_stat[[i,j]]<-colnames(datAgTwo2.s2)[j] } else { 
					if(j==0){ ta_ag2_stat[[i,j]]<-rownames(datAgTwo2.s2)[i] } else {
					tem<-datAgTwo2.s2[i,j]
					ta_ag2_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				datAgTwo2.s3[[tag2_i]]<<-datAgTwo2.s2
				
				ta_ag2_stat2[[tag2_i]]<<-ta_ag2_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_Ag2<-topTable(datAgTwo2.f),silent=TRUE)	
				err<-try(DE_Ag2All[[tag2_i]]<<-topTable(datAgTwo2.f,number=nrow(datAgTwo2.f)),silent=TRUE)
				ta_ag2_dge=NULL
				ta_ag2_dge<-tclArray()
				for(i in 0:dim(DE_Ag2)[1]){ for(j in 0:dim(DE_Ag2)[2]){ 
					if(i==0){ ta_ag2_dge[[i,j]]<-colnames(DE_Ag2)[j] } else { 
					if(j==0){ ta_ag2_dge[[i,j]]<-rownames(DE_Ag2)[i] } else {
					tem<-DE_Ag2[i,j]
					ta_ag2_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_Ag2_2[[tag2_i]]<<-DE_Ag2
				
				ta_ag2_dge2[[tag2_i]]<<-ta_ag2_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_Ag2[[tag2_i]]<<-prcomp(t(datAgTwo2.m))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_Ag2[[tag2_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Ag2[[tag2_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(datAgTwo2.m)) 
				sample.dist_Ag2<-as.dist(1-cor(datAgTwo2.m,method="pearson"))
				sample.clust_Ag2[[tag2_i]]<<-hclust(sample.dist_Ag2,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_Ag2[[tag2_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Ag2[[tag2_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_Ag2<-as.matrix(DE_Ag2),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_Ag2)<-as.character(DE_Ag2)
					use_DE_Ag2<-as.matrix(DE_Ag2)
				}
				rnames<-rownames(use_DE_Ag2)
				clas_Ag2[[tag2_i]]<<-as.matrix(use.datAgTwo2.m[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_Ag2[[tag2_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_Ag2[[tag2_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
agtwo()
