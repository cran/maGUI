dmo<-function(h,...){
	tonl_i<<-tonl_i+1
	w_imp1<-tktoplevel()
	tkwm.title(w_imp1,"Online Data")
	tkwm.resizable(w_imp1,0,0)
	frame1_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
	tkpack(ttklabel(frame1_imp1,text='Enter a GSE number:'),side="left",padx=12)
	exp_file<-tclVar("")
	c1<-ttkentry(frame1_imp1,textvariable=exp_file,width=25)
	tkpack(c1,side="left",anchor="e")
	tkpack(frame1_imp1,expand=TRUE,fill="both",side="top")
	frame2_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
	r_but1<-ttkbutton(frame2_imp1,text="Ok")
	tkpack(r_but1,side="top",anchor="e",pady=5)
	tkpack(frame2_imp1,expand=TRUE,fill="both",side="top")
	tkconfigure(r_but1,command=function(){
		tkwm.withdraw(w_imp1)
		if(length(grep("GSE",tclvalue(exp_file))!=0)){
			gse_no<-tclvalue(exp_file)
			folder<-tclvalue(tkchooseDirectory(title="Destination Folder"))
			folder_O[[tonl_i]]<<-folder
			setwd(folder)
			new_data.matrix_online=NULL;data.matrix_online=NULL
			try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				new_gse=NULL
				new_gse<-getGEO(gse_no,GSEMatrix=TRUE,destdir=folder_O[[tonl_i]])
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				gse<-new_gse[[1]]
				ann_O[[tonl_i]]<<-annotation(gse)
				ORF<-featureData(gse)@data$ORF
				use_probe<-which(is.na(ORF)==F & match(ORF,"",nomatch=0)==0)
				new_data.matrix_online<-exprs(gse)[use_probe,]
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
				data.matrix_online<-new_data.matrix_online
				rownames(new_data.matrix_online)<-ORF[use_probe]
				if(length(new_data.matrix_online)!=0){data.matrix_online=NULL;data.matrix_online<-new_data.matrix_online}
		 	}),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
			if(length(new_data.matrix_online)==0)tkmessageBox(title="Error",message="Could not load this type of data!",icon="error",type="ok")
		}
		if(length(data.matrix_online)!=0){
			if(length(new_data.matrix_online)!=0){
	
			prnt<-paste("Online-Data",tonl_i,sep="_")
			tree_prnts<<-c(tree_prnts,prnt)
			l_data.matrix_online[[tonl_i]]<<-data.matrix_online
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))
	
			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				Med<-median(data.matrix_online,na.rm=TRUE)
				if(Med>16)
				{
					data.matrix_onlineLog<-log2(data.matrix_online)
				} else {
					data.matrix_onlineLog<-data.matrix_online
				}
				l_data.matrix_onlineLog[[tonl_i]]<<-data.matrix_onlineLog
				na.length<-length(which(is.na(data.matrix_onlineLog)==TRUE))
				if(na.length > 0)l_data.matrix_onlineImp[[tonl_i]]<<-impute.knn(data.matrix_onlineLog,k=10,rowmax=0.5,colmax=0.3)$data
				if(na.length <= 0)l_data.matrix_onlineImp[[tonl_i]]<<-data.matrix_onlineLog
	
				data.matrix_onlineNorm=NULL;
				data.matrix_onlineNorm.m=NULL;use.data.matrix_onlineNorm.m=NULL;xf=NULL;data.matrix_onlineNorm.f=NULL;data.matrix_onlineNorm.s=NULL;
				design=NULL;design_O=NULL;ttx=NULL;DE_O=NULL;DE_O_2=NULL;use_DE_O=NULL;
				sample.dist_O=NULL;
#			pca_O=NULL;
#			sample.clust_O=NULL;	
#			clas_O=NULL;
				try(({data.matrix_onlineImp<-l_data.matrix_onlineImp[[tonl_i]];heatcol<-heatcol;}),silent=TRUE)
				if(length(data.matrix_onlineImp)!=0){
					data.matrix_onlineNorm<-normalizeBetweenArrays(data.matrix_onlineImp,method="quantile")
					tmp<-aggregate(data.matrix_onlineNorm,list(rownames(data.matrix_onlineNorm)),median)
					rownames(tmp)<-tmp[,1]
					data.matrix_onlineNorm.m<-as.matrix(tmp[,-1])
					try({
						if(length(rownames(data.matrix_onlineNorm.m))==0){
							rownames(data.matrix_onlineNorm.m)<-1:length(data.matrix_onlineNorm.m[,1])
						} else {
							rownames(data.matrix_onlineNorm.m)<-as.character(rownames(data.matrix_onlineNorm.m))
						}
					},silent=TRUE)
				}
				ta_onl_norm<-tclArray()
				for(i in 0:dim(data.matrix_onlineNorm.m)[1]){ for(j in 0:dim(data.matrix_onlineNorm.m)[2]){ 
					if(i==0){ ta_onl_norm[[i,j]]<-colnames(data.matrix_onlineNorm.m)[j] } else { 
					if(j==0){ ta_onl_norm[[i,j]]<-rownames(data.matrix_onlineNorm.m)[i] } else {
					tem<-data.matrix_onlineNorm.m[i,j]
					ta_onl_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
#			data.matrix_onlineNorm.m<<-data.matrix_onlineNorm.m
				data.matrix_onlineNorm.m2[[tonl_i]]<<-data.matrix_onlineNorm.m
				
				ta_onl_norm2[[tonl_i]]<<-ta_onl_norm
#				
	
p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
					if(is.null(ma_img)==TRUE){
						try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
						try(tkwm.withdraw(w_legend),silent=TRUE)
						ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...)	{boxplot(data.matrix_onlineNorm.m)},hscale=1.7,vscale=1.07)
						ma_img2<<-ma_img
						tkpack(ma_img2)
					} else {
						ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_onlineNorm.m)},hscale=1.7,vscale=1.07)
					}
p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
					use.data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m
					try(({
						fff<-pOverA(A=1,p=0.05)
						i<-genefilter(use.data.matrix_onlineNorm.m,fff)
						dat.fo<-use.data.matrix_onlineNorm.m[i,]
						i<-genefilter(-use.data.matrix_onlineNorm.m,fff)
						dat.fu<-use.data.matrix_onlineNorm.m[i,]
						dat.f<-unique(rbind(dat.fo,dat.fu))
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_O[[tonl_i]]<<-fit2
						data.matrix_onlineNorm.f<-fit2
					}),silent=TRUE)
					if(length(data.matrix_onlineNorm.f)==0)
					{
						try(({
							rsd<-rowSds(use.data.matrix_onlineNorm.m)
							i<-rsd>=2
							dat.f<-use.data.matrix_onlineNorm.m[i,]
							fit<-lmFit(dat.f)
							fit2<-eBayes(fit)
							fit2_O[[tonl_i]]<<-fit2
							data.matrix_onlineNorm.f<-fit2
						}),silent=TRUE)
					}
					ps1<-as.data.frame(data.matrix_onlineNorm.f)	
					rownames(ps1)<-rownames(data.matrix_onlineNorm.f)
					data.matrix_onlineNorm.f2=NULL
					data.matrix_onlineNorm.f2<-ps1
					ta_onl_filt=NULL
					ta_onl_filt<-tclArray()
					for(i in 0:dim(data.matrix_onlineNorm.f2)[1]){ for(j in 0:dim(data.matrix_onlineNorm.f2)[2]){ 
						if(i==0){ ta_onl_filt[[i,j]]<-colnames(data.matrix_onlineNorm.f2)[j] } else { 
						if(j==0){ ta_onl_filt[[i,j]]<-rownames(data.matrix_onlineNorm.f2)[i] } else {
						tem<-data.matrix_onlineNorm.f2[i,j]
						ta_onl_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
						} }
					} }
					
					data.matrix_onlineNorm.f3[[tonl_i]]<<-data.matrix_onlineNorm.f2
					
					ta_onl_filt2[[tonl_i]]<<-ta_onl_filt
	
p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
					if(length(data.matrix_onlineNorm.f)!=0){
						err<-try({
							ttx<-topTable(data.matrix_onlineNorm.f,number=nrow(use.data.matrix_onlineNorm.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_onlineNorm.m)==TRUE)]
							data.matrix_onlineNorm.s<-use.data.matrix_onlineNorm.m[rn,]
						},silent=TRUE)
						if(length(grep("Error",err))!=0)
						{
							ttx<-topTable(data.matrix_onlineNorm.f,number=nrow(use.data.matrix_onlineNorm.m))
							rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_onlineNorm.m)==TRUE)])
							data.matrix_onlineNorm.s<-data.matrix_onlineNorm.m[rownames(ttx)[ttx$P.Value<=0.01],]
						}
					}
					data.matrix_onlineNorm.s2=NULL
					data.matrix_onlineNorm.s2<-data.matrix_onlineNorm.s
					ta_onl_stat=NULL
					ta_onl_stat<-tclArray()
					for(i in 0:dim(data.matrix_onlineNorm.s2)[1]){ for(j in 0:dim(data.matrix_onlineNorm.s2)[2]){ 
						if(i==0){ ta_onl_stat[[i,j]]<-colnames(data.matrix_onlineNorm.s2)[j] } else { 
						if(j==0){ ta_onl_stat[[i,j]]<-rownames(data.matrix_onlineNorm.s2)[i] } else {
						tem<-data.matrix_onlineNorm.s2[i,j]
						ta_onl_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
						} }
					} }
					
					data.matrix_onlineNorm.s3[[tonl_i]]<<-data.matrix_onlineNorm.s2
					
					ta_onl_stat2[[tonl_i]]<<-ta_onl_stat
	
p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
					err<-try(DE_O<-topTable(data.matrix_onlineNorm.f),silent=TRUE)	
					err<-try(DE_OAll[[tonl_i]]<<-topTable(data.matrix_onlineNorm.f,number=nrow(data.matrix_onlineNorm.f)),silent=TRUE)
					ta_onl_dge=NULL
					ta_onl_dge<-tclArray()
					for(i in 0:dim(DE_O)[1]){ for(j in 0:dim(DE_O)[2]){ 
						if(i==0){ ta_onl_dge[[i,j]]<-colnames(DE_O)[j] } else { 
						if(j==0){ ta_onl_dge[[i,j]]<-rownames(DE_O)[i] } else {
						tem<-DE_O[i,j]
						ta_onl_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
						} }
					} }
					
					DE_O_2[[tonl_i]]<<-DE_O
					
					ta_onl_dge2[[tonl_i]]<<-ta_onl_dge
	
p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
					pca_O[[tonl_i]]<<-prcomp(t(data.matrix_onlineNorm.m))
					if(is.null(ma_img)==TRUE){
						try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
						try(tkwm.withdraw(w_legend),silent=TRUE)
						ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_O[[tonl_i]])},hscale=1.7,vscale=1.07)
						ma_img2<<-ma_img
						tkpack(ma_img2)			
					} else {
						ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_O[[tonl_i]])},hscale=1.7,vscale=1.07)
					}
	
p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
					sample.dist_c1<-dist(t(data.matrix_onlineNorm.m)) 
					sample.dist_O<-as.dist(1-cor(data.matrix_onlineNorm.m,method="pearson"))
					sample.clust_O[[tonl_i]]<<-hclust(sample.dist_O,method="complete")
					if(is.null(ma_img)==TRUE){
						try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
						try(tkwm.withdraw(w_legend),silent=TRUE)
						ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_O[[tonl_i]])},hscale=1.7,vscale=1.07)
						ma_img2<<-ma_img
						tkpack(ma_img2)			
					}else{
						ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_O[[tonl_i]])},hscale=1.7,vscale=1.07)
					}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
					er_x<-try(use_DE_O<-as.matrix(DE_O),silent=TRUE)
					if(length(grep("Error",er_x))!=0)
					{
						rownames(DE_O)<-as.character(DE_O)
						use_DE_O<-as.matrix(DE_O)
					}
					rnames<-rownames(use_DE_O)
					clas_O[[tonl_i]]<<-as.matrix(data.matrix_onlineNorm.m[rnames,])
					if(is.null(ma_img)==TRUE){
						try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
						try(tkwm.withdraw(w_legend),silent=TRUE)
						ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
							heatmap(clas_O[[tonl_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
						},hscale=1.7,vscale=1.07)
						ma_img2<<-ma_img
						tkpack(ma_img2)			
					}else{
						ma_img2<<-tkrreplot(ma_img,function(...){
							heatmap(clas_O[[tonl_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
						},hscale=1.7,vscale=1.07)
					}
p_enter_class<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Classification"))
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
				}
			}
		}
		
	})
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
if(!requireNamespace("limma",quietly=TRUE)){BiocManager::install("limma",ask=FALSE,update=FALSE);library(limma)} else {library(limma)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("GEOquery",quietly=TRUE)){BiocManager::install("GEOquery",ask=FALSE,update=FALSE);library(GEOquery)} else {library(GEOquery)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("impute",quietly=TRUE)){BiocManager::install("impute",ask=FALSE,update=FALSE);library(impute)} else {library(impute)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("stats",quietly=TRUE)){install.packages("stats",ask=FALSE,update=FALSE);library(stats)} else {library(stats)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk",ask=FALSE,update=FALSE);library(tcltk)} else {library(tcltk)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
dmo()
