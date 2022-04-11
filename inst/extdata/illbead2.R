illbead<-function(h,...){
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
	tilb_i<<-tilb_i+1
	folder_Il_B[[tilb_i]]<<-dirname(file)
	if(file==""){
		folder_Il_B<-NULL
	} else {
		try(({
pb_ma<-tkProgressBar(title="Pre-processing...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
			new_datIllBA=NULL
			setwd(folder_Il_B[[tilb_i]])
			probes<-c("ProbeID","TargetID","ProbeId","TargetId","Probe_ID","Probe_Id","Target_ID","Target_Id","ID_REF","Id_Ref")
			x<-read.table(file,sep="\t",fill=TRUE)
			y<-which(x[,4]!="")
			ToSkip<-y[1]-1
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			x2<-read.table(file,sep="\t",header=TRUE,skip=ToSkip,fill=TRUE)
			ProbeID=probes[which(probes %in% colnames(x2))]
#			colv1<-grep("X4068838021_H.AVG_Signal",colnames(x2))
#			list(exprs = "AVG_Signal", se.exprs="BEAD_STDERR",nObservations = "Avg_NBEADS", Detection="Detection Pval"),
			new_datIllBA<-readBeadSummaryData(file,sep="\t",skip=ToSkip,ProbeID=ProbeID,columns=list(exprs="AVG_Signal"))
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			if(length(new_datIllBA)!=0){datIllBA=NULL;datIllBA<-new_datIllBA}
			ann_Il_B[[tilb_i]]<<-annotation(datIllBA)
		}),silent=TRUE)
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
		if(length(new_datIllBA)==0)tkmessageBox(title="Error",message="Could not load this type of data! Try Series_matrix or Online_Data",icon="error",type="ok")
	 }
	if(length(datIllBA)!=0){
		if(length(new_datIllBA)!=0){

		prnt<-paste("Illumina-Beadarray",tilb_i,sep="_")
		tree_prnts<<-c(tree_prnts,prnt)
		l_datIllBA[[tilb_i]]<<-datIllBA
p_enter<-tcl(treeview,"insert","","end",values=as.tclObj(prnt))

			autopa<-tkmessageBox(title="Confirm",message="Do you want to pre-process and analyze automatically!",icon="question",type="okcancel")
			if(tclvalue(autopa)=="ok"){
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0,width=500)
				datIllBA2.m1=NULL;use.datIllBA2.m1=NULL;xf=NULL;datIllBA2.f=NULL;datIllBA2.s=NULL;
				design=NULL;design_Il_B=NULL;ttx=NULL;DE_Il_B=NULL;DE_Il_B_2=NULL;use_DE_Il_B=NULL;
				sample.dist_Il_B=NULL;
#				pca_Il_B=NULL;
#				sample.clust_Il_B=NULL;
#				clas_Il_B=NULL;
				try(({folder_Il_B[[tilb_i]]<-folder_Il_B[[tilb_i]];datIllBA<-datIllBA;heatcol<-heatcol;}),silent=TRUE)
				setwd(folder_Il_B[[tilb_i]])	
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
				datIllBA2<-normaliseIllumina(datIllBA,method="quantile")
				datIllBA2.m<-exprs(datIllBA2)
				datIllBA2.m1<-log2(datIllBA2.m)
				datIllBA2.m1<-as.matrix(datIllBA2.m1)
				try({
					if(length(rownames(datIllBA2.m1))==0){
						rownames(datIllBA2.m1)<-1:length(datIllBA2.m1[,1])
					} else {
						rownames(datIllBA2.m1)<-as.character(rownames(datIllBA2.m1))
					}
				},silent=TRUE)
				ta_ilb_norm=NULL
				ta_ilb_norm<-tclArray()
				for(i in 0:dim(datIllBA2.m1)[1]){ for(j in 0:dim(datIllBA2.m1)[2]){ 
					if(i==0){ ta_ilb_norm[[i,j]]<-colnames(datIllBA2.m1)[j] } else { 
					if(j==0){ ta_ilb_norm[[i,j]]<-rownames(datIllBA2.m1)[i] } else {
					tem<-datIllBA2.m1[i,j]
					ta_ilb_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }

#				datIllBA2.m1[[tilb_i]]<<-datIllBA2.m1[[tilb_i]]
				datIllBA2.m2[[tilb_i]]<<-datIllBA2.m1
				
				ta_ilb_norm2[[tilb_i]]<<-ta_ilb_norm
#				

p_enter_norm<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Normalization"))
setTkProgressBar(pb_ma,0.3,title="QC...",label="")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datIllBA2.m2[[tilb_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datIllBA2.m2[[tilb_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_qc<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...QC"))
setTkProgressBar(pb_ma,0.35,title="Filtering...",label="")
				use.datIllBA2.m1<-datIllBA2.m1
				try(({
					fff<-pOverA(A=1,p=0.05)
					i<-genefilter(use.datIllBA2.m1,fff)
					dat.fo<-use.datIllBA2.m1[i,]
					i<-genefilter(-use.datIllBA2.m1,fff)
					dat.fu<-use.datIllBA2.m1[i,]
					dat.f<-unique(rbind(dat.fo,dat.fu))
					fit<-lmFit(dat.f)
					fit2<-eBayes(fit)
					fit2_Il_B[[tilb_i]]<<-fit2
					datIllBA2.f<-fit2
				}),silent=TRUE)
				if(length(datIllBA2.f)==0)
				{
					try(({
						rsd<-rowSds(use.datIllBA2.m1)
						i<-rsd>=2
						dat.f<-use.datIllBA2.m1[i,]
						fit<-lmFit(dat.f)
						fit2<-eBayes(fit)
						fit2_Il_B[[tilb_i]]<<-fit2
						datIllBA2.f<-fit2
					}),silent=TRUE)
				}
				ps1<-as.data.frame(datIllBA2.f)	
				rownames(ps1)<-rownames(datIllBA2.f)
				datIllBA2.f2=NULL
				datIllBA2.f2<-ps1
				ta_ilb_filt=NULL
				ta_ilb_filt<-tclArray()
				for(i in 0:dim(datIllBA2.f2)[1]){ for(j in 0:dim(datIllBA2.f2)[2]){ 
					if(i==0){ ta_ilb_filt[[i,j]]<-colnames(datIllBA2.f2)[j] } else { 
					if(j==0){ ta_ilb_filt[[i,j]]<-rownames(datIllBA2.f2)[i] } else {
					tem<-datIllBA2.f2[i,j]
					ta_ilb_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				datIllBA2.f3[[tilb_i]]<<-datIllBA2.f2
				
				ta_ilb_filt2[[tilb_i]]<<-ta_ilb_filt

p_enter_filter<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Filtering"))
setTkProgressBar(pb_ma,0.55,title="Statistical_Significant...",label="")
				if(length(datIllBA2.f)!=0){
					err<-try({
							ttx<-topTable(datIllBA2.f,number=nrow(use.datIllBA2.m1))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datIllBA2.m1)==TRUE)]
							datIllBA2.s<-use.datIllBA2.m1[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datIllBA2.f,number=nrow(use.datIllBA2.m1))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datIllBA2.m1)==TRUE)])
						datIllBA2.s<-datIllBA2.m1[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datIllBA2.s2=NULL
				datIllBA2.s2<-datIllBA2.s
				ta_ilb_stat=NULL
				ta_ilb_stat<-tclArray()
				for(i in 0:dim(datIllBA2.s2)[1]){ for(j in 0:dim(datIllBA2.s2)[2]){ 
					if(i==0){ ta_ilb_stat[[i,j]]<-colnames(datIllBA2.s2)[j] } else { 
					if(j==0){ ta_ilb_stat[[i,j]]<-rownames(datIllBA2.s2)[i] } else {
					tem<-datIllBA2.s2[i,j]
					ta_ilb_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				datIllBA2.s3[[tilb_i]]<<-datIllBA2.s2
				
				ta_ilb_stat2[[tilb_i]]<<-ta_ilb_stat

p_enter_stat_sign<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Statistical_Significant"))
setTkProgressBar(pb_ma,0.75,title="DGE...",label="")
				err<-try(DE_Il_B<-topTable(datIllBA2.f),silent=TRUE)	
				err<-try(DE_Il_BAll[[tilb_i]]<<-topTable(datIllBA2.f,number=nrow(datIllBA2.f)),silent=TRUE)
				ta_ilb_dge=NULL
				ta_ilb_dge<-tclArray()
				for(i in 0:dim(DE_Il_B)[1]){ for(j in 0:dim(DE_Il_B)[2]){ 
					if(i==0){ ta_ilb_dge[[i,j]]<-colnames(DE_Il_B)[j] } else { 
					if(j==0){ ta_ilb_dge[[i,j]]<-rownames(DE_Il_B)[i] } else {
					tem<-DE_Il_B[i,j]
					ta_ilb_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				
				DE_Il_B_2[[tilb_i]]<<-DE_Il_B
				
				ta_ilb_dge2[[tilb_i]]<<-ta_ilb_dge

p_enter_dge<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...DGE"))
setTkProgressBar(pb_ma,0.85,title="PCA...",label="")
				pca_Il_B[[tilb_i]]<<-prcomp(t(datIllBA2.m1))
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(pca_Il_B[[tilb_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				} else {
					ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Il_B[[tilb_i]])},hscale=1.7,vscale=1.07)
				}

p_enter_pca<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...PCA"))
setTkProgressBar(pb_ma,0.9,title="Clustering...",label="")
				sample.dist_c1<-dist(t(datIllBA2.m1)) 
				sample.dist_Il_B<-as.dist(1-cor(datIllBA2.m1,method="pearson"))
				sample.clust_Il_B[[tilb_i]]<<-hclust(sample.dist_Il_B,method="complete")
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(sample.clust_Il_B[[tilb_i]])},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Il_B[[tilb_i]])},hscale=1.7,vscale=1.07)
				}
p_enter_clust<-tcl(treeview,"insert",p_enter,"end",values=as.tclObj("...Clustering"))
setTkProgressBar(pb_ma,0.95,title="Classification...",label="")
				er_x<-try(use_DE_Il_B<-as.matrix(DE_Il_B),silent=TRUE)
				if(length(grep("Error",er_x))!=0)
				{
					rownames(DE_Il_B)<-as.character(DE_Il_B)
					use_DE_Il_B<-as.matrix(DE_Il_B)
				}
				rnames<-rownames(use_DE_Il_B)
				clas_Il_B[[tilb_i]]<<-as.matrix(datIllBA2.m1[rnames,])
				if(is.null(ma_img)==TRUE){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
						heatmap(clas_Il_B[[tilb_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
					},hscale=1.7,vscale=1.07)
					ma_img2<<-ma_img
					tkpack(ma_img2)			
				}else{
					ma_img2<<-tkrreplot(ma_img,function(...){
						heatmap(clas_Il_B[[tilb_i]],Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
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
if(!requireNamespace("beadarray",quietly=TRUE)){BiocManager::install("beadarray",ask=FALSE,update=FALSE);library(beadarray)} else {library(beadarray)}
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
illbead()
