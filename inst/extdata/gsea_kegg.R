gsea_kegg<-function(h,...){
pre_rs<-c("GPL32","GPL33","GPL34","GPL71","GPL72","GPL74","GPL75","GPL76","GPL77","GPL78","GPL79","GPL80","GPL81","GPL82","GPL83","GPL85","GPL86","GPL87","GPL88","GPL89","GPL90","GPL91","GPL92","GPL93","GPL94","GPL95","GPL96","GPL97","GPL98","GPL99","GPL100","GPL101","GPL198","GPL199","GPL200","GPL201","GPL339","GPL340","GPL341","GPL342","GPL570","GPL571","GPL886","GPL887","GPL1261","GPL1318","GPL1319","GPL1322","GPL1352","GPL1355","GPL1708","GPL2112","GPL2529","GPL2891","GPL2898","GPL3154","GPL3213","GPL3533","GPL3738","GPL3921","GPL3979","GPL4032","GPL4191","GPL5689","GPL6097","GPL6102","GPL6244","GPL6947","GPL8300","GPL8490","GPL10558","GPL11532","GPL13497","GPL13534","GPL13667","GPL15380","GPL15396","GPL17897","mgu74a","mgu74b","mgu74c","ag","drosgenome1","hcg110","mu11ksuba","mu11ksubb","mu19ksuba","mu19ksubb","mu19ksubc","hu6800","mgu74av2","mgu74bv2","mgu74cv2","rgu34a","rgu34b","rgu34c","rnu34","rtu34","ygs98","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","hgu133a","hgu133b","hu35ksuba","hu35ksubb","hu35ksubc","hu35ksubd","ath1121501","ecoli2","celegans","hgfocus","moe430a","mouse4302","rae230a","rae230b","hgu133plus2","hgu133a2","hgug4111a","hgug4110b","mouse430a2","xenopuslaevis","zebrafish","drosophila2","u133x3p","rat2302","hgug4112a","bovine","yeast2","h20kcod","adme16cod","ecoli2","chicken","porcine","canine2","hthgu133a","canine","","h10kcod","hgug4100a","illuminaHumanv1","illuminaHumanv2","hugene10sttranscriptcluster","illuminaHumanv3","hgu95av2","IlluminaHumanMethylation27k","illuminaHumanv4","hugene11sttranscriptcluster","HsAgilentDesign026652","IlluminaHumanMethylation450k","hgu219","GGHumanMethCancerPanelv1","hthgu133b","hthgu133a")

	rs<-NULL;
	con<-NULL;
	folder_Ann<-NULL;
	rm(rs,con,folder_Ann)
	loc=gconfirm("Annotation needs GEOmetadb.sqlite database",icon="info")
	if(loc==TRUE)
	{
#		choose_folder()
#		folderchoose=NULL;
		folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
#		folderchoose<<-folder
#		folder_Ann<<-folderchoose
		folder_Ann<<-folder
#		folderchoose=NULL
		if(length(folder_Ann)!=0)
		{
			setwd(folder_Ann)
			err=NULL;
			try(con<<-dbConnect(SQLite(),'GEOmetadb.sqlite'),silent=TRUE)
			try(geo_tables<-dbListTables(con),silent=TRUE)
			try(rs<<-dbGetQuery(con,'select gpl,bioc_package from gpl'),silent=TRUE)
			try(colnames(rs)<<-c("gpl","bioc_package"),silent=TRUE)
		}
	}	
	else
	{
		loc2<-gconfirm("Annotation from inbuilt database",icon="question")
		if(loc2==TRUE)
		{
			rs<<-matrix(pre_rs,ncol=2)
			try(colnames(rs)<<-c("gpl","bioc_package"),silent=TRUE)
		}
	}
	if(length(rs)!=0)
	{
		try(({dat2Affy.s<-dat2Affy.s;datAgOne2.s<-datAgOne2.s;datAgTwo2.s<-datAgTwo2.s;
		datIllBA2.s<-datIllBA2.s;lumi_NQ.s<-lumi_NQ.s;data.matrix_Nimblegen2.s<-data.matrix_Nimblegen2.s;
		data.matrixNorm.s<-data.matrixNorm.s;data.matrix_onlineNorm.s<-data.matrix_onlineNorm.s;g1_1<-g1_1;l<-l;tree<-tree;
			}),silent=TRUE)
	platforms=NULL
	aa=0;bb=0;cc=0;dd=0;ee=0;fff=0;gg=0;hh=0;
	try(({
		if(exists("dat2Affy.s"))aa=length(dat2Affy.s)
		if(exists("datAgOne2.s"))bb=length(datAgOne2.s)
		if(exists("datAgTwo2.s"))cc=length(datAgTwo2.s)
		if(exists("datIllBA2.s"))dd=length(datIllBA2.s)
		if(exists("lumi_NQ.s"))ee=length(lumi_NQ.s)
		if(exists("data.matrix_Nimblegen2.s"))fff=length(data.matrix_Nimblegen2.s)
		if(exists("data.matrixNorm.s"))gg=length(data.matrixNorm.s)
		if(exists("data.matrix_onlineNorm.s"))hh=length(data.matrix_onlineNorm.s)
		}),silent=TRUE)
	if(aa!=0)platforms=c(platforms,"Affymetrix")
	if(bb!=0)platforms=c(platforms,"Agilent_OneColor")
	if(cc!=0)platforms=c(platforms,"Agilent_TwoColor")
	if(dd!=0)platforms=c(platforms,"Illumina_Beadarray")
	if(ee!=0)platforms=c(platforms,"Illumina_Lumi")
	if(fff!=0)platforms=c(platforms,"Nimblegen")
	if(gg!=0)platforms=c(platforms,"Series_Matrix")
	if(hh!=0)platforms=c(platforms,"Online_Data")

	KEGGresult_Affy=NULL;KEGGresult_Ag1=NULL;KEGGresult_Ag2=NULL;KEGGresult_Il_B=NULL;KEGGresult_Il_L=NULL;
	KEGGresult_N=NULL;KEGGresult_S=NULL;KEGGresult_O=NULL;p_v=NULL;
	
	rm(KEGGresult_Affy,KEGGresult_Ag1,KEGGresult_Ag2,KEGGresult_Il_B,KEGGresult_Il_L,
	KEGGresult_N,KEGGresult_S,KEGGresult_O,p_v)

	x=NULL;
	f<-function(h,...){
		x<<-svalue(h$obj)
		}
	p_value=c(0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.5,1)
	w_gsea<-gwindow("Select p-value",horizontal=FALSE,height=100,width=100)
	gp_gsea<-ggroup(container=w_gsea,horizontal=FALSE)
	glabel("p-value",container=gp_gsea)
	cb_gsea<-gcombobox(p_value,editable=TRUE,selected=7,container=gp_gsea,handler=function(h,...){
		z<-svalue(h$obj)
		p_v<<-as.numeric(z)
	}
	)

	gp_gsea2<-ggroup(container=w_gsea)
	gbutton("CANCEL",border=TRUE,handler=function(h,...){
		p_v<<-0.05
		svalue(sb)<-"Done"
		dispose(w_gsea)
	},container=gp_gsea2,anchor=c(1,-1))
	gbutton("OK",border=TRUE,handler=function(h,...){
		z<-svalue(cb_gsea)
		p_v<<-as.numeric(z)
		dispose(w_gsea)
		svalue(sb)<-"				Please wait while KEGG pathways.."
		w_dge<-gwindow("Select your data",width=260,height=280,visible=FALSE,horizontal=FALSE)
		gp_dge<-ggroup(container=w_dge,horizontal=FALSE)
		cbg_dge<-gcheckboxgroup(platforms,container=gp_dge,handler=f)
		svalue(cbg_dge,index=FALSE)<-1:8
		gp2_dge<-ggroup(container=gp_dge,width=30,height=15,horizontal=TRUE)
		addSpring(gp2_dge)
		y<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
			svalue(sb)<-"Done"
			dispose(w_dge)
			},container=gp2_dge,anchor=c(1,-1))
		y2<-gbutton("OK",border=TRUE,handler=function(h,...){
			if(length(x)!=0){
			dispose(w_dge);
			if(length(which(x=="Affymetrix"))!=0)
			{
				if(length(ann_Affy)!=0 && is.na(ann_Affy)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_Affy)==TRUE)
						{
							ann_Affy<<-rs[which(rs[,1]==ann_Affy),2]
							if(length(ann_Affy)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Affy=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_Affy)==FALSE)
							{
								ann_Affy<<-ann_Affy
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_Affy<<-h$input
						
						if(grep("GPL[0-9]",ann_Affy)==TRUE)
						{
							ann_Affy<<-rs[which(rs[,1]==ann_Affy),2]
							
							if(length(ann_Affy)==0)
							{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Affy=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_Affy)!=0 && is.na(ann_Affy)==FALSE)
				{	
					db<-annPkgName(ann_Affy)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_Affy,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(dat2Affy.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_Affy<<-res2

					if(dim(summary(KEGGresult_Affy))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_Affy)!=0){
						l$Affymetrix$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
					}
					display()	
				}
			}
			if(length(which(x=="Agilent_OneColor"))!=0)
			{
				if(length(ann_Ag1)!=0 && is.na(ann_Ag1)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_Ag1)==TRUE)
						{
							ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
							if(length(ann_Ag1)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Ag1=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_Ag1)==FALSE)
							{
								ann_Ag1<<-ann_Ag1
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_Ag1<<-h$input
						if(grep("GPL[0-9]",ann_Ag1)==TRUE)
						{
							ann_Ag1<<-rs[which(rs[,1]==ann_Ag1),2]
							
							if(length(ann_Ag1)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Ag1=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_Ag1)!=0 && is.na(ann_Ag1)==FALSE)
				{	
					db<-annPkgName(ann_Ag1)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_Ag1,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(datAgOne2.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_Ag1<<-res2
					if(dim(summary(KEGGresult_Ag1))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_Ag1)!=0){
					l$Agilent_OneColor$GSEA_KEGG<<-list()
					tr<<-gtree(offspring=tree,container=g1_1)
					size(tr)<-c(300,480)
					}
					display()	
				}
			}
			if(length(which(x=="Agilent_TwoColor"))!=0)
			{
				if(length(ann_Ag2)!=0 && is.na(ann_Ag2)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_Ag2)==TRUE)
						{
							ann_Ag2<<-rs[which(rs[,1]==ann_Ag2),2]
							if(length(ann_Ag2)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Ag2=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_Ag2)==FALSE)
							{
								ann_Ag2<<-ann_Ag2
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_Ag2<<-h$input
						
						if(grep("GPL[0-9]",ann_Ag2)==TRUE)
						{
							ann_Ag2<<-rs[which(rs[,1]==ann_Ag2),2]
							
							if(length(ann_Ag2)==0)
							{
						gmessage("Annotation package not available for this platform",icon="warning")
						ann_Ag2=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_Ag2)!=0 && is.na(ann_Ag2)==FALSE)
				{	
					db<-annPkgName(ann_Ag2)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_Ag2,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(datAgTwo2.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
						z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_Ag2<<-res2
					if(dim(summary(KEGGresult_Ag2))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_Ag2)!=0){
						l$Agilent_TwoColor$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
						}
					display()	
				}
			}
			if(length(which(x=="Illumina_Beadarray"))!=0)
			{
				if(length(ann_Il_B)!=0 && is.na(ann_Il_B)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_Il_B)==TRUE)
						{
							ann_Il_B<<-rs[which(rs[,1]==ann_Il_B),2]
							if(length(ann_Il_B)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Il_B=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_Il_B)==FALSE)
							{
								ann_Il_B<<-ann_Il_B
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_Il_B<<-h$input
						if(grep("GPL[0-9]",ann_Il_B)==TRUE)
						{
							ann_Il_B<<-rs[which(rs[,1]==ann_Il_B),2]
							
							if(length(ann_Il_B)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Il_B=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_Il_B)!=0 && is.na(ann_Il_B)==FALSE)
				{	
					db<-annPkgName(ann_Il_B)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_Il_B,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(datIllBA2.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_Il_B<<-res2
					if(dim(summary(KEGGresult_Il_B))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_Il_B)!=0){
						l$Illumina_Beadarray$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
						}
					display()	
				}
			}
			if(length(which(x=="Illumina_Lumi"))!=0)
			{
				if(length(ann_Il_L)!=0 && is.na(ann_Il_L)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_Il_L)==TRUE)
						{
							ann_Il_L<<-rs[which(rs[,1]==ann_Il_L),2]
							if(length(ann_Il_L)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Il_L=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_Il_L)==FALSE)
							{
								ann_Il_L<<-ann_Il_L
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_Il_L<<-h$input
						if(grep("GPL[0-9]",ann_Il_L)==TRUE)
						{
							ann_Il_L<<-rs[which(rs[,1]==ann_Il_L),2]
							if(length(ann_Il_L)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_Il_L=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_Il_L)!=0 && is.na(ann_Il_L)==FALSE)
				{	
					db<-annPkgName(ann_Il_L)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_Il_L,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(lumi_NQ.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_Il_L<<-res2
					if(dim(summary(KEGGresult_Il_L))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_Il_L)!=0){
						l$Illumina_Lumi$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
						}
					display()	
				}
			}
			if(length(which(x=="Nimblegen"))!=0)
			{
				if(length(ann_N)!=0 && is.na(ann_N)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_N)==TRUE)
						{
							ann_N<<-rs[which(rs[,1]==ann_N),2]
							if(length(ann_N)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_N=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_N)==FALSE)
							{
								ann_N<<-ann_N
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_N<<-h$input
						if(grep("GPL[0-9]",ann_N)==TRUE)
						{
							ann_N<<-rs[which(rs[,1]==ann_N),2]
							if(length(ann_N)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_N=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_N)!=0 && is.na(ann_N)==FALSE)
				{	
					db<-annPkgName(ann_N)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_N,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(data.matrix_Nimblegen2.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
						z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_N<<-res2
					if(dim(summary(KEGGresult_N))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_N)!=0){
						l$Nimblegen$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
						}
					display()	
		  		}
	  		}
			if(length(which(x=="Series_Matrix"))!=0)
			{
				if(length(ann_S)!=0 && is.na(ann_S)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_S)==TRUE)
						{
							ann_S<<-rs[which(rs[,1]==ann_S),2]
							if(length(ann_S)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_S=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_S)==FALSE)
							{
								ann_S<<-ann_S
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_S<<-h$input
						if(grep("GPL[0-9]",ann_S)==TRUE)
						{
							ann_S<<-rs[which(rs[,1]==ann_S),2]
							if(length(ann_S)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_S=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_S)!=0 && is.na(ann_S)==FALSE)
				{	
					db<-annPkgName(ann_S)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_S,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(data.matrixNorm.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_S<<-res2
					if(dim(summary(KEGGresult_S))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_S)!=0){
						l$Series_Matrix$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
						}
					display()	
				}
			}
			if(length(which(x=="Online_Data"))!=0)
			{
				if(length(ann_O)!=0 && is.na(ann_O)==FALSE)
				{
					try(({
						if(grep("GPL[0-9]",ann_O)==TRUE)
						{
							ann_O<<-rs[which(rs[,1]==ann_O),2]
							if(length(ann_O)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_O=NULL
							}
						}
						else
						{
							if(grep("GPL[0-9]",ann_O)==FALSE)
							{
								ann_O<<-ann_O
							}
						}
					}),silent=TRUE)
				}
				else
				{
					ginput("Please provide GPL name",icon="question",handler=function(h,...)
					{
						ann_O<<-h$input
						if(grep("GPL[0-9]",ann_O)==TRUE)
						{
							ann_O<<-rs[which(rs[,1]==ann_O),2]
							if(length(ann_O)==0)
							{
								gmessage("Annotation package not available for this platform",icon="warning")
								ann_O=NULL
							}
						}
					}
					)	
				}	
				if(length(ann_O)!=0 && is.na(ann_O)==FALSE)
				{	
					db<-annPkgName(ann_O)
					err<-try(library(db,character.only=TRUE),silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						biocLite<-function()
						{
						    .Deprecated("BiocManager")
						}
						BiocInstaller<-function()
						{
						    .Deprecated("BiocManager")
						}
#						source("http://bioconductor.org/biocLite.R")
#						biocLite(db,dependencies=TRUE,suppressUpdates=TRUE)
						library(db,character.only=TRUE)
					}
					entrez_id<-paste(ann_O,"ENTREZID",sep="")
					allg<-get(entrez_id)
					allg<-as.data.frame(unlist(as.list(allg)))
					colnames(allg)<-"gene_id"
					y<-allg
					y[,2]<-rownames(y)
					colnames(y)<-c("gene_id","probe_id")
					myids<-unique(allg[rownames(data.matrix_onlineNorm.s),])
					params<-new("KEGGHyperGParams",geneIds=myids,annotation=db,pvalueCutoff=p_v,testDirection="over")
					KEGGresult<-hyperGTest(params)
					res<-summary(KEGGresult)
					rownames(res)<-res[,1]
					res2<-res
					ids<-geneIdsByCategory(KEGGresult,catids=NULL)
					for(i in 1:dim(res2)[1])
					{
						x<-c()
						for(j in 1:length(ids[[rownames(res2)[i]]]))
						{
							z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
							z2<-c()
							for(k in 1:length(z))
							{
								z2<-paste(z[k],z2,sep=",")
							}
							z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
							x<-paste(z3,x,sep=", ")
						}
						res2$genes[i]<-x
					}
					KEGGresult_O<<-res2
					if(dim(summary(KEGGresult_O))[1]==0)gmessage("No results met the specified criteria",icon="warning") 	
					if(length(KEGGresult_O)!=0){
						l$Online_Data$GSEA_KEGG<<-list()
						tr<<-gtree(offspring=tree,container=g1_1)
						size(tr)<-c(300,480)
					}
					display()	
				}
			}
			svalue(sb)<-"Done"
			dispose(w_dge)
		} else{
			gmessage("Plz select the data for GSEA","Select Data")
		}
		},container=gp2_dge,anchor=c(1,-1)
		)
			visible(w_dge)<-TRUE
	},container=gp_gsea2,anchor=c(1,-1)
	)
	}
}

if (!requireNamespace("annotate", quietly = TRUE))
{
	install.packages("annotate")
} else {
	library(annotate)
}
if (!requireNamespace("Category", quietly = TRUE))
{
	install.packages("Category")
} else {
	library(Category)
}
if (!requireNamespace("DBI", quietly = TRUE))
{
	install.packages("DBI")
} else {
	library(DBI)
}
if (!requireNamespace("GOstats", quietly = TRUE))
{
	install.packages("GOstats")
} else {
	library(GOstats)
}
if (!requireNamespace("RSQLite", quietly = TRUE))
{
	install.packages("RSQLite")
} else {
	library(RSQLite)
}
gsea_kegg()
