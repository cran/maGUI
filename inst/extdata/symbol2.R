symbol<-function(h,...){
	r_choice="._."
	children<-as.character(tcl(treeview,"children",""))
	w_analz<-tktoplevel()
	w_analz_frame<-ttkframe(w_analz,padding=c(3,3,50,20),borderwidth=1,relief="groove")
	tkpack(w_analz_frame,expand=TRUE,fill="both",side="left")
	var<-tclVar(tree_prnts[1])
	sapply(tree_prnts,function(i){
		r_but_vals<-ttkradiobutton(w_analz_frame,variable=var,text=i,value=i)
		tkpack(r_but_vals,side="top",anchor="w")
	})
	tkpack(ttklabel(w_analz_frame,text=''))
	tkpack(ttklabel(w_analz_frame,text=''))
	r_but<-ttkbutton(w_analz_frame,text="Ok")
	tkpack(r_but,pady=2)
	tkconfigure(r_but,command=function(){
		r_choice<-tclvalue(var)
		tkwm.withdraw(w_analz)
		err<-try({
			if(!file.exists("GEOmetadb.sqlite")){
				sqlfile<-getSQLiteFile(destdir=cur_dir,destfile="GEOmetadb.sqlite.gz",type="normal")
			} else {
				sqlfile<-"GEOmetadb.sqlite"
			}
			con<-dbConnect(SQLite(),sqlfile)
			geo_tables<-dbListTables(con)
			rs<-dbGetQuery(con,'select gpl,bioc_package from gpl')
			rs2<-rs[!is.na(rs$bioc_package),]
			rs3<-rs2[!(rs2$bioc_package==""),]
			pre_rs<-rs3
		},silent=TRUE)
		if(length(grep("Error",err))!=0)
		{
			rs<-c("GPL71","GPL198","GPL2112","GPL3979","GPL3738","GPL200","GPL3213","GPL72","GPL1322","GPL199",
			"GPL3154","GPL74","GPL201","GPL96","GPL571","GPL97","GPL14877","GPL570","GPL13667","GPL8300",
			"GPL91","GPL92","GPL93","GPL94","GPL95","GPL887","GPL1708","GPL13497","GPL17897","GPL3921",
			"GPL15396","GPL80","GPL5188","GPL17556","GPL6244","GPL11532","GPL18190","GPL8490","GPL13534","GPL6102",
			"GPL6947","GPL10558","GPL32","GPL81","GPL33","GPL82","GPL34","GPL83","GPL339","GPL6246",
			"GPL340","GPL1261","GPL8321","GPL75","GPL76","GPL3533","GPL341","GPL342","GPL1355","GPL85",
			"GPL86","GPL87","GPL88","GPL1352","GPL1318","GPL2529","GPL90","GPL1319",
			"ag","ath1121501","bovine","canine","canine2","celegans","chicken","drosgenome1","drosophila2","ecoli2",
			"ecoli2","hcg110","hgfocus","hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu133plus2","hgu219","hgu95av2",
			"hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","hgug4110b","hgug4112a","HsAgilentDesign026652","hthgu133a","hthgu133a",
			"hthgu133b","hu6800","huex10sttranscriptcluster","hugene10sttranscriptcluster","hugene10sttranscriptcluster",
			"hugene11sttranscriptcluster","hugene11sttranscriptcluster","IlluminaHumanMethylation27k","IlluminaHumanMethylation450k","illuminaHumanv2",
			"illuminaHumanv3","illuminaHumanv4","mgu74a","mgu74av2","mgu74b","mgu74bv2","mgu74c","mgu74cv2","moe430a","mogene10sttranscriptcluster",
			"mouse4302","mouse430a2","mouse430a2","mu11ksuba","mu11ksubb","porcine","rae230a","rae230b","rat2302","rgu34a",
			"rgu34b","rgu34c","rnu34","u133x3p","xenopuslaevis","yeast2","ygs98","zebrafish")
			rs2<-matrix(rs,ncol=2)
			colnames(rs2)<-c("gpl","bioc_package")
			pre_rs<-data.frame(rs2)
		}
		tree_slc<-strsplit(r_choice[1],"_")
		ts_name<-tree_slc[[1]][1]
		ts_val<-as.numeric(tree_slc[[1]][2])
		if(ts_name=="Affymetrix"){
			anAffy<-ann_Affy[[ts_val]];print(anAffy)
			if(length(grep("GPL[0-9]",anAffy))==0){
				ann_pk<-anAffy
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anAffy),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_Affy[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(dat2Affy.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(dat2Affy.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_Affy[[ts_val]]<<-genes_t
				ta_affy_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Agilent-OneColor"){
			anAg1<-ann_Ag1[[ts_val]]
			if(length(grep("GPL[0-9]",anAg1))==0){
				ann_pk<-anAg1
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anAg1),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_Ag1[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(datAgOne2.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(datAgOne2.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_Ag1[[ts_val]]<<-genes_t
				ta_ag1_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Agilent-TwoColor"){
			anAg2<-ann_Ag2[[ts_val]]
			if(length(grep("GPL[0-9]",anAg2))==0){
				ann_pk<-anAg2
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anAg2),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_Ag2[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(datAgTwo2.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(datAgTwo2.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_Ag2[[ts_val]]<<-genes_t
				ta_ag2_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Illumina-Beadarray"){
			anIl_B<-ann_Il_B[[ts_val]]
			if(length(grep("GPL[0-9]",anIl_B))==0){
				ann_pk<-anIl_B
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anIl_B),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_Il_B[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(datIllBA2.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(datIllBA2.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_Il_B[[ts_val]]<<-genes_t
				ta_ilb_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Illumina-Lumi"){
			anIl_L<-ann_Il_L[[ts_val]]
			if(length(grep("GPL[0-9]",anIl_L))==0){
				ann_pk<-anIl_L
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anIl_L),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_Il_L[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(lumi_NQ.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(lumi_NQ.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_Il_L[[ts_val]]<<-genes_t
				ta_ill_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Nimblegen"){
			anN<-ann_N[[ts_val]]
			if(length(grep("GPL[0-9]",anN))==0){
				ann_pk<-anN
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anN),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_N[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(data.matrix_Nimblegen2.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(data.matrix_Nimblegen2.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_N[[ts_val]]<<-genes_t
				ta_nbl_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Series-Matrix"){
			anS<-ann_S[[ts_val]]
			if(length(grep("GPL[0-9]",anS))==0){
				ann_pk<-anS
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anS),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_S[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(data.matrixNorm.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(data.matrixNorm.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_S[[ts_val]]<<-genes_t
				ta_smt_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
				}
			}		
		}
		if(ts_name=="Online-Data"){
			anO<-ann_O[[ts_val]]
			if(length(grep("GPL[0-9]",anO))==0){
				ann_pk<-anO
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anO),2]
				if(length(ann_pk)==0){
					tkmessageBox(title="Error",message="Annotation package not available for this platform",icon="warning",type="ok")
				} else {
					ann_O[[ts_val]]<<-ann_pk
				}
			}
			if(ann_pk!=""){
				db<-annPkgName(ann_pk)
				if(!requireNamespace(db,quietly=TRUE)){
					BiocManager::install(db,ask=FALSE,update=FALSE)
				} else {
					library(db,character.only=TRUE)
				}
				Identifier<-rownames(data.matrix_onlineNorm.m2[[ts_val]])
				Symbol<-getSYMBOL(rownames(data.matrix_onlineNorm.m2[[ts_val]]),db)
				genes<-as.data.frame(Symbol)
				x<-as.data.frame(cbind(rownames(genes),genes))
				y<-x[which(is.na(x$Symbol)==FALSE),]
				genes_t<-as.data.frame(y$Symbol,rownames(y))
				colnames(genes_t)<-"Symbol"
				gn_sy<-tclArray()
				for(i in 0:dim(genes_t)[1]){ for(j in 0:dim(genes_t)[2]){ 
					if(i==0){ gn_sy[[i,j]]<-colnames(genes_t)[j] } else { 
					if(j==0){ gn_sy[[i,j]]<-rownames(genes_t)[i] } else {
					tem<-genes_t[i,j]
					gn_sy[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				genes_O[[ts_val]]<<-genes_t
				ta_onl_genes[[ts_val]]<<-gn_sy
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Gene_Symbol"))}
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
if(!requireNamespace("GEOmetadb",quietly=TRUE)){BiocManager::install("GEOmetadb",ask=FALSE,update=FALSE);library(GEOmetadb)} else {library(GEOmetadb)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("annotate",quietly=TRUE)){BiocManager::install("annotate",ask=FALSE,update=FALSE);library(annotate)} else {library(annotate)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("DBI",quietly=TRUE)){BiocManager::install("DBI",ask=FALSE,update=FALSE);library(DBI)} else {library(DBI)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("RSQLite",quietly=TRUE)){BiocManager::install("RSQLite",ask=FALSE,update=FALSE);library(RSQLite)} else {library(RSQLite)}
setTkProgressBar(pb_ma,1,title="Done...",label="")		
close(pb_ma)
symbol()

