gsea_gobp<-function(h,...){

	gsea_goBP_params<-function(h,h2,h3,h4,h5,h6,...){
print(h,h2)
		run_params<-function(h,h2,h3,h4,h5,h6,h7,h8,...){
			ts_name=h;ts_val=h2;ann_pk<-h3;db<-h4;stat_sig<-h5;p_v<-h6;children=h7;r_choice=h8;
pb_ma<-tkProgressBar(title="GSEA GO Biological Process...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			entrez_id<-paste(ann_pk,"ENTREZID",sep="")
			allg<-get(entrez_id)
			allg<-as.data.frame(unlist(as.list(allg)))
#			colnames(allg)<-"gene_id"
			myids<-unique(allg[rownames(stat_sig),])
			params<-new("GOHyperGParams",geneIds=myids,annotation=c(db),ontology="BP",pvalueCutoff=p_v,conditional=FALSE,testDirection="over")
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			GOresultBP<-hyperGTest(params)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
			res<-summary(GOresultBP)
			rownames(res)<-res[,1]
			res2<-res
			ids<-geneIdsByCategory(GOresultBP,catids=NULL)
			y<-allg
			y[,2]<-rownames(y)
			colnames(y)<-c("gene_id","probe_id")
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			for(i in 1:dim(res2)[1]){
				x<-c()
				for(j in 1:length(ids[[rownames(res2)[i]]])){
					z<-rownames(y[which((y$gene_id==ids[[rownames(res)[i]]][j])==TRUE),])
					z2<-c()
					for(k in 1:length(z)){
						z2<-paste(z[k],z2,sep=",")
					}
					z3<-paste(ids[[rownames(res)[i]]][j],paste("(",")",sep=z2),sep="")
					x<-paste(z3,x,sep=", ")
				}
				res2$genes[i]<-x
			}
			res2<-res2[,-1]
			res2<-res2[res2$Count>=3,]
			
			if(dim(res2)[1]==0)tkmessageBox(title="Warning",message="No results met the specified criteria",icon="warning",type="ok")
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			gsea_gobp<-tclArray()
			for(i in 0:dim(res2)[1]){ for(j in 0:dim(res2)[2]){ 
				if(i==0){ gsea_gobp[[i,j]]<-colnames(res2)[j] } else { 
				if(j==0){ gsea_gobp[[i,j]]<-rownames(res2)[i] } else {
				tem<-res2[i,j]
				gsea_gobp[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
				} }
			} }
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(ts_name=="Affymetrix"){
				GOresultBP_Affy[[ts_val]]<<-res2
				ta_affy_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Agilent-OneColor"){
				GOresultBP_Ag1[[ts_val]]<<-res2
				ta_ag1_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Agilent-TwoColor"){
				GOresultBP_Ag2[[ts_val]]<<-res2
				ta_ag2_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Illumina-Beadarray"){
				GOresultBP_Il_B[[ts_val]]<<-res2
				ta_ilb_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Illumina-Lumi"){
				GOresultBP_Il_L[[ts_val]]<<-res2
				ta_ill_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Nimblegen"){
				GOresultBP_N[[ts_val]]<<-res2
				ta_nbl_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Series-Matrix"){
				GOresultBP_S[[ts_val]]<<-res2
				ta_smt_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
			if(ts_name=="Online-Data"){
				GOresultBP_O[[ts_val]]<<-res2
				ta_onl_gsea_gobp[[ts_val]]<<-gsea_gobp
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOBP"))}
				}
			}
		}

		ts_name=h;ts_val=h2;pre_rs=h3;p_value_s=h4;children=h5;r_choice=h6;
		if(ts_name=="Affymetrix"){
			anAffy<-ann_Affy[[ts_val]]
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
				run_params(ts_name,ts_val,anAffy,db,DE_Affy_2[[ts_val]],p_value_s,children,r_choice)
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
				run_params(ts_name,ts_val,anAg1,db,DE_Ag1_2[[ts_val]],p_value_s,children,r_choice)
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
				run_params(ts_name,ts_val,anAg2,db,DE_Ag2_2[[ts_val]],p_value_s,children,r_choice)
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
				run_params(ts_name,ts_val,anIl_B,db,DE_Il_B_2[[ts_val]],p_value_s,children,r_choice)
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
				run_params(ts_name,ts_val,anIl_L,db,DE_Il_L_2[[ts_val]],p_value_s,children,r_choice)
			}
		}
		if(ts_name=="Nimblegen"){
			anNbl<-ann_N[[ts_val]]
			if(length(grep("GPL[0-9]",anNbl))==0){
				ann_pk<-anNbl
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anNbl),2]
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
				run_params(ts_name,ts_val,anNbl,db,DE_N_2[[ts_val]],p_value_s,children,r_choice)
			}
		}
		if(ts_name=="Series-Matrix"){
			anSmt<-ann_S[[ts_val]]
			if(length(grep("GPL[0-9]",anSmt))==0){
				ann_pk<-anSmt
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anSmt),2]
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
				run_params(ts_name,ts_val,anSmt,db,DE_S_2[[ts_val]],p_value_s,children,r_choice)
			}
		}
		if(ts_name=="Online-Data"){
			anOnl<-ann_O[[ts_val]]
			if(length(grep("GPL[0-9]",anOnl))==0){
				ann_pk<-anOnl
			} else {
				ann_pk<-pre_rs[which(pre_rs[,1]==anOnl),2]
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
				run_params(ts_name,ts_val,anOnl,db,DE_O_2[[ts_val]],p_value_s,children,r_choice)
			}
		}
	}

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
		tree_slc<-strsplit(r_choice[1],"_")
		ts_name<-tree_slc[[1]][1]
		ts_val<-as.numeric(tree_slc[[1]][2])
		w_imp<-tktoplevel()
		tkwm.title(w_imp,"Select P-value")
		tkwm.resizable(w_imp,0,0)
		frame_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame_imp,text='P-value'),side="left",padx=12)
		p_list<-c(1,0.5,0.1,0.05,0.01,0.001,0.0001,0.00001,0.0000000001)
		p_value<-tclVar(p_list[4])
		p_comb<-ttkcombobox(frame_imp,values=p_list,textvariable=p_value,state="normal",justify="left",width=7)
		tkpack(p_comb,side="left",padx=33)
		tkpack(frame_imp,expand=TRUE,fill="both",side="top")
		frame2_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		r_but<-ttkbutton(frame2_imp,text="Ok")
		tkpack(r_but,side="right",padx=12)
		tkpack(frame2_imp,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but,command=function(){
			tkwm.withdraw(w_imp)
			p_value_s<-as.numeric(tclvalue(p_value))
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
			gsea_goBP_params(ts_name,ts_val,pre_rs,p_value_s,children,r_choice)
		})
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
if(!requireNamespace("annotate",quietly=TRUE)){BiocManager::install("annotate",ask=FALSE,update=FALSE);library(annotate)} else {library(annotate)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("GEOmetadb",quietly=TRUE)){BiocManager::install("GEOmetadb",ask=FALSE,update=FALSE);library(GEOmetadb)} else {library(GEOmetadb)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("GOstats",quietly=TRUE)){BiocManager::install("GOstats",ask=FALSE,update=FALSE);library(GOstats)} else {library(GOstats)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("Category",quietly=TRUE)){BiocManager::install("Category",ask=FALSE,update=FALSE);library(Category)} else {library(Category)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("RSQLite",quietly=TRUE)){install.packages("RSQLite",ask=FALSE,update=FALSE);library(RSQLite)} else {library(RSQLite)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("DBI",quietly=TRUE)){install.packages("DBI",ask=FALSE,update=FALSE);library(DBI)} else {library(DBI)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
gsea_gobp()
