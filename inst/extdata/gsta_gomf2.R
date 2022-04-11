gsta_gomf<-function(h,...){

	gsta_goMF_params<-function(h,h2,h3,h4,h5,h6,h7,...){

		run_params<-function(h,h2,h3,h4,h5,h6,h7,h8,h9,...){
			ts_name=h;ts_val=h2;ann_pk<-h3;db<-h4;stat_sig<-h5;ctls<-h6;tsts<-h7;children=h8;r_choice=h9;
pb_ma<-tkProgressBar(title="GSTA GO Molecular Function...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
			test.goMF<-gtGO(groups,t(stat_sig[,samp_gsta]),annotation=db,ontology="MF",sort=TRUE)
			res<-result(test.goMF)	
			names(res)[7]<-"Cov"
			res1<-res[which(res$Cov!=0),]
			res2<-res1
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
			for(i in 1:dim(res1)[1])
			{
				x<-c()
				for(j in 1:length(subsets(test.goMF)[[rownames(res1)[i]]]))
				{
					x<-paste(subsets(test.goMF)[[rownames(res1)[i]]][j],x,sep=",")
				}
				res2$genes[i]<-x
			}
			res2<-res2[,-1]
			if(dim(res2)[1]==0)tkmessageBox(title="Warning",message="No results met the specified criteria",icon="warning",type="ok")
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
			gsta_gomf<-tclArray()
			for(i in 0:dim(res2)[1]){ for(j in 0:dim(res2)[2]){ 
				if(i==0){ gsta_gomf[[i,j]]<-colnames(res2)[j] } else { 
				if(j==0){ gsta_gomf[[i,j]]<-rownames(res2)[i] } else {
				tem<-res2[i,j]
				gsta_gomf[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
				} }
			} }
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(ts_name=="Affymetrix"){
				GOtable.outMF_Affy[[ts_val]]<<-res2
				ta_affy_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Agilent-OneColor"){
				GOtable.outMF_Ag1[[ts_val]]<<-res2
				ta_ag1_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Agilent-TwoColor"){
				GOtable.outMF_Ag2[[ts_val]]<<-res2
				ta_ag2_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Illumina-Beadarray"){
				GOtable.outMF_Il_B[[ts_val]]<<-res2
				ta_ilb_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Illumina-Lumi"){
				GOtable.outMF_Il_L[[ts_val]]<<-res2
				ta_ill_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Nimblegen"){
				GOtable.outMF_N[[ts_val]]<<-res2
				ta_nbl_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Series-Matrix"){
				GOtable.outMF_S[[ts_val]]<<-res2
				ta_smt_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
			if(ts_name=="Online-Data"){
				GOtable.outMF_O[[ts_val]]<<-res2
				ta_onl_gsta_gomf[[ts_val]]<<-gsta_gomf
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSTA_GOMF"))}
				}
			}
		}

		ts_name=h;ts_val=h2;pre_rs=h3;samp_gsta=h4;groups=h5;children=h6;r_choice=h7;
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
				run_params(ts_name,ts_val,anAffy,db,dat2Affy.m2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anAg1,db,DE_Ag1_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anAg2,db,DE_Ag2_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anIl_B,db,DE_Il_B_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anIl_L,db,DE_Il_L_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anNbl,db,DE_N_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anSmt,db,DE_S_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
				run_params(ts_name,ts_val,anOnl,db,DE_O_2[[ts_val]],samp_gsta,groups,children,r_choice)
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
		tkwm.title(w_imp,"Select Samples")
		tkwm.resizable(w_imp,0,0)
		frame_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame_imp,text='Control Samples'),side="left",padx=45)
		tkpack(ttklabel(frame_imp,text='Test Samples'),side="right",padx=62)
		tkpack(frame_imp,expand=TRUE,fill="both",side="top")
		frame2_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		controls<-tclVar("")
		tkpack(ttkentry(frame2_imp,textvariable=controls,width=25),side="left",anchor="e")
		tests<-tclVar("")
		tkpack(ttkentry(frame2_imp,textvariable=tests,width=25),side="right",anchor="e")
		tkpack(frame2_imp,expand=TRUE,fill="both",side="top")
		frame3_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		r_but<-ttkbutton(frame3_imp,text="Ok")
		tkpack(r_but,side="right",padx=12)
		tkpack(frame3_imp,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but,command=function(){
			tkwm.withdraw(w_imp)
			ctls<-tclvalue(controls);tsts<-tclvalue(tests);
			xf<-strsplit(ctls,",")
			yf<-strsplit(tsts,",")
			samp_gsta<-c(xf[[1]],yf[[1]])
			groups<-c(rep(0,length(xf[[1]])),rep(1,length(yf[[1]])))
			
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
			gsta_goMF_params(ts_name,ts_val,pre_rs,samp_gsta,groups,children,r_choice)
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
if(!requireNamespace("globaltest",quietly=TRUE)){BiocManager::install("globaltest",ask=FALSE,update=FALSE);library(globaltest)} else {library(globaltest)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("RSQLite",quietly=TRUE)){install.packages("RSQLite",ask=FALSE,update=FALSE);library(RSQLite)} else {library(RSQLite)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("DBI",quietly=TRUE)){install.packages("DBI",ask=FALSE,update=FALSE);library(DBI)} else {library(DBI)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
gsta_gomf()
