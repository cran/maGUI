graph_gsea_kegg<-function(h,...){

	graph_gsea_kegg_params<-function(h,h2,h3,h4,h5,h6,h7,...){
		ts_name=h;ts_val=h2;ann_pk=h3;kegg_id=h4;DE=h5;children=h6;r_choice=h7;
pb_ma<-tkProgressBar(title="Graph GSEA KEGG Pathway...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
		ar<-18
		cols<-colorRampPalette(c("Green", "Red"))(ar)
		db<-"graph"
		try(library(db,character.only=TRUE),silent=TRUE)
		org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_pk),]$db_org
		org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
		tmp_kegg<-paste(tempfile(),"xml",sep=".")
		tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
		err<-try(mapkG<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
		try({
			if(length(grep("Error",err))!=0){
				mapkG<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
			}
		},silent=TRUE)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
		ids<-rownames(DE)
		entrez_id<-paste(ann_pk,"ENTREZID",sep="")
		allg<-get(entrez_id)
		allg<-as.data.frame(unlist(as.list(allg)))
		myids<-allg[ids,]
		ENTREZ<-data.frame(cbind(ids,myids))
		colnames(ENTREZ)<-c("PROBEID","ENTREZID")
		top.Cell<-cbind(DE,ENTREZ=ENTREZ$ENTREZID)
		top.Cell<-top.Cell[which(!(is.na(top.Cell$ENTREZ))==TRUE),]
		BT_all<-as.character(top.Cell$ENTREZ)
		tg<-top.Cell[top.Cell$adj.P.Val<0.05,]
		BT_sig<-tg$logFC
		names(BT_sig)<-as.vector(tg$ENTREZ)
		isDiffExp<-gsub(paste(org_code,":",sep=""),"",nodes(mapkG)) %in% names(BT_sig)
		logfcs<-BT_sig[match(gsub(paste(org_code,":",sep=""),"",nodes(mapkG)),names(BT_sig))]
		names(logfcs)<-nodes(mapkG)
		logfcs[is.na(logfcs)]<-0
		incol<-round((logfcs+5)*3)
		incol[incol>ar]<-ar
		undetected<-!gsub(paste(org_code,":",sep=""),"",nodes(mapkG)) %in% BT_all
		logcol<-cols[incol]
		logcol[logfcs==0]<-"darkgrey"
		logcol[undetected]<-"white"
		names(logcol)<-names(logfcs)
		mapkG_gsea_kegg<-mapkG
		logcol_gsea_kegg<-logcol
		nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
		oldpar<-par(no.readonly=TRUE)
		on.exit(par(oldpar))
		par(mar=c(3,5,0,5),mgp=c(0,0,0))
		layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
		graph_gsea_kegg<-plot(mapkG,"dot",nodeAttrs=nA)
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
		image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
		mtext("down-regulation",side=1,at=0,line=1)
		mtext("up-regulation",side=1,at=1,line=1)
		legend<-data.frame(cbind(logcol,logfcs))
		legend_gsea_kegg<-legend
		colnames(legend_gsea_kegg)<-c("Node_Color","logFC")		
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
		ta_legend_gsea_kegg<-tclArray()
		for(i in 0:dim(legend_gsea_kegg)[1]){ for(j in 0:dim(legend_gsea_kegg)[2]){ 
			if(i==0){ ta_legend_gsea_kegg[[i,j]]<-colnames(legend_gsea_kegg)[j] } else { 
			if(j==0){ ta_legend_gsea_kegg[[i,j]]<-rownames(legend_gsea_kegg)[i] } else {
			tem<-legend_gsea_kegg[i,j]
			ta_legend_gsea_kegg[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
			} }
		} }
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
		if(ts_name=="Affymetrix"){
			graph_gsea_KEGG_Affy[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_Affy[[ts_val]]<<-logcol
			legend_gsea_KEGG_Affy[[ts_val]]<<-legend_gsea_kegg
			ta_affy_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_Affy[[ts_val]],nodeAttrs=nA)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_Affy[[ts_val]],nodeAttrs=nA)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_affy_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Affy[[ts_val]]),ncol(legend_gsea_KEGG_Affy[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Agilent-OneColor"){
			graph_gsea_KEGG_Ag1[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_Ag1[[ts_val]]<<-logcol
			legend_gsea_KEGG_Ag1[[ts_val]]<<-legend_gsea_kegg
			ta_ag1_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_Ag1[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_Ag1[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_ag1_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Ag1[[ts_val]]),ncol(legend_gsea_KEGG_Ag1[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Agilent-TwoColor"){
			graph_gsea_KEGG_Ag2[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_Ag2[[ts_val]]<<-logcol
			legend_gsea_KEGG_Ag2[[ts_val]]<<-legend_gsea_kegg
			ta_ag2_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_Ag2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_Ag2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_ag2_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Ag2[[ts_val]]),ncol(legend_gsea_KEGG_Ag2[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Illumina-Beadarray"){
			graph_gsea_KEGG_Il_B[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_Il_B[[ts_val]]<<-logcol
			legend_gsea_KEGG_Il_B[[ts_val]]<<-legend_gsea_kegg
			ta_ilb_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_Il_B[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_Il_B[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_ilb_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Il_B[[ts_val]]),ncol(legend_gsea_KEGG_Il_B[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Illumina-Lumi"){
			graph_gsea_KEGG_Il_L[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_Il_L[[ts_val]]<<-logcol
			legend_gsea_KEGG_Il_L[[ts_val]]<<-legend_gsea_kegg
			ta_ill_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_Il_L[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_Il_L[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_ill_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Il_L[[ts_val]]),ncol(legend_gsea_KEGG_Il_L[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Nimblegen"){
			graph_gsea_KEGG_N[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_N[[ts_val]]<<-logcol
			legend_gsea_KEGG_N[[ts_val]]<<-legend_gsea_kegg
			ta_nbl_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_N[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_N[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_nbl_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_N[[ts_val]]),ncol(legend_gsea_KEGG_N[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Series-Matrix"){
			graph_gsea_KEGG_S[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_S[[ts_val]]<<-logcol
			legend_gsea_KEGG_S[[ts_val]]<<-legend_gsea_kegg
			ta_smt_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_S[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_S[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_smt_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_S[[ts_val]]),ncol(legend_gsea_KEGG_S[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
		if(ts_name=="Online-Data"){
			graph_gsea_KEGG_O[[ts_val]]<<-mapkG
			nodefill_gsea_KEGG_O[[ts_val]]<<-logcol
			legend_gsea_KEGG_O[[ts_val]]<<-legend_gsea_kegg
			ta_onl_legend_gsea_kegg[[ts_val]]<<-ta_legend_gsea_kegg
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_KEGG_O[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_KEGG_O[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_KEGG_Graph"))
			}
			legend_tab(ta_onl_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_O[[ts_val]]),ncol(legend_gsea_KEGG_O[[ts_val]]),"GSEA_KEGG_Graph_Legend")

		}
	}
		
	db_bioc<-c("mgu74a","mgu74b","mgu74c","ag","drosgenome1","hcg110","mu11ksuba","mu11ksubb","mu19ksuba","mu19ksubb","mu19ksubc","hu6800","mgu74av2","mgu74bv2","mgu74cv2","rgu34a","rgu34b","rgu34c","rnu34","rtu34","ygs98","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","hgu133a","hgu133b","hu35ksuba","hu35ksubb","hu35ksubc","hu35ksubd","ath1121501","ecoli2","celegans","hgfocus","moe430a","mouse4302","rae230a","rae230b","hgu133plus2","hgu133a2","hgug4111a","hgug4110b","mouse430a2","xenopuslaevis","zebrafish","drosophila2","u133x3p","rat2302","hgug4112a","bovine","yeast2","yeast2","h20kcod","adme16cod","ecoli2","chicken","porcine","canine2","hthgu133a","canine","h10kcod","hgug4100a","illuminaHumanv1","illuminaHumanv2","hugene10sttranscriptcluster","illuminaHumanv3","hgu95av2","IlluminaHumanMethylation27k","illuminaHumanv4","hugene11sttranscriptcluster","HsAgilentDesign026652","IlluminaHumanMethylation450k","hgu219","GGHumanMethCancerPanelv1","hthgu133b","hthgu133a")
	db_org<-c("Mus musculus","Mus musculus","Mus musculus","Arabidopsis thaliana","Drosophila melanogaster","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Mus musculus","Homo sapiens","Mus musculus","Mus musculus","Mus musculus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Rattus norvegicus","Saccharomyces cerevisiae","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Arabidopsis thaliana","Escherichia coli K-12","Caenorhabditis elegans","Homo sapiens","Mus musculus","Mus musculus","Rattus norvegicus","Rattus norvegicus","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus","Xenopus laevis","Danio rerio","Drosophila melanogaster","Homo sapiens","Rattus norvegicus","Homo sapiens","Bos taurus","Schizosaccharomyces pombe","Saccharomyces cerevisiae","Homo sapiens","Rattus norvegicus","Escherichia coli","Gallus gallus","Sus scrofa","Canis lupus familiaris","Homo sapiens","Canis lupus familiaris","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens","Homo sapiens")
	kegg_info<-c("Mus musculus (mouse)","Arabidopsis thaliana (thale cress)","Drosophila melanogaster (fruit fly)","Wolbachia wMel (Drosophila melanogaster)","Homo sapiens (human)","Rattus norvegicus (rat)","Saccharomyces cerevisiae (budding yeast)","Escherichia coli K-12 MG1655","Escherichia coli K-12 W3110","Escherichia coli K-12 DH10B","Escherichia coli K-12 MDS42","Caenorhabditis elegans (nematode)","Xenopus laevis (African clawed frog)","Danio rerio (zebrafish)","Bos taurus (cow)","Schizosaccharomyces pombe (fission yeast)","Escherichia coli BW2952","Escherichia coli O157:H7 EDL933 (EHEC)","Escherichia coli O157:H7 Sakai (EHEC)","Escherichia coli O157:H7 EC4115 (EHEC)","Escherichia coli O157:H7 TW14359 (EHEC)","Escherichia coli O157:H7 Xuzhou21 (EHEC)","Escherichia coli O26:H11 11368 (EHEC)","Escherichia coli O111:H- 11128 (EHEC)","Escherichia coli O103:H2 12009 (EHEC)","Escherichia coli O127:H6 E2348/69 (EPEC)","Escherichia coli O55:H7 CB9615 (EPEC)","Escherichia coli O55:H7 RM12579 (EPEC)","Escherichia coli O6:K2:H1 CFT073 (UPEC)","Escherichia coli O6:K15:H31 536 (UPEC)","Escherichia coli O18:K1:H7 UTI89 (UPEC)","Escherichia coli APEC O1 (APEC)","Escherichia coli O9 HS (commensal)","Escherichia coli E24377A (ETEC)","Escherichia coli SMS-3-5 (environmental)","Escherichia coli O152:H28 SE11 (commensal)","Escherichia coli O8 IAI1 (commensal)","Escherichia coli O81 ED1a (commensal)","Escherichia coli 55989 (EAEC)","Escherichia coli O7:K1 IAI39 (ExPEC)","Escherichia coli O7:K1 CE10","Escherichia coli O17:K52:H18 UMN026 (ExPEC)","Escherichia coli O45:K1:H7 S88 (ExPEC)","Escherichia coli O44:H18 042 (EAEC)","Escherichia coli O83:H1 NRG 857C","Escherichia coli O78:H11:K80 H10407 (ETEC)","Escherichia coli O150:H5 SE15 (commensal)","Escherichia coli O104:H4 2009EL-2071","Escherichia coli O104:H4 2009EL-2050","Escherichia coli O104:H4 2011C-3493","Escherichia coli ATCC 8739","Escherichia coli B REL606","Escherichia coli BL21-Gold(DE3)pLysS AG","Escherichia coli KO11FL","Escherichia coli KO11FL","Escherichia coli ABU 83972","Escherichia coli DH1","Escherichia coli DH1","Escherichia coli IHE3034","Escherichia coli NA114 (UPEC)","Escherichia coli UM146","Escherichia coli UMNK88","Escherichia coli W","Escherichia coli W","Escherichia coli clone D i14","Escherichia coli clone D i2","Escherichia coli P12b","Escherichia coli BL21(DE3)","Escherichia coli BL21(DE3)","Escherichia coli LF82","Escherichia coli APEC O78","Escherichia coli LY180","Escherichia coli PMV-1","Escherichia coli JJ1886","Escherichia coli O145:H28 RM13514 (EHEC)","Escherichia coli O145:H28 RM13516 (EHEC)","Escherichia coli O25b:K100:H4-ST131 EC958 (UPEC)","Gallus gallus (chicken)","Sus scrofa (pig)","mmu","ath","dme","wol","hsa","rno","sce","eco","ecj","ecd","ecok","cel","xla","dre","bta","spo","ebw","ece","ecs",
"ecf","etw","elx","eoj","eoi","eoh","ecg","eok","elr","ecc","ecp","eci","ecv","ecx","ecw","ecm","ecy","ecr","ecq","eck","ect","eoc","eum","ecz","elo","eln",
"elh","ese","eso","esm","esl","ecl","ebr","ebd","eko","ekf","eab","edh","edj","eih","ena","elu","eun","elw","ell","elc","eld","elp","ebl","ebe","elf","ecoa",
"ecol","ecoi","ecoj","ecoo","ecoh","ecos","gga","ssc")
	db_code<-data.frame(kegg_info[1:79],kegg_info[80:158])
	colnames(db_code)<-c("Organism Name","Organism Code")	
	db_bioc_org<-data.frame(cbind(db_bioc,db_org))
	ar<-18
	cols<-colorRampPalette(c("Green", "Red"))(ar)

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
		p_list<-c(1,0.5,0.1,0.05,0.01,0.001,0.0001,0.00001,0.0000000001,0)
		p_value<-tclVar(p_list[4])
		p_comb<-ttkcombobox(frame_imp,values=p_list,textvariable=p_value,state="normal",justify="left",width=7)
		tkpack(p_comb,side="left",padx=33)
		tkpack(frame_imp,expand=TRUE,fill="both",side="top")
		frame2_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		r_but1<-ttkbutton(frame2_imp,text="Ok")
		tkpack(r_but1,side="right",padx=12)
		tkpack(frame2_imp,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but1,command=function(){
			tkwm.withdraw(w_imp)
			p_v<-as.numeric(tclvalue(p_value))
			if(ts_name=="Affymetrix"){
				KEGG.vec<-rownames(KEGGresult_Affy[[ts_val]])[KEGGresult_Affy[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_Affy[[ts_val]],kegg_id,DE_Affy_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Agilent-OneColor"){
				KEGG.vec<-rownames(KEGGresult_Ag1[[ts_val]])[KEGGresult_Ag1[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_Ag1[[ts_val]],kegg_id,DE_Ag1_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Agilent-TwoColor"){
				KEGG.vec<-rownames(KEGGresult_Ag2[[ts_val]])[KEGGresult_Ag2[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_Ag2[[ts_val]],kegg_id,DE_Ag2_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Illumina-Beadarray"){
				KEGG.vec<-rownames(KEGGresult_Il_B[[ts_val]])[KEGGresult_Il_B[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_Il_B[[ts_val]],kegg_id,DE_Il_B_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Illumina-Lumi"){
				KEGG.vec<-rownames(KEGGresult_Il_L[[ts_val]])[KEGGresult_Il_L[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_Il_L[[ts_val]],kegg_id,DE_Il_L_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Nimblegen"){
				KEGG.vec<-rownames(KEGGresult_N[[ts_val]])[KEGGresult_N[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_N[[ts_val]],kegg_id,DE_N_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Series-Matrix"){
				KEGG.vec<-rownames(KEGGresult_S[[ts_val]])[KEGGresult_S[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_S[[ts_val]],kegg_id,DE_S_2[[ts_val]],children,r_choice)
				})
			}
			if(ts_name=="Online-Data"){
				KEGG.vec<-rownames(KEGGresult_O[[ts_val]])[KEGGresult_O[[ts_val]]$Pvalue<=p_v]
				w_imp2<-tktoplevel()
				tkwm.title(w_imp2,"Select KEGG ID")
				tkwm.resizable(w_imp2,0,0)
				frame_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				tkpack(ttklabel(frame_imp2,text='KEGG ID'),side="left",padx=12)
				kids<-tclVar(KEGG.vec[1])
				k_comb<-ttkcombobox(frame_imp2,values=KEGG.vec,textvariable=kids,state="normal",justify="left",width=7)
				tkpack(k_comb,side="left",padx=33)
				tkpack(frame_imp2,expand=TRUE,fill="both",side="top")
				frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
				r_but2<-ttkbutton(frame2_imp2,text="Ok")
				tkpack(r_but2,side="right",padx=12)
				tkpack(frame2_imp2,expand=TRUE,fill="both",side="top")
				tkconfigure(r_but2,command=function(){
					tkwm.withdraw(w_imp2)
					kegg_id<-tclvalue(kids)
					graph_gsea_kegg_params(ts_name,ts_val,ann_O[[ts_val]],kegg_id,DE_O_2[[ts_val]],children,r_choice)
				})
			}	
		})
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("graph",quietly=TRUE)){ BiocManager::install("graph",ask=FALSE,update=FALSE);library(graph) }else{ library(graph) }
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("KEGGgraph",quietly=TRUE)){ BiocManager::install("KEGGgraph",ask=FALSE,update=FALSE);library(KEGGgraph) }else{ library(KEGGgraph) }
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
graph_gsea_kegg()

