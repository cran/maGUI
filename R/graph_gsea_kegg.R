graph_gsea_kegg<-function(h,...){
	try(({KEGGresult_Affy<-KEGGresult_Affy;KEGGresult_Ag1<-KEGGresult_Ag1;KEGGresult_Ag2<-KEGGresult_Ag2;
		KEGGresult_Il_B<-KEGGresult_Il_B;KEGGresult_Il_L<-KEGGresult_Il_L;KEGGresult_N<-KEGGresult_N;
		KEGGresult_S<-KEGGresult_S;KEGGresult_O<-KEGGresult_O;
		ann_Affy<-ann_Affy;ann_Ag1<-ann_Ag1;ann_Ag2<-ann_Ag2;ann_Il_B<-ann_Il_B;ann_Il_L<-ann_Il_L;
		ann_N<-ann_N;ann_S<-ann_S;ann_O<-ann_O;
		DE_Affy2<-DE_Affy2;DE_Ag1_2<-DE_Ag1_2;DE_Ag2_2<-DE_Ag2_2;DE_Il_B2<-DE_Il_B2;DE_Il_L2<-DE_Il_L2;
		DE_N2<-DE_N2;DE_S2<-DE_S2;DE_O2<-DE_O2;
		g1_1<-g1_1;l<-l;tree<-tree;
	}),silent=TRUE)
	platforms=NULL;
	aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;
	try(({
		if(exists("KEGGresult_Affy"))aa=length(KEGGresult_Affy)
		if(exists("KEGGresult_Ag1"))bb=length(KEGGresult_Ag1)
		if(exists("KEGGresult_Ag2"))cc=length(KEGGresult_Ag2)
		if(exists("KEGGresult_Il_B"))dd=length(KEGGresult_Il_B)
		if(exists("KEGGresult_Il_L"))ee=length(KEGGresult_Il_L)
		if(exists("KEGGresult_N"))ff=length(KEGGresult_N)
		if(exists("KEGGresult_S"))gg=length(KEGGresult_S)
		if(exists("KEGGresult_O"))hh=length(KEGGresult_O)
		}),silent=TRUE)
	if(aa!=0)platforms=c(platforms,"Affymetrix")
	if(bb!=0)platforms=c(platforms,"Agilent_OneColor")
	if(cc!=0)platforms=c(platforms,"Agilent_TwoColor")
	if(dd!=0)platforms=c(platforms,"Illumina_Beadarray")
	if(ee!=0)platforms=c(platforms,"Illumina_Lumi")
	if(ff!=0)platforms=c(platforms,"Nimblegen")
	if(gg!=0)platforms=c(platforms,"Series_Matrix")
	if(hh!=0)platforms=c(platforms,"Online_Data")

	graph_gsea_kegg_Affy=NULL;graph_gsea_kegg_Ag1=NULL;graph_gsea_kegg_Ag2=NULL;graph_gsea_kegg_Il_B=NULL;graph_gsea_kegg_Il_L=NULL;
	graph_gsea_kegg_N=NULL;graph_gsea_kegg_S=NULL;graph_gsea_kegg_O=NULL;
	mapkG_gsea_kegg_Affy=NULL;logcol_gsea_kegg_Affy=NULL;mapkG_gsea_kegg_Ag1=NULL;logcol_gsea_kegg_Ag1=NULL;mapkG_gsea_kegg_Ag2=NULL;logcol_gsea_kegg_Ag2=NULL;
	mapkG_gsea_kegg_Il_B=NULL;logcol_gsea_kegg_Il_B=NULL;mapkG_gsea_kegg_Il_L=NULL;logcol_gsea_kegg_Il_L=NULL;mapkG_gsea_kegg_N=NULL;logcol_gsea_kegg_N=NULL;
	mapkG_gsea_kegg_S=NULL;logcol_gsea_kegg_S=NULL;mapkG_gsea_kegg_O=NULL;logcol_gsea_kegg_O=NULL;
	legend_gsea_kegg_Affy=NULL;legend_gsea_kegg_Ag1=NULL;legend_gsea_kegg_Ag2=NULL;legend_gsea_kegg_Il_B=NULL;legend_gsea_kegg_Il_L=NULL;
	legend_gsea_kegg_N=NULL;legend_gsea_kegg_S=NULL;legend_gsea_kegg_O=NULL;
	p_v=NULL;view_ww=NULL;mapkG=NULL;

	rm(graph_gsea_kegg_Affy,graph_gsea_kegg_Ag1,graph_gsea_kegg_Ag2,graph_gsea_kegg_Il_B,graph_gsea_kegg_Il_L,
	graph_gsea_kegg_N,graph_gsea_kegg_S,graph_gsea_kegg_O,
	mapkG_gsea_kegg_Affy,logcol_gsea_kegg_Affy,mapkG_gsea_kegg_Ag1,logcol_gsea_kegg_Ag1,mapkG_gsea_kegg_Ag2,logcol_gsea_kegg_Ag2,
	mapkG_gsea_kegg_Il_B,logcol_gsea_kegg_Il_B,mapkG_gsea_kegg_Il_L,logcol_gsea_kegg_Il_L,mapkG_gsea_kegg_N,logcol_gsea_kegg_N,
	mapkG_gsea_kegg_S,logcol_gsea_kegg_S,mapkG_gsea_kegg_O,logcol_gsea_kegg_O,
	legend_gsea_kegg_Affy,legend_gsea_kegg_Ag1,legend_gsea_kegg_Ag2,legend_gsea_kegg_Il_B,legend_gsea_kegg_Il_L,
	legend_gsea_kegg_N,legend_gsea_kegg_S,legend_gsea_kegg_O,p_v,view_ww,mapkG)

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
	x=NULL
	f<-function(h,...){
		x<<-svalue(h$obj)
		}
	z=NULL
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
#		svalue(sb)<-"				Please wait while Graphs.."
		w_dge<-gwindow("Select your data",width=260,height=280,visible=FALSE,horizontal=FALSE)
		gp_dge<-ggroup(container=w_dge,horizontal=FALSE)
		cbg_dge<-gcheckboxgroup(platforms,container=gp_dge,handler=f)
		svalue(cbg_dge,index=FALSE)<-1:8
		gp2_dge<-ggroup(container=gp_dge,width=30,height=15,horizontal=TRUE)
		addSpring(gp2_dge)
		y<-gbutton("CANCEL",border=TRUE,handler=function(h,...){
			dispose(w_dge)
		},container=gp2_dge,anchor=c(1,-1))
		y2<-gbutton("OK",border=TRUE,handler=function(h,...){
			svalue(sb)<-"				Please wait while Graphs.."
			if(length(x)!=0){
				dispose(w_dge)
				if(length(which(x=="Affymetrix"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_Affy)[KEGGresult_Affy$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						db<-"graph"
						try(library(db,character.only=TRUE),silent=TRUE)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_Affy),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_Affy2)
						entrez_id<-paste(ann_Affy,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_Affy2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_Affy<<-mapkG
						logcol_gsea_kegg_Affy<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_Affy<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_Affy<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_Affy)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_Affy,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_Affy)!=0){
							visible(g1_1)<-FALSE
							l$Affymetrix$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Agilent_OneColor"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_Ag1)[KEGGresult_Ag1$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_Ag1),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_Ag1_2)
						entrez_id<-paste(ann_Ag1,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_Ag1_2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_Ag1<<-mapkG
						logcol_gsea_kegg_Ag1<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_Ag1<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_Ag1<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_Ag1)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_Ag1,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_Ag1)!=0){
							visible(g1_1)<-FALSE
							l$Agilent_OneColor$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Agilent_TwoColor"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_Ag2)[KEGGresult_Ag2$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_Ag2),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_Ag2_2)
						entrez_id<-paste(ann_Ag2,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_Ag2_2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_Ag2<<-mapkG
						logcol_gsea_kegg_Ag2<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_Ag2<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_Ag2<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_Ag2)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_Ag2,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_Ag2)!=0){
							visible(g1_1)<-FALSE
							l$Agilent_TwoColor$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
					display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Illumina_Beadarray"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_Il_B)[KEGGresult_Il_B$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_Il_B),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_Il_B2)
						entrez_id<-paste(ann_Il_B,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_Il_B2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_Il_B<<-mapkG
						logcol_gsea_kegg_Il_B<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_Il_B<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_Il_B<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_Il_B)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_Il_B,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_Il_B)!=0){
							visible(g1_1)<-FALSE
							l$Illumina_Beadarray$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Illumina_Lumi"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_Il_L)[KEGGresult_Il_L$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_Il_L),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_Il_L2)
						entrez_id<-paste(ann_Il_L,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_Il_L2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_Il_L<<-mapkG
						logcol_gsea_kegg_Il_L<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_Il_L<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_Il_L<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_Il_L)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_Il_L,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_Il_L)!=0){
							visible(g1_1)<-FALSE
							l$Illumina_Lumi$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Nimblegen"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_N)[KEGGresult_N$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_N),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_N2)
						entrez_id<-paste(ann_N,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_N2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_N<<-mapkG
						logcol_gsea_kegg_N<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_N<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_N<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_N)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_N,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_N)!=0){
							visible(g1_1)<-FALSE
							l$Nimblegen$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Series_Matrix"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_S)[KEGGresult_S$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_S),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_S2)
						entrez_id<-paste(ann_S,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_S2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_S<<-mapkG
						logcol_gsea_kegg_S<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_S<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_S<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_S)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_S,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_S)!=0){
							visible(g1_1)<-FALSE
							l$Series_Matrix$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				if(length(which(x=="Online_Data"))!=0)
				{
					try(dispose(view_ww),silent=TRUE)
					KEGG.vec<-rownames(KEGGresult_O)[KEGGresult_O$Pvalue<=p_v]
				
					z=NULL
					w_kegg<-gwindow("Select KEGG ID",horizontal=FALSE,height=100,width=100)
					gp_kegg<-ggroup(container=w_kegg,horizontal=FALSE)
					glabel("KEGG Id",container=gp_kegg)
					cb_kegg<-gcombobox(KEGG.vec,editable=TRUE,selected=1,container=gp_kegg,handler=function(h,...){
						z<-svalue(h$obj)
						kegg_id<-z
					}
					)
					
					gp_kegg2<-ggroup(container=w_kegg)
					gbutton("CANCEL",border=TRUE,handler=function(h,...){
						svalue(sb)<-"Done"
						dispose(w_kegg)
						},container=gp_kegg2,anchor=c(1,-1))
					gbutton("OK",border=TRUE,handler=function(h,...){
						z<-svalue(cb_kegg)
						kegg_id<-z
						dispose(w_kegg)
						org_name<-db_bioc_org[which(db_bioc_org$db_bioc==ann_O),]$db_org
						org_code<-db_code[grep(org_name,db_code$"Organism Name"),]$"Organism Code"
						tmp_kegg<-paste(tempfile(),"xml",sep=".")
						tmp_kegg2<-retrieveKGML(kegg_id,organism=org_code,destfile=tmp_kegg,method="wget",quiet=TRUE)
						err<-try(mapkG<<-parseKGML2Graph(tmp_kegg,expandGenes=TRUE,genesOnly=TRUE),silent=TRUE)
						try(
						({
							if(length(grep("Error",err))!=0)
							{
								mapkG<<-parseKGML2Graph(tmp_kegg2,expandGenes=TRUE,genesOnly=TRUE)
							}
						}),silent=TRUE)
						ids<-rownames(DE_O2)
						entrez_id<-paste(ann_O,"ENTREZID",sep="")
						allg<-get(entrez_id)
						allg<-as.data.frame(unlist(as.list(allg)))
						myids<-allg[ids,]
						ENTREZ<-data.frame(cbind(ids,myids))
						colnames(ENTREZ)<-c("PROBEID","ENTREZID")
						top.Cell<-cbind(DE_O2,ENTREZ=ENTREZ$ENTREZID)
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
						mapkG_gsea_kegg_O<<-mapkG
						logcol_gsea_kegg_O<<-logcol
						nA<-makeNodeAttrs(mapkG,fillcolor=logcol,width=10,height=1.2)
						par(mar=c(3,5,0,5),mgp=c(0,0,0))
						layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
						graph_gsea_kegg_O<<-plot(mapkG,"dot",nodeAttrs=nA)
						image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
						mtext("down-regulation",side=1,at=0,line=1)
						mtext("up-regulation",side=1,at=1,line=1)
						legend<-data.frame(cbind(logcol,logfcs))
						legend_gsea_kegg_O<<-data.frame(rownames(legend),legend)
						colnames(legend_gsea_kegg_O)<<-c("Node","Node Color","logFC")
						view_ww<<-gwindow("Graph_Legend",visible=FALSE,height=400,width=600)
						gtable(legend_gsea_kegg_O,container=view_ww)
						visible(view_ww)<<-TRUE
						if(length(graph_gsea_kegg_O)!=0){
							visible(g1_1)<-FALSE
							l$Online_Data$Graph_GSEA_KEGG<<-list()
							tr<<-gtree(offspring=tree,container=g1_1)
							size(tr)<-c(300,400)
							visible(g1_1)<-TRUE
						}
						display()	
					},container=gp_kegg2,anchor=c(1,-1))
				}
				svalue(sb)<-"Done"
				dispose(w_dge)
				}else{
				gmessage("Plz select the data for generating Co-expression network","Select Data")
			}
		},container=gp2_dge,anchor=c(1,-1)
		)
		visible(w_dge)<-TRUE
	},container=gp_gsea2,anchor=c(1,-1)
	)
}
