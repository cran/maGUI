display_mod<-function(h,...){
	treeview<-treeview
	tree_w<<-NULL
	tree_mod_n<<-1
	tree_mods<<-rbind(c(tree_mod_n,"",""))
	tkbind(treeview,"<Button-3>",function(W,x,y){
		try(tkdestroy(tree_w),silent=TRUE)
		mapapa<-tcl(treeview,"parent",tkselect(treeview))
		mapapa_name<-as.character(tcl(treeview,"item",mapapa,"-values"))
		betabeti_name<-as.character(tcl(treeview,"item",tkselect(treeview),"-values"))
		tree_detdel<-function(h,...){
			tree_w<<-tktoplevel()
			tree_frame<-ttkframe(tree_w)
			tkpack(tree_frame,expand=TRUE,fill="both")
			detdelbut_frame<-ttkframe(tree_frame)
			tkpack(detdelbut_frame,anchor="ne",fill="both")
			detbut<-ttkbutton(detdelbut_frame,text="Remove")
				tkpack(detbut,side="top")
			tkconfigure(detbut,command=function(){
				d_p_affy<-tcl(treeview,"detach",tkselect(treeview))
				tree_mod_n<<-tree_mod_n+1
				tree_mods<<-rbind(tree_mods,c(tree_mod_n,mapapa_name,betabeti_name))
				tkdestroy(tree_w)
			})
			delbut<-ttkbutton(detdelbut_frame,text="Delete")
			tkpack(delbut,side="bottom")
			tkconfigure(delbut,command=function(){
				d_p_affy<-tcl(treeview,"delete",tkselect(treeview))
				tree_slc<-strsplit(mapapa_name[1],"_")
				ts_val<-as.numeric(tree_slc[[1]][2])
				if(!is.na(ts_val) && mapapa_name==paste("Affymetrix",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_affy_norm2[[ts_val]]=NA; dat2Affy.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_affy_filt2[[ts_val]]=NA; dat2Affy.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_affy_stat2[[ts_val]]=NA; dat2Affy.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_affy_dge2[[ts_val]]=NA; DE_Affy_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_Affy[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_Affy[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_affy_gsea_gobp[[ts_val]]=NA;GOresultBP_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_affy_gsea_gomf[[ts_val]]=NA;GOresultMF_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_affy_gsea_gocc[[ts_val]]=NA;GOresultCC_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_affy_gsea_kegg[[ts_val]]=NA;KEGGresult_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_affy_gsta_gobp[[ts_val]]=NA;GOtable.outBP_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_affy_gsta_gomf[[ts_val]]=NA;GOtable.outMF_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_affy_gsta_gocc[[ts_val]]=NA;GOtable.outCC_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_affy_gsta_kegg[[ts_val]]=NA;KEGGtable.out_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_Affy[[ts_val]]=NA;ta_affy_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_Affy[[ts_val]]=NA;ta_affy_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_Affy[[ts_val]]=NA;ta_affy_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_Affy[[ts_val]]=NA;ta_affy_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_Affy[[ts_val]]=NA;ta_affy_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_Affy[[ts_val]]=NA;ta_affy_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_Affy[[ts_val]]=NA;ta_affy_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_Affy[[ts_val]]=NA;ta_affy_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_affy_genes[[ts_val]]=NA;genes_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...Coexprs_Network"){ myGraph_Affy[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_Affy[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Agilent-OneColor",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_ag1_norm2[[ts_val]]=NA;datAgOne2.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_ag1_filt2[[ts_val]]=NA; datAgOne2.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_ag1_stat2[[ts_val]]=NA; datAgOne2.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_ag1_dge2[[ts_val]]=NA; DE_Ag1_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_Ag1[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_Ag1[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_ag1_gsea_gobp[[ts_val]]=NA;GOresultBP_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_ag1_gsea_gomf[[ts_val]]=NA;GOresultMF_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_ag1_gsea_gocc[[ts_val]]=NA;GOresultCC_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_ag1_gsea_kegg[[ts_val]]=NA;KEGGresult_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_ag1_gsta_gobp[[ts_val]]=NA;GOtable.outBP_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_ag1_gsta_gomf[[ts_val]]=NA;GOtable.outMF_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_ag1_gsta_gocc[[ts_val]]=NA;GOtable.outCC_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_ag1_gsta_kegg[[ts_val]]=NA;KEGGtable.out_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_Ag1[[ts_val]]=NA;ta_ag1_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_Ag1[[ts_val]]=NA;ta_ag1_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_Ag1[[ts_val]]=NA;ta_ag1_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_Ag1[[ts_val]]=NA;ta_ag1_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_Ag1[[ts_val]]=NA;ta_ag1_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_Ag1[[ts_val]]=NA;ta_ag1_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_Ag1[[ts_val]]=NA;ta_ag1_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_Ag1[[ts_val]]=NA;ta_ag1_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_ag1_genes[[ts_val]]=NA;genes_Ag1[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_Ag1[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_Ag1[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Agilent-TwoColor",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_ag2_norm2[[ts_val]]=NA;datAgTwo2.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_ag2_filt2[[ts_val]]=NA; datAgTwo2.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_ag2_stat2[[ts_val]]=NA; datAgTwo2.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_ag2_dge2[[ts_val]]=NA; DE_Ag2_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_Ag2[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_Ag2[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_ag2_gsea_gobp[[ts_val]]=NA;GOresultBP_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_ag2_gsea_gomf[[ts_val]]=NA;GOresultMF_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_ag2_gsea_gocc[[ts_val]]=NA;GOresultCC_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_ag2_gsea_kegg[[ts_val]]=NA;KEGGresult_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_ag2_gsta_gobp[[ts_val]]=NA;GOtable.outBP_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_ag2_gsta_gomf[[ts_val]]=NA;GOtable.outMF_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_ag2_gsta_gocc[[ts_val]]=NA;GOtable.outCC_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_ag2_gsta_kegg[[ts_val]]=NA;KEGGtable.out_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_Ag2[[ts_val]]=NA;ta_ag2_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_Ag2[[ts_val]]=NA;ta_ag2_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_Ag2[[ts_val]]=NA;ta_ag2_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_Ag2[[ts_val]]=NA;ta_ag2_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_Ag2[[ts_val]]=NA;ta_ag2_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_Ag2[[ts_val]]=NA;ta_ag2_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_Ag2[[ts_val]]=NA;ta_ag2_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_Ag2[[ts_val]]=NA;ta_ag2_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_ag2_genes[[ts_val]]=NA;genes_Ag2[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_Ag2[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_Ag2[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Illumina-Beadarray",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_ilb_norm2[[ts_val]]=NA;datIllBA2.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_ilb_filt2[[ts_val]]=NA; datIllBA2.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_ilb_stat2[[ts_val]]=NA; datIllBA2.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_ilb_dge2[[ts_val]]=NA; DE_Il_B_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_Il_B[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_Il_B[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_ilb_gsea_gobp[[ts_val]]=NA;GOresultBP_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_ilb_gsea_gomf[[ts_val]]=NA;GOresultMF_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_ilb_gsea_gocc[[ts_val]]=NA;GOresultCC_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_ilb_gsea_kegg[[ts_val]]=NA;KEGGresult_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_ilb_gsta_gobp[[ts_val]]=NA;GOtable.outBP_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_ilb_gsta_gomf[[ts_val]]=NA;GOtable.outMF_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_ilb_gsta_gocc[[ts_val]]=NA;GOtable.outCC_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_ilb_gsta_kegg[[ts_val]]=NA;KEGGtable.out_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_Il_B[[ts_val]]=NA;ta_ilb_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_Il_B[[ts_val]]=NA;ta_ilb_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_Il_B[[ts_val]]=NA;ta_ilb_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_Il_B[[ts_val]]=NA;ta_ilb_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_Il_B[[ts_val]]=NA;ta_ilb_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_Il_B[[ts_val]]=NA;ta_ilb_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_Il_B[[ts_val]]=NA;ta_ilb_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_Il_B[[ts_val]]=NA;ta_ilb_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_ilb_genes[[ts_val]]=NA;genes_Il_B[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_Il_B[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_Il_B[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Illumina-Lumi",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_ill_norm2[[ts_val]]=NA;lumi_NQ.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_ill_filt2[[ts_val]]=NA; lumi_NQ.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_ill_stat2[[ts_val]]=NA; lumi_NQ.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_ill_dge2[[ts_val]]=NA; DE_Il_L_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_Il_L[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_Il_L[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_ill_gsea_gobp[[ts_val]]=NA;GOresultBP_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_ill_gsea_gomf[[ts_val]]=NA;GOresultMF_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_ill_gsea_gocc[[ts_val]]=NA;GOresultCC_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_ill_gsea_kegg[[ts_val]]=NA;KEGGresult_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_ill_gsta_gobp[[ts_val]]=NA;GOtable.outBP_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_ill_gsta_gomf[[ts_val]]=NA;GOtable.outMF_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_ill_gsta_gocc[[ts_val]]=NA;GOtable.outCC_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_ill_gsta_kegg[[ts_val]]=NA;KEGGtable.out_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_Il_L[[ts_val]]=NA;ta_ill_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_Il_L[[ts_val]]=NA;ta_ill_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_Il_L[[ts_val]]=NA;ta_ill_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_Il_L[[ts_val]]=NA;ta_ill_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_Il_L[[ts_val]]=NA;ta_ill_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_Il_L[[ts_val]]=NA;ta_ill_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_Il_L[[ts_val]]=NA;ta_ill_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_Il_L[[ts_val]]=NA;ta_ill_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_ill_genes[[ts_val]]=NA;genes_Il_L[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_Il_L[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_Il_L[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Nimblegen",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_nbl_norm2[[ts_val]]=NA;data.matrix_Nimblegen2.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_nbl_filt2[[ts_val]]=NA; data.matrix_Nimblegen2.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_nbl_stat2[[ts_val]]=NA; data.matrix_Nimblegen2.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_nbl_dge2[[ts_val]]=NA; DE_N_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_N[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_N[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_nbl_gsea_gobp[[ts_val]]=NA;GOresultBP_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_nbl_gsea_gomf[[ts_val]]=NA;GOresultMF_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_nbl_gsea_gocc[[ts_val]]=NA;GOresultCC_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_nbl_gsea_kegg[[ts_val]]=NA;KEGGresult_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_nbl_gsta_gobp[[ts_val]]=NA;GOtable.outBP_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_nbl_gsta_gomf[[ts_val]]=NA;GOtable.outMF_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_nbl_gsta_gocc[[ts_val]]=NA;GOtable.outCC_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_nbl_gsta_kegg[[ts_val]]=NA;KEGGtable.out_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_N[[ts_val]]=NA;ta_nbl_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_N[[ts_val]]=NA;ta_nbl_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_N[[ts_val]]=NA;ta_nbl_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_N[[ts_val]]=NA;ta_nbl_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_N[[ts_val]]=NA;ta_nbl_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_N[[ts_val]]=NA;ta_nbl_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_N[[ts_val]]=NA;ta_nbl_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_N[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_N[[ts_val]]=NA;ta_nbl_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_N[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_nbl_genes[[ts_val]]=NA;genes_N[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_N[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_N[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Series-Matrix",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_smt_norm2[[ts_val]]=NA;data.matrixNorm.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_smt_filt2[[ts_val]]=NA; data.matrixNorm.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_smt_stat2[[ts_val]]=NA; data.matrixNorm.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_smt_dge2[[ts_val]]=NA; DE_S_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_S[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_S[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_smt_gsea_gobp[[ts_val]]=NA;GOresultBP_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_smt_gsea_gomf[[ts_val]]=NA;GOresultMF_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_smt_gsea_gocc[[ts_val]]=NA;GOresultCC_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_smt_gsea_kegg[[ts_val]]=NA;KEGGresult_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_smt_gsta_gobp[[ts_val]]=NA;GOtable.outBP_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_smt_gsta_gomf[[ts_val]]=NA;GOtable.outMF_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_smt_gsta_gocc[[ts_val]]=NA;GOtable.outCC_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_smt_gsta_kegg[[ts_val]]=NA;KEGGtable.out_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_S[[ts_val]]=NA;ta_smt_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_S[[ts_val]]=NA;ta_smt_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_S[[ts_val]]=NA;ta_smt_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_S[[ts_val]]=NA;ta_smt_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_S[[ts_val]]=NA;ta_smt_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_S[[ts_val]]=NA;ta_smt_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_S[[ts_val]]=NA;ta_smt_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_S[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_S[[ts_val]]=NA;ta_smt_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_S[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_smt_genes[[ts_val]]=NA;genes_S[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_S[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_S[[ts_val]]=NA } else {}
				} else 
				if(!is.na(ts_val) && mapapa_name==paste("Agilent-TwoColor",ts_val,sep="_"))
				{
					if(betabeti_name=="...Normalization"){ ta_onl_norm2[[ts_val]]=NA;data.matrix_onlineNorm.m2[[ts_val]]=NA } else 
					if(betabeti_name=="...QC"){ } else 
					if(betabeti_name=="...Filtering"){ ta_onl_filt2[[ts_val]]=NA; data.matrix_onlineNorm.f3[[ts_val]]=NA } else 
					if(betabeti_name=="...Statistical_Significant"){ ta_onl_stat2[[ts_val]]=NA; data.matrix_onlineNorm.s3[[ts_val]]=NA } else 
					if(betabeti_name=="...DGE"){ ta_onl_dge2[[ts_val]]=NA; DE_O_2[[ts_val]]=NA } else 
					if(betabeti_name=="...PCA"){ pca_O[[ts_val]]=NA } else 
					if(betabeti_name=="...Clustering"){ sample.clust_O[[ts_val]]=NA } else 
					if(betabeti_name=="...Classification"){ clas_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP"){ ta_onl_gsea_gobp[[ts_val]]=NA;GOresultBP_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF"){ ta_onl_gsea_gomf[[ts_val]]=NA;GOresultMF_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC"){ ta_onl_gsea_gocc[[ts_val]]=NA;GOresultCC_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG"){ ta_onl_gsea_kegg[[ts_val]]=NA;KEGGresult_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP"){ ta_onl_gsta_gobp[[ts_val]]=NA;GOtable.outBP_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF"){ ta_onl_gsta_gomf[[ts_val]]=NA;GOtable.outMF_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC"){ ta_onl_gsta_gocc[[ts_val]]=NA;GOtable.outCC_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG"){ ta_onl_gsta_kegg[[ts_val]]=NA;KEGGtable.out_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOBP_Graph"){ graph_gsea_goBP_O[[ts_val]]=NA;ta_onl_legend_gsea_gobp[[ts_val]]=NA;legend_gsea_goBP_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOMF_Graph"){ graph_gsea_goMF_O[[ts_val]]=NA;ta_onl_legend_gsea_gomf[[ts_val]]=NA;legend_gsea_goMF_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_GOCC_Graph"){ graph_gsea_goCC_O[[ts_val]]=NA;ta_onl_legend_gsea_gocc[[ts_val]]=NA;legend_gsea_goCC_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSEA_KEGG_Graph"){ graph_gsea_KEGG_O[[ts_val]]=NA;ta_onl_legend_gsea_kegg[[ts_val]]=NA;legend_gsea_KEGG_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOBP_Graph"){ graph_gsta_goBP_O[[ts_val]]=NA;ta_onl_legend_gsta_gobp[[ts_val]]=NA;legend_gsta_goBP_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOMF_Graph"){ graph_gsta_goMF_O[[ts_val]]=NA;ta_onl_legend_gsta_gomf[[ts_val]]=NA;legend_gsta_goMF_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_GOCC_Graph"){ graph_gsta_goCC_O[[ts_val]]=NA;ta_onl_legend_gsta_gocc[[ts_val]]=NA;legend_gsta_goCC_O[[ts_val]]=NA } else
					if(betabeti_name=="...GSTA_KEGG_Graph"){ graph_gsta_KEGG_O[[ts_val]]=NA;ta_onl_legend_gsta_kegg[[ts_val]]=NA;legend_gsta_KEGG_O[[ts_val]]=NA } else
					if(betabeti_name=="...Gene_Symbol"){ ta_onl_genes[[ts_val]]=NA;genes_O[[ts_val]]=NA }
					if(betabeti_name=="...Coexprs_Network"){ myGraph_O[[ts_val]]=NA } else
					if(betabeti_name=="...SSE"){ ssize_O[[ts_val]]=NA } else {}
				}
				tkdestroy(tree_w)
			})
		}

		tree_add<-function(h,...){
			if(dim(tree_mods)[1]>1){
				add_chksb=NULL
				chksb=NULL
				tree_mods2<<-NULL
				chksb_val=NULL
				tree_w<<-tktoplevel()
				tree_frame<-ttkframe(tree_w)
				tkpack(tree_frame,expand=TRUE,fill="both")
				chksbut<-ttkbutton(tree_frame,text="Add")
				tkpack(chksbut)
				tkconfigure(chksbut,command=function(){
					for(i in 1:length(chksb)){
						tree_rein<-tcl(treeview,"insert",tkselect(treeview),"end",values=chksb[i])
						tree_rein_n<-tree_mods2[(which(tree_mods2[,3]==chksb[i])),1]
						tree_mods<<-tree_mods[-(which(tree_mods[,1]==tree_rein_n)),]					
						tree_mods2<<-tree_mods2[-(which(tree_mods2[,3]==chksb[i])),]
					}
					tkdestroy(tree_w)
				})
				chksb_fun<-function(h,...){
					tkconfigure(h,command=function(){
						chksb_name<-tclvalue(tkcget(h,"-text"))
						chksb_status<-tclvalue(tkcget(h,"-variable"))
						if(tclvalue(chksb_status)==1)chksb<<-c(chksb,chksb_name)
					})
				}
				tree_part_name<-as.character(tcl(treeview,"item",tkselect(treeview),"-values"))
				tree_mods2<<-tree_mods[which(tree_mods[,2]==tree_part_name),]
				tree_mods2<<-rbind(c(1,"",""),tree_mods2)
				
				for(i in 2:dim(tree_mods2)[1]){
					chks_txt<-tclVar(tree_mods2[i,3])
					chksb_val[[i]]<-tclVar(FALSE)
					add_chksb[[i]]<-ttkcheckbutton(tree_frame,variable=chksb_val[[i]],textvariable=chks_txt)
					tkpack(add_chksb[[i]],anchor="sw")
				}
				for(i in 2:dim(tree_mods2)[1]){
					chksb_fun(add_chksb[[i]])
				}
			}
		}
		if(length(mapapa_name)!=0){try(tree_detdel(),silent=TRUE)}
		if(length(mapapa_name)==0){try(tree_add(),silent=TRUE)}
	})
}
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk");library(tcltk)} else {library(tcltk)}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
display_mod()
