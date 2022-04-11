ma_export<<-function(h,...){
	ma_tab<<-function(h1,...){
		w_imp1<-tktoplevel()
		tkwm.title(w_imp1,'Export table')
		tkwm.resizable(w_imp1,0,0)
		frame1_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
		tkpack(ttklabel(frame1_imp1,text='Give a name:'),side='left',padx=12)
		exp_file<-tclVar('')
		c1<-ttkentry(frame1_imp1,textvariable=exp_file,width=25)
		tkpack(c1,side='left',anchor='e')
		tkpack(frame1_imp1,expand=TRUE,fill='both',side='top')
		frame2_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
		r_but1<-ttkbutton(frame2_imp1,text='Ok')
		tkpack(r_but1,side='top',anchor='e',pady=5)
		tkpack(frame2_imp1,expand=TRUE,fill='both',side='top')
		tkconfigure(r_but1,command=function(){
			tkwm.withdraw(w_imp1)
			out_name<-paste(tclvalue(exp_file),'txt',sep='.')
			folder<-tclvalue(tkchooseDirectory(title='Destination Folder'))
			setwd(folder)
			write.table(h1,file=out_name,sep='\t',quote=FALSE)
		})
	}
	ma_fig<<-function(h2,h3,h4,...){
		w_imp1<-tktoplevel()
		tkwm.title(w_imp1,'Export figure')
		tkwm.resizable(w_imp1,0,0)
		frame1_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
		tkpack(ttklabel(frame1_imp1,text='Give a name:'),side='left',padx=12)
		exp_file<-tclVar('')
		c1<-ttkentry(frame1_imp1,textvariable=exp_file,width=25)
		tkpack(c1,side='left',anchor='e')
		tkpack(frame1_imp1,expand=TRUE,fill='both',side='top')
		frame2_imp1<-ttkframe(w_imp1,padding=c(3,3,12,12))
		r_but1<-ttkbutton(frame2_imp1,text='Ok')
		tkpack(r_but1,side='top',anchor='e',pady=5)
		tkpack(frame2_imp1,expand=TRUE,fill='both',side='top')
		tkconfigure(r_but1,command=function(){
			tkwm.withdraw(w_imp1)
			out_name<-paste(tclvalue(exp_file),'png',sep='.')
			folder<-tclvalue(tkchooseDirectory(title='Destination Folder'))
			setwd(folder)
			if(h3=='QC'){
				png(out_name,width=750,height=750)
				boxplot(h2)
				dev.off()
								
			} else 
			if(h3=='PCA'){
				png(out_name,width=750,height=750)
				plot(h2)
				dev.off()
				
			} else
			if(h3=='Clustered_Data'){
				png(out_name,width=750,height=750)
				plot(h2)
				dev.off()
				
			} else
			if(h3=='Classification_Data'){
				png(out_name,width=750,height=750)
				heatmap(h2,Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				dev.off()
				
			} else
			if(h3=='Graph_GSEA_Data_GO_BP'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSEA_Data_GO_MF'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSEA_Data_GO_CC'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSEA_Data_KEGG_Pathways'){
				png(out_name,width=750,height=750)
				nA<-makeNodeAttrs(h2,fillcolor=h4,width=10,height=1.2)
				oldpar<-par(no.readonly=TRUE)
				on.exit(par(oldpar))
				par(mar=c(3,5,0,5),mgp=c(0,0,0))
				layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
				plot(h2,'dot',nodeAttrs=nA)
				ar<-18
				cols<-colorRampPalette(c('Green', 'Red'))(ar)
				image(as.matrix(seq(1,ar)),col=cols,yaxt='n',xaxt='n')
				mtext('down-regulation',side=1,at=0,line=1)
				mtext('up-regulation',side=1,at=1,line=1)
				dev.off()
				
			} else 
			if(h3=='Graph_GSTA_Data_GO_BP'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSTA_Data_GO_MF'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSTA_Data_GO_CC'){
				png(out_name,width=750,height=750)
				nAttrs<-list()
				nAttrs$fillcolor<-h4
				plot(h2,nodeAttrs=nAttrs)
				dev.off()
				
			} else 
			if(h3=='Graph_GSTA_Data_KEGG_Pathways'){
				png(out_name,width=750,height=750)
				nA<-makeNodeAttrs(h2,fillcolor=h4,width=10,height=1.2)
				oldpar<-par(no.readonly=TRUE)
				on.exit(par(oldpar))
				par(mar=c(3,5,0,5),mgp=c(0,0,0))
				layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
				plot(h2,'dot',nodeAttrs=nA)
				ar<-18
				cols<-colorRampPalette(c('Green', 'Red'))(ar)
				image(as.matrix(seq(1,ar)),col=cols,yaxt='n',xaxt='n')
				mtext('down-regulation',side=1,at=0,line=1)
				mtext('up-regulation',side=1,at=1,line=1)
				dev.off()
				
			} else 
			if(h3=='Co-expression_Network_Data'){
				png(out_name,width=750,height=750)
				plot(h2,nodeAttrs=makeNodeAttrs(h2,fontsize=18,fillcolor='grey'))
				dev.off()
				
			} else 
			if(h3=='SSE_Data'){
				png(out_name,width=750,height=750)
				ssize.plot(h2,xlim=c(0,20),main=paste('Sample size to detect 2-fold change',sep=''))
				dev.off()
				
			} else {}

		})
	}
	export_data<-h
	r_choice='._.'
	children<-as.character(tcl(treeview,'children',''))
	w_analz<-tktoplevel()
	w_analz_frame<-ttkframe(w_analz,padding=c(3,3,50,20),borderwidth=1,relief='groove')
	tkpack(w_analz_frame,expand=TRUE,fill='both',side='left')
	var<-tclVar(tree_prnts[1])
	sapply(tree_prnts,function(i){
		r_but_vals<-ttkradiobutton(w_analz_frame,variable=var,text=i,value=i)
		tkpack(r_but_vals,side='top',anchor='w')
	})
	tkpack(ttklabel(w_analz_frame,text=''))
	tkpack(ttklabel(w_analz_frame,text=''))
	r_but<-ttkbutton(w_analz_frame,text='Ok')
	tkpack(r_but,pady=2)
	tkconfigure(r_but,command=function(){
		r_choice<-tclvalue(var)
		tkwm.withdraw(w_analz)
		tree_slc<-strsplit(r_choice[1],'_')
		ts_name<-tree_slc[[1]][1]
		ts_val<-as.numeric(tree_slc[[1]][2])
		if(ts_name=='Affymetrix'){
			if(export_data=='Normalized_Data'){
				ma_tab(dat2Affy.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(dat2Affy.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(dat2Affy.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(dat2Affy.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_Affy_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_Affy[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_Affy[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_Affy[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_Affy[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_Affy[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_Affy[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_Affy[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_Affy[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_Affy[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_Affy[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_Affy[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_Affy[[ts_val]])
				ma_tab(legend_gsea_goBP_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_Affy[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_Affy[[ts_val]])
				ma_tab(legend_gsea_goMF_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_Affy[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_Affy[[ts_val]])
				ma_tab(legend_gsea_goCC_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_Affy[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_Affy[[ts_val]])
				ma_tab(legend_gsea_KEGG_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_Affy[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_Affy[[ts_val]])
				ma_tab(legend_gsta_goBP_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_Affy[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_Affy[[ts_val]])
				ma_tab(legend_gsta_goMF_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_Affy[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_Affy[[ts_val]])
				ma_tab(legend_gsta_goCC_Affy[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_Affy[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_Affy[[ts_val]])
				ma_tab(legend_gsta_KEGG_Affy[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_Affy[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_Affy[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_Affy[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Agilent-OneColor'){
			if(export_data=='Normalized_Data'){
				ma_tab(datAgOne2.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(datAgOne2.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(datAgOne2.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(datAgOne2.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_Ag1_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_Ag1[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_Ag1[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_Ag1[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_Ag1[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_Ag1[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_Ag1[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_Ag1[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_Ag1[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_Ag1[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_Ag1[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_Ag1[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_Ag1[[ts_val]])
				ma_tab(legend_gsea_goBP_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_Ag1[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_Ag1[[ts_val]])
				ma_tab(legend_gsea_goMF_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_Ag1[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_Ag1[[ts_val]])
				ma_tab(legend_gsea_goCC_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_Ag1[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_Ag1[[ts_val]])
				ma_tab(legend_gsea_KEGG_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_Ag1[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_Ag1[[ts_val]])
				ma_tab(legend_gsta_goBP_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_Ag1[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_Ag1[[ts_val]])
				ma_tab(legend_gsta_goMF_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_Ag1[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_Ag1[[ts_val]])
				ma_tab(legend_gsta_goCC_Ag1[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_Ag1[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_Ag1[[ts_val]])
				ma_tab(legend_gsta_KEGG_Ag1[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_Ag1[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_Ag1[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_Ag1[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Agilent-TwoColor'){
			if(export_data=='Normalized_Data'){
				ma_tab(datAgTwo2.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(datAgTwo2.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(datAgTwo2.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(datAgTwo2.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_Ag2_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_Ag2[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_Ag2[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_Ag2[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_Ag2[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_Ag2[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_Ag2[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_Ag2[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_Ag2[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_Ag2[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_Ag2[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_Ag2[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_Ag2[[ts_val]])
				ma_tab(legend_gsea_goBP_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_Ag2[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_Ag2[[ts_val]])
				ma_tab(legend_gsea_goMF_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_Ag2[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_Ag2[[ts_val]])
				ma_tab(legend_gsea_goCC_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_Ag2[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_Ag2[[ts_val]])
				ma_tab(legend_gsea_KEGG_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_Ag2[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_Ag2[[ts_val]])
				ma_tab(legend_gsta_goBP_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_Ag2[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_Ag2[[ts_val]])
				ma_tab(legend_gsta_goMF_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_Ag2[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_Ag2[[ts_val]])
				ma_tab(legend_gsta_goCC_Ag2[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_Ag2[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_Ag2[[ts_val]])
				ma_tab(legend_gsta_KEGG_Ag2[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_Ag2[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_Ag2[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_Ag2[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Illumina-Beadarray'){
			if(export_data=='Normalized_Data'){
				ma_tab(datIllBA2.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(datIllBA2.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(datIllBA2.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(datIllBA2.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_Il_B_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_Il_B[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_Il_B[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_Il_B[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_Il_B[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_Il_B[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_Il_B[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_Il_B[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_Il_B[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_Il_B[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_Il_B[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_Il_B[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_Il_B[[ts_val]])
				ma_tab(legend_gsea_goBP_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_Il_B[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_Il_B[[ts_val]])
				ma_tab(legend_gsea_goMF_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_Il_B[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_Il_B[[ts_val]])
				ma_tab(legend_gsea_goCC_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_Il_B[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_Il_B[[ts_val]])
				ma_tab(legend_gsea_KEGG_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_Il_B[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_Il_B[[ts_val]])
				ma_tab(legend_gsta_goBP_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_Il_B[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_Il_B[[ts_val]])
				ma_tab(legend_gsta_goMF_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_Il_B[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_Il_B[[ts_val]])
				ma_tab(legend_gsta_goCC_Il_B[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_Il_B[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_Il_B[[ts_val]])
				ma_tab(legend_gsta_KEGG_Il_B[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_Il_B[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_Il_B[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_Il_B[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Illumina-Lumi'){
			if(export_data=='Normalized_Data'){
				ma_tab(lumi_NQ.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(lumi_NQ.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(lumi_NQ.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(lumi_NQ.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_Il_L_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_Il_L[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_Il_L[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_Il_L[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_Il_L[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_Il_L[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_Il_L[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_Il_L[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_Il_L[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_Il_L[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_Il_L[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_Il_L[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_Il_L[[ts_val]])
				ma_tab(legend_gsea_goBP_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_Il_L[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_Il_L[[ts_val]])
				ma_tab(legend_gsea_goMF_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_Il_L[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_Il_L[[ts_val]])
				ma_tab(legend_gsea_goCC_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_Il_L[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_Il_L[[ts_val]])
				ma_tab(legend_gsea_KEGG_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_Il_L[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_Il_L[[ts_val]])
				ma_tab(legend_gsta_goBP_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_Il_L[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_Il_L[[ts_val]])
				ma_tab(legend_gsta_goMF_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_Il_L[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_Il_L[[ts_val]])
				ma_tab(legend_gsta_goCC_Il_L[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_Il_L[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_Il_L[[ts_val]])
				ma_tab(legend_gsta_KEGG_Il_L[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_Il_L[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_Il_L[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_Il_L[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Nimblegen'){
			if(export_data=='Normalized_Data'){
				ma_tab(data.matrix_Nimblegen2.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(data.matrix_Nimblegen2.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(data.matrix_Nimblegen2.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(data.matrix_Nimblegen2.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_N_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_N[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_N[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_N[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_N[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_N[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_N[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_N[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_N[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_N[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_N[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_N[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_N[[ts_val]])
				ma_tab(legend_gsea_goBP_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_N[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_N[[ts_val]])
				ma_tab(legend_gsea_goMF_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_N[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_N[[ts_val]])
				ma_tab(legend_gsea_goCC_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_N[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_N[[ts_val]])
				ma_tab(legend_gsea_KEGG_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_N[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_N[[ts_val]])
				ma_tab(legend_gsta_goBP_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_N[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_N[[ts_val]])
				ma_tab(legend_gsta_goMF_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_N[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_N[[ts_val]])
				ma_tab(legend_gsta_goCC_N[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_N[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_N[[ts_val]])
				ma_tab(legend_gsta_KEGG_N[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_N[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_N[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_N[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Series-Matrix'){
			if(export_data=='Normalized_Data'){
				ma_tab(data.matrixNorm.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(data.matrixNorm.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(data.matrixNorm.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(data.matrixNorm.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_S_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_S[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_S[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_S[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_S[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_S[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_S[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_S[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_S[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_S[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_S[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_S[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_S[[ts_val]])
				ma_tab(legend_gsea_goBP_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_S[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_S[[ts_val]])
				ma_tab(legend_gsea_goMF_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_S[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_S[[ts_val]])
				ma_tab(legend_gsea_goCC_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_S[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_S[[ts_val]])
				ma_tab(legend_gsea_KEGG_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_S[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_S[[ts_val]])
				ma_tab(legend_gsta_goBP_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_S[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_S[[ts_val]])
				ma_tab(legend_gsta_goMF_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_S[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_S[[ts_val]])
				ma_tab(legend_gsta_goCC_S[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_S[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_S[[ts_val]])
				ma_tab(legend_gsta_KEGG_S[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_S[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_S[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_S[[ts_val]],'SSE_Data','')
			} else {}
		}		
		if(ts_name=='Online-Data'){
			if(export_data=='Normalized_Data'){
				ma_tab(data.matrix_onlineNorm.m2[[ts_val]])
			} else 
			if(export_data=='QC_Plot'){
				ma_fig(data.matrix_onlineNorm.m2[[ts_val]],'QC','')
			} else 
			if(export_data=='Filtered_Data'){
				ma_tab(data.matrix_onlineNorm.f3[[ts_val]])
			} else 
			if(export_data=='Stat_Sign_Data'){
				ma_tab(data.matrix_onlineNorm.s3[[ts_val]])
			} else 
			if(export_data=='DGE_Data'){
				ma_tab(DE_O_2[[ts_val]])
			} else 
			if(export_data=='PCA_Data'){
				ma_fig(pca_O[[ts_val]],'PCA','')
			} else 
			if(export_data=='Clustered_Data'){
				ma_fig(sample.clust_O[[ts_val]],'Clustered_Data','')
			} else 
			if(export_data=='Classification_Data'){
				ma_fig(clas_O[[ts_val]],'Classification_Data','')
			} else 
			if(export_data=='GSEA_Data_GO_BP'){
				ma_tab(GOresultBP_O[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_MF'){
				ma_tab(GOresultMF_O[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_GO_CC'){
				ma_tab(GOresultCC_O[[ts_val]])
			} else 
			if(export_data=='GSEA_Data_KEGG_Pathways'){
				ma_tab(KEGGresult_O[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_BP'){
				ma_tab(GOtable.outBP_O[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_MF'){
				ma_tab(GOtable.outMF_O[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_GO_CC'){
				ma_tab(GOtable.outCC_O[[ts_val]])
			} else 
			if(export_data=='GSTA_Data_KEGG_Pathways'){
				ma_tab(KEGGtable.out_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_BP'){
				ma_fig(graph_gsea_goBP_O[[ts_val]],'Graph_GSEA_Data_GO_BP',nodefill_gsea_goBP_O[[ts_val]])
				ma_tab(legend_gsea_goBP_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_MF'){
				ma_fig(graph_gsea_goMF_O[[ts_val]],'Graph_GSEA_Data_GO_MF',nodefill_gsea_goMF_O[[ts_val]])
				ma_tab(legend_gsea_goMF_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_GO_CC'){
				ma_fig(graph_gsea_goCC_O[[ts_val]],'Graph_GSEA_Data_GO_CC',nodefill_gsea_goCC_O[[ts_val]])
				ma_tab(legend_gsea_goCC_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSEA_Data_KEGG_Pathways'){
				ma_fig(graph_gsea_KEGG_O[[ts_val]],'Graph_GSEA_Data_KEGG_Pathways',nodefill_gsea_KEGG_O[[ts_val]])
				ma_tab(legend_gsea_KEGG_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_BP'){
				ma_fig(graph_gsta_goBP_O[[ts_val]],'Graph_GSTA_Data_GO_BP',nodefill_gsta_goBP_O[[ts_val]])
				ma_tab(legend_gsta_goBP_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_MF'){
				ma_fig(graph_gsta_goMF_O[[ts_val]],'Graph_GSTA_Data_GO_MF',nodefill_gsta_goMF_O[[ts_val]])
				ma_tab(legend_gsta_goMF_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_GO_CC'){
				ma_fig(graph_gsta_goCC_O[[ts_val]],'Graph_GSTA_Data_GO_CC',nodefill_gsta_goCC_O[[ts_val]])
				ma_tab(legend_gsta_goCC_O[[ts_val]])
			} else 
			if(export_data=='Graph_GSTA_Data_KEGG_Pathways'){
				ma_fig(graph_gsta_KEGG_O[[ts_val]],'Graph_GSTA_Data_KEGG_Pathways',nodefill_gsta_KEGG_O[[ts_val]])
				ma_tab(legend_gsta_KEGG_O[[ts_val]])
			} else 
			if(export_data=='Gene_Symbol_Data'){
				ma_tab(genes_O[[ts_val]])
			} else 
			if(export_data=='Co-expression_Network_Data'){
				ma_fig(myGraph_O[[ts_val]],'Co-expression_Network_Data','')
			} else 
			if(export_data=='SSE_Data'){
				ma_fig(ssize_O[[ts_val]],'SSE_Data','')
			} else {}
		}		
	})
}
pb_ma<-tkProgressBar(title='Loading packages...',label='',min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label='')
if(!requireNamespace('gWidgets2tcltk',quietly=TRUE)){install.packages('gWidgets2tcltk');library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.3,title=NULL,label='')
if(!requireNamespace('graph',quietly=TRUE)){install.packages('graph');library(graph)} else {library(graph)}
setTkProgressBar(pb_ma,0.5,title=NULL,label='')
if(!requireNamespace('Rgraphviz',quietly=TRUE)){install.packages('Rgraphviz');library(Rgraphviz)} else {library(Rgraphviz)}
setTkProgressBar(pb_ma,0.8,title=NULL,label='')
if(!requireNamespace('tkrplot',quietly=TRUE)){install.packages('tkrplot');library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,1,title='Done...',label='')				
close(pb_ma)
ma_export()
