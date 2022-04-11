display<-function(h,...){
#	treeview<-treeview
	tkbind(treeview,"<Double-Button-1>",function(W,x,y){
		mapapa<-tcl(treeview,"parent",tkselect(treeview))
		mapapa_name<<-as.character(tcl(treeview,"item",mapapa,"-values"))
		betabeti_name<<-as.character(tcl(treeview,"item",tkselect(treeview),"-values"))
		print(mapapa_name)
		print(betabeti_name)
		tabfun<-function(h,tnr,tnc,...){
			try(tkpack.forget(scrX,scrY,tableData,ma_img),silent=TRUE)
			scrX<<-ttkscrollbar(ma_frame2,orient="horizont",command=function(...)tcl(tableData,'xview',...))
			scrY<<-ttkscrollbar(ma_frame2,orient="vertical",command=function(...)tcl(tableData,'yview',...))
			fxscroll<-function(...){tcl(scrX,'set',...)}
			fyscroll<-function(...){tcl(scrY,'set',...)}
			tableData<<-tkwidget(ma_frame2,'table',rows=tnr+1,cols=tnc+1,height=-1,width=-1,ellipsis='...',insertofftime=0,
				flashmode=TRUE,flashtime=1,anchor='w',resizeborders='col',wrap=FALSE,font='{Times}',padx=5,pady=2,ipadx=3,ipady=1,
				rowstretchmode='unset',colstretchmode='fill',multiline=FALSE,cache=TRUE,background="white",foreground="black",selectmode="extended",
				selecttitle=TRUE,relief='groove',borderwidth=c(0,1,0,1),drawmode='compatible',#colwidth=12,
				highlightcolor="gray",highlightbackground="white",highlightthickness=1,xscrollcommand=fxscroll,yscrollcommand=fyscroll,
				rowseparator='\n',colseparator='\t',validate=TRUE,variable=h)
			tcl(tableData,"tag","celltag","ZeroZero","0,0")
			tcl(tableData,"tag","rowtag","rowtitle","0")
			tcl(tableData,"tag","coltag","coltitle","0")
			tcl(tableData,"tag","configure","rowtitle",bg='lightgray',relief='groove',anchor='w')#,borderwidth=c(1,1,1,1))
			tcl(tableData,"tag","configure","coltitle",bg='lightgray',relief='groove',anchor='w')#,borderwidth=c(1,1,1,1))
			tcl(tableData,"tag","configure","active",fg='green',bg='gray90',relief='solid',borderwidth=c(1,1,1,1))
			tcl(tableData,"width","0","12")
			tkpack(scrX,side="bottom",fill="x",expand=FALSE)	# x scrollbar
			tkpack(scrY,side="right",fill="y",expand=FALSE,pady=c(0,18)) #y scrollbar
#			tcl('pack',tableData,expand=TRUE,fill="both",side='left',anchor='e')
			tkpack(tableData,expand=TRUE,fill="both",side='left',anchor='e')
		}
		
		legend_tab<<-function(h1,h2,h3,h4,...){
			w_legend<<-tktoplevel()
			tkwm.title(w_legend,h4)
			legend_frame<-ttkframe(w_legend,padding=c(3,3,20,20),borderwidth=1,relief="groove")
			tkpack(legend_frame,expand=TRUE,fill="both",side="left")
			tabfun2<-function(h,tnr,tnc,...){
				scrX2<-ttkscrollbar(legend_frame,orient="horizont",command=function(...)tcl(tableData2,'xview',...))
				scrY2<-ttkscrollbar(legend_frame,orient="vertical",command=function(...)tcl(tableData2,'yview',...))
				fxscroll<-function(...){tcl(scrX2,'set',...)}
				fyscroll<-function(...){tcl(scrY2,'set',...)}
				tableData2<-tkwidget(legend_frame,'table',rows=tnr+1,cols=tnc+1,height=-1,width=-1,ellipsis='...',insertofftime=0,
					flashmode=TRUE,flashtime=1,anchor='w',resizeborders='col',wrap=FALSE,font='{Times}',padx=5,pady=2,ipadx=3,ipady=1,
						rowstretchmode='unset',colstretchmode='fill',multiline=FALSE,cache=TRUE,background="white",foreground="black",selectmode="extended",
					selecttitle=TRUE,relief='groove',borderwidth=c(0,1,0,1),drawmode='compatible',#colwidth=12,
					highlightcolor="gray",highlightbackground="white",highlightthickness=1,xscrollcommand=fxscroll,yscrollcommand=fyscroll,
					rowseparator='\n',colseparator='\t',validate=TRUE,variable=h)
				tcl(tableData2,"tag","celltag","ZeroZero","0,0")
				tcl(tableData2,"tag","rowtag","rowtitle","0")
				tcl(tableData2,"tag","coltag","coltitle","0")
				tcl(tableData2,"tag","configure","rowtitle",bg='lightgray',relief='groove',anchor='w')#,borderwidth=c(1,1,1,1))
				tcl(tableData2,"tag","configure","coltitle",bg='lightgray',relief='groove',anchor='w')#,borderwidth=c(1,1,1,1))
				tcl(tableData2,"tag","configure","active",fg='green',bg='gray90',relief='solid',borderwidth=c(1,1,1,1))
				tcl(tableData2,"width","0","12")
				tkpack(scrX2,side="bottom",fill="x",expand=FALSE)	# x scrollbar
				tkpack(scrY2,side="right",fill="y",expand=FALSE,pady=c(0,18)) #y scrollbar
	#			tcl('pack',tableData2,expand=TRUE,fill="both",side='left',anchor='e')
				tkpack(tableData2,expand=TRUE,fill="both",side='left',anchor='e')
			}
			tabfun2(h1,h2,h3)
		}


		tree_slc<-strsplit(mapapa_name[1],"_")
		ts_val<-as.numeric(tree_slc[[1]][2])
		if(!is.na(ts_val) && mapapa_name==paste("Affymetrix",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_norm2[[ts_val]],nrow(dat2Affy.m2[[ts_val]]),ncol(dat2Affy.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(dat2Affy.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_filt2[[ts_val]],nrow(dat2Affy.f3[[ts_val]]),ncol(dat2Affy.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_stat2[[ts_val]],nrow(dat2Affy.s3[[ts_val]]),ncol(dat2Affy.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_dge2[[ts_val]],nrow(DE_Affy_2[[ts_val]]),ncol(DE_Affy_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Affy[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Affy[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_Affy[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsea_gobp[[ts_val]],nrow(GOresultBP_Affy[[ts_val]]),ncol(GOresultBP_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsea_gomf[[ts_val]],nrow(GOresultMF_Affy[[ts_val]]),ncol(GOresultMF_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsea_gocc[[ts_val]],nrow(GOresultCC_Affy[[ts_val]]),ncol(GOresultCC_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsea_kegg[[ts_val]],nrow(KEGGresult_Affy[[ts_val]]),ncol(KEGGresult_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsta_gobp[[ts_val]],nrow(GOtable.outBP_Affy[[ts_val]]),ncol(GOtable.outBP_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsta_gomf[[ts_val]],nrow(GOtable.outMF_Affy[[ts_val]]),ncol(GOtable.outMF_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsta_gocc[[ts_val]],nrow(GOtable.outCC_Affy[[ts_val]]),ncol(GOtable.outCC_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_gsta_kegg[[ts_val]],nrow(KEGGtable.out_Affy[[ts_val]]),ncol(KEGGtable.out_Affy[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_Affy[[ts_val]]
					plot(graph_gsea_goBP_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_Affy[[ts_val]]),ncol(legend_gsea_goBP_Affy[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_Affy[[ts_val]]
					plot(graph_gsea_goMF_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_Affy[[ts_val]]),ncol(legend_gsea_goMF_Affy[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_Affy[[ts_val]]
					plot(graph_gsea_goCC_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Affy[[ts_val]]),ncol(legend_gsea_goCC_Affy[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_Affy[[ts_val]],fillcolor=nodefill_gsea_KEGG_Affy[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_Affy[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Affy[[ts_val]]),ncol(legend_gsea_KEGG_Affy[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_Affy[[ts_val]]
					plot(graph_gsta_goBP_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_Affy[[ts_val]]),ncol(legend_gsta_goBP_Affy[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_Affy[[ts_val]]
					plot(graph_gsta_goMF_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_Affy[[ts_val]]),ncol(legend_gsta_goMF_Affy[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_Affy[[ts_val]]
					plot(graph_gsta_goCC_Affy[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_Affy[[ts_val]]),ncol(legend_gsta_goCC_Affy[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_Affy[[ts_val]],fillcolor=nodefill_gsta_KEGG_Affy[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_Affy[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_affy_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_Affy[[ts_val]]),ncol(legend_gsta_KEGG_Affy[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_affy_genes[[ts_val]],nrow(genes_Affy[[ts_val]]),ncol(genes_Affy[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_Affy[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_Affy[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Affy[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Agilent-OneColor",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_norm2[[ts_val]],nrow(datAgOne2.m2[[ts_val]]),ncol(datAgOne2.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgOne2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_filt2[[ts_val]],nrow(datAgOne2.f3[[ts_val]]),ncol(datAgOne2.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_stat2[[ts_val]],nrow(datAgOne2.s3[[ts_val]]),ncol(datAgOne2.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_dge2[[ts_val]],nrow(DE_Ag1_2[[ts_val]]),ncol(DE_Ag1_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Ag1[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Ag1[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_Ag1[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsea_gobp[[ts_val]],nrow(GOresultBP_Ag1[[ts_val]]),ncol(GOresultBP_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsea_gomf[[ts_val]],nrow(GOresultMF_Ag1[[ts_val]]),ncol(GOresultMF_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsea_gocc[[ts_val]],nrow(GOresultCC_Ag1[[ts_val]]),ncol(GOresultCC_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsea_kegg[[ts_val]],nrow(KEGGresult_Ag1[[ts_val]]),ncol(KEGGresult_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsta_gobp[[ts_val]],nrow(GOtable.outBP_Ag1[[ts_val]]),ncol(GOtable.outBP_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsta_gomf[[ts_val]],nrow(GOtable.outMF_Ag1[[ts_val]]),ncol(GOtable.outMF_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsta_gocc[[ts_val]],nrow(GOtable.outCC_Ag1[[ts_val]]),ncol(GOtable.outCC_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_gsta_kegg[[ts_val]],nrow(KEGGtable.out_Ag1[[ts_val]]),ncol(KEGGtable.out_Ag1[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_Ag1[[ts_val]]
					plot(graph_gsea_goBP_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_Ag1[[ts_val]]),ncol(legend_gsea_goBP_Ag1[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_Ag1[[ts_val]]
					plot(graph_gsea_goMF_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_Ag1[[ts_val]]),ncol(legend_gsea_goMF_Ag1[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_Ag1[[ts_val]]
					plot(graph_gsea_goCC_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Ag1[[ts_val]]),ncol(legend_gsea_goCC_Ag1[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_Ag1[[ts_val]],fillcolor=nodefill_gsea_KEGG_Ag1[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_Ag1[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Ag1[[ts_val]]),ncol(legend_gsea_KEGG_Ag1[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_Ag1[[ts_val]]
					plot(graph_gsta_goBP_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_Ag1[[ts_val]]),ncol(legend_gsta_goBP_Ag1[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_Ag1[[ts_val]]
					plot(graph_gsta_goMF_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_Ag1[[ts_val]]),ncol(legend_gsta_goMF_Ag1[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_Ag1[[ts_val]]
					plot(graph_gsta_goCC_Ag1[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_Ag1[[ts_val]]),ncol(legend_gsta_goCC_Ag1[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_Ag1[[ts_val]],fillcolor=nodefill_gsta_KEGG_Ag1[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_Ag1[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag1_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_Ag1[[ts_val]]),ncol(legend_gsta_KEGG_Ag1[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag1_genes[[ts_val]],nrow(genes_Ag1[[ts_val]]),ncol(genes_Ag1[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_Ag1[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_Ag1[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Ag1[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Agilent-TwoColor",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_norm2[[ts_val]],nrow(datAgTwo2.m2[[ts_val]]),ncol(datAgTwo2.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgTwo2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_filt2[[ts_val]],nrow(datAgTwo2.f3[[ts_val]]),ncol(datAgTwo2.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_stat2[[ts_val]],nrow(datAgTwo2.s3[[ts_val]]),ncol(datAgTwo2.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_dge2[[ts_val]],nrow(DE_Ag2_2[[ts_val]]),ncol(DE_Ag2_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Ag2[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Ag2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_Ag2[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsea_gobp[[ts_val]],nrow(GOresultBP_Ag2[[ts_val]]),ncol(GOresultBP_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsea_gomf[[ts_val]],nrow(GOresultMF_Ag2[[ts_val]]),ncol(GOresultMF_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsea_gocc[[ts_val]],nrow(GOresultCC_Ag2[[ts_val]]),ncol(GOresultCC_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsea_kegg[[ts_val]],nrow(KEGGresult_Ag2[[ts_val]]),ncol(KEGGresult_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsta_gobp[[ts_val]],nrow(GOtable.outBP_Ag2[[ts_val]]),ncol(GOtable.outBP_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsta_gomf[[ts_val]],nrow(GOtable.outMF_Ag2[[ts_val]]),ncol(GOtable.outMF_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsta_gocc[[ts_val]],nrow(GOtable.outCC_Ag2[[ts_val]]),ncol(GOtable.outCC_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_gsta_kegg[[ts_val]],nrow(KEGGtable.out_Ag2[[ts_val]]),ncol(KEGGtable.out_Ag2[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_Ag2[[ts_val]]
					plot(graph_gsea_goBP_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_Ag2[[ts_val]]),ncol(legend_gsea_goBP_Ag2[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_Ag2[[ts_val]]
					plot(graph_gsea_goMF_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_Ag2[[ts_val]]),ncol(legend_gsea_goMF_Ag2[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_Ag2[[ts_val]]
					plot(graph_gsea_goCC_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Ag2[[ts_val]]),ncol(legend_gsea_goCC_Ag2[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_Ag2[[ts_val]],fillcolor=nodefill_gsea_KEGG_Ag2[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_Ag2[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Ag2[[ts_val]]),ncol(legend_gsea_KEGG_Ag2[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_Ag2[[ts_val]]
					plot(graph_gsta_goBP_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_Ag2[[ts_val]]),ncol(legend_gsta_goBP_Ag2[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_Ag2[[ts_val]]
					plot(graph_gsta_goMF_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_Ag2[[ts_val]]),ncol(legend_gsta_goMF_Ag2[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_Ag2[[ts_val]]
					plot(graph_gsta_goCC_Ag2[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_Ag2[[ts_val]]),ncol(legend_gsta_goCC_Ag2[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_Ag2[[ts_val]],fillcolor=nodefill_gsta_KEGG_Ag2[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_Ag2[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ag2_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_Ag2[[ts_val]]),ncol(legend_gsta_KEGG_Ag2[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ag2_genes[[ts_val]],nrow(genes_Ag2[[ts_val]]),ncol(genes_Ag2[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_Ag2[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_Ag2[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Ag2[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Illumina-Beadarray",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_norm2[[ts_val]],nrow(datIllBA2.m2[[ts_val]]),ncol(datIllBA2.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datIllBA2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_filt2[[ts_val]],nrow(datIllBA2.f3[[ts_val]]),ncol(datIllBA2.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_stat2[[ts_val]],nrow(datIllBA2.s3[[ts_val]]),ncol(datIllBA2.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_dge2[[ts_val]],nrow(DE_Il_B_2[[ts_val]]),ncol(DE_Il_B_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Il_B[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Il_B[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_Il_B[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsea_gobp[[ts_val]],nrow(GOresultBP_Il_B[[ts_val]]),ncol(GOresultBP_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsea_gomf[[ts_val]],nrow(GOresultMF_Il_B[[ts_val]]),ncol(GOresultMF_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsea_gocc[[ts_val]],nrow(GOresultCC_Il_B[[ts_val]]),ncol(GOresultCC_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsea_kegg[[ts_val]],nrow(KEGGresult_Il_B[[ts_val]]),ncol(KEGGresult_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsta_gobp[[ts_val]],nrow(GOtable.outBP_Il_B[[ts_val]]),ncol(GOtable.outBP_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsta_gomf[[ts_val]],nrow(GOtable.outMF_Il_B[[ts_val]]),ncol(GOtable.outMF_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsta_gocc[[ts_val]],nrow(GOtable.outCC_Il_B[[ts_val]]),ncol(GOtable.outCC_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_gsta_kegg[[ts_val]],nrow(KEGGtable.out_Il_B[[ts_val]]),ncol(KEGGtable.out_Il_B[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_Il_B[[ts_val]]
					plot(graph_gsea_goBP_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_Il_B[[ts_val]]),ncol(legend_gsea_goBP_Il_B[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_Il_B[[ts_val]]
					plot(graph_gsea_goMF_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_Il_B[[ts_val]]),ncol(legend_gsea_goMF_Il_B[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_Il_B[[ts_val]]
					plot(graph_gsea_goCC_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Il_B[[ts_val]]),ncol(legend_gsea_goCC_Il_B[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_Il_B[[ts_val]],fillcolor=nodefill_gsea_KEGG_Il_B[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_Il_B[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Il_B[[ts_val]]),ncol(legend_gsea_KEGG_Il_B[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_Il_B[[ts_val]]
					plot(graph_gsta_goBP_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_Il_B[[ts_val]]),ncol(legend_gsta_goBP_Il_B[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_Il_B[[ts_val]]
					plot(graph_gsta_goMF_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_Il_B[[ts_val]]),ncol(legend_gsta_goMF_Il_B[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_Il_B[[ts_val]]
					plot(graph_gsta_goCC_Il_B[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_Il_B[[ts_val]]),ncol(legend_gsta_goCC_Il_B[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_Il_B[[ts_val]],fillcolor=nodefill_gsta_KEGG_Il_B[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_Il_B[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ilb_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_Il_B[[ts_val]]),ncol(legend_gsta_KEGG_Il_B[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ilb_genes[[ts_val]],nrow(genes_Il_B[[ts_val]]),ncol(genes_Il_B[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_Il_B[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_Il_B[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Il_B[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Illumina-Lumi",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_norm2[[ts_val]],nrow(lumi_NQ.m2[[ts_val]]),ncol(lumi_NQ.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(lumi_NQ.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_filt2[[ts_val]],nrow(lumi_NQ.f3[[ts_val]]),ncol(lumi_NQ.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_stat2[[ts_val]],nrow(lumi_NQ.s3[[ts_val]]),ncol(lumi_NQ.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_dge2[[ts_val]],nrow(DE_Il_L_2[[ts_val]]),ncol(DE_Il_L_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_Il_L[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_Il_L[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_Il_L[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsea_gobp[[ts_val]],nrow(GOresultBP_Il_L[[ts_val]]),ncol(GOresultBP_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsea_gomf[[ts_val]],nrow(GOresultMF_Il_L[[ts_val]]),ncol(GOresultMF_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsea_gocc[[ts_val]],nrow(GOresultCC_Il_L[[ts_val]]),ncol(GOresultCC_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsea_kegg[[ts_val]],nrow(KEGGresult_Il_L[[ts_val]]),ncol(KEGGresult_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsta_gobp[[ts_val]],nrow(GOtable.outBP_Il_L[[ts_val]]),ncol(GOtable.outBP_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsta_gomf[[ts_val]],nrow(GOtable.outMF_Il_L[[ts_val]]),ncol(GOtable.outMF_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsta_gocc[[ts_val]],nrow(GOtable.outCC_Il_L[[ts_val]]),ncol(GOtable.outCC_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_gsta_kegg[[ts_val]],nrow(KEGGtable.out_Il_L[[ts_val]]),ncol(KEGGtable.out_Il_L[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_Il_L[[ts_val]]
					plot(graph_gsea_goBP_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_Il_L[[ts_val]]),ncol(legend_gsea_goBP_Il_L[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_Il_L[[ts_val]]
					plot(graph_gsea_goMF_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_Il_L[[ts_val]]),ncol(legend_gsea_goMF_Il_L[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_Il_L[[ts_val]]
					plot(graph_gsea_goCC_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Il_L[[ts_val]]),ncol(legend_gsea_goCC_Il_L[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_Il_L[[ts_val]],fillcolor=nodefill_gsea_KEGG_Il_L[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_Il_L[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_Il_L[[ts_val]]),ncol(legend_gsea_KEGG_Il_L[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_Il_L[[ts_val]]
					plot(graph_gsta_goBP_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_Il_L[[ts_val]]),ncol(legend_gsta_goBP_Il_L[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_Il_L[[ts_val]]
					plot(graph_gsta_goMF_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_Il_L[[ts_val]]),ncol(legend_gsta_goMF_Il_L[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_Il_L[[ts_val]]
					plot(graph_gsta_goCC_Il_L[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_Il_L[[ts_val]]),ncol(legend_gsta_goCC_Il_L[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_Il_L[[ts_val]],fillcolor=nodefill_gsta_KEGG_Il_L[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_Il_L[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_ill_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_Il_L[[ts_val]]),ncol(legend_gsta_KEGG_Il_L[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_ill_genes[[ts_val]],nrow(genes_Il_L[[ts_val]]),ncol(genes_Il_L[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_Il_L[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_Il_L[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Il_L[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Nimblegen",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_norm2[[ts_val]],nrow(data.matrix_Nimblegen2.m2[[ts_val]]),ncol(data.matrix_Nimblegen2.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_Nimblegen2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_filt2[[ts_val]],nrow(data.matrix_Nimblegen2.f3[[ts_val]]),ncol(data.matrix_Nimblegen2.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_stat2[[ts_val]],nrow(data.matrix_Nimblegen2.s3[[ts_val]]),ncol(data.matrix_Nimblegen2.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_dge2[[ts_val]],nrow(DE_N_2[[ts_val]]),ncol(DE_N_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_N[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_N[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_N[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsea_gobp[[ts_val]],nrow(GOresultBP_N[[ts_val]]),ncol(GOresultBP_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsea_gomf[[ts_val]],nrow(GOresultMF_N[[ts_val]]),ncol(GOresultMF_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsea_gocc[[ts_val]],nrow(GOresultCC_N[[ts_val]]),ncol(GOresultCC_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsea_kegg[[ts_val]],nrow(KEGGresult_N[[ts_val]]),ncol(KEGGresult_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsta_gobp[[ts_val]],nrow(GOtable.outBP_N[[ts_val]]),ncol(GOtable.outBP_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsta_gomf[[ts_val]],nrow(GOtable.outMF_N[[ts_val]]),ncol(GOtable.outMF_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsta_gocc[[ts_val]],nrow(GOtable.outCC_N[[ts_val]]),ncol(GOtable.outCC_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_gsta_kegg[[ts_val]],nrow(KEGGtable.out_N[[ts_val]]),ncol(KEGGtable.out_N[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_N[[ts_val]]
					plot(graph_gsea_goBP_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_N[[ts_val]]),ncol(legend_gsea_goBP_N[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_N[[ts_val]]
					plot(graph_gsea_goMF_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_N[[ts_val]]),ncol(legend_gsea_goMF_N[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_N[[ts_val]]
					plot(graph_gsea_goCC_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_N[[ts_val]]),ncol(legend_gsea_goCC_N[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_N[[ts_val]],fillcolor=nodefill_gsea_KEGG_N[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_N[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_N[[ts_val]]),ncol(legend_gsea_KEGG_N[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_N[[ts_val]]
					plot(graph_gsta_goBP_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_N[[ts_val]]),ncol(legend_gsta_goBP_N[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_N[[ts_val]]
					plot(graph_gsta_goMF_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_N[[ts_val]]),ncol(legend_gsta_goMF_N[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_N[[ts_val]]
					plot(graph_gsta_goCC_N[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_N[[ts_val]]),ncol(legend_gsta_goCC_N[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_N[[ts_val]],fillcolor=nodefill_gsta_KEGG_N[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_N[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_nbl_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_N[[ts_val]]),ncol(legend_gsta_KEGG_N[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_nbl_genes[[ts_val]],nrow(genes_N[[ts_val]]),ncol(genes_N[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_N[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_N[[ts_val]],fontsize=18,fnblcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_N[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Series-Matrix",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_norm2[[ts_val]],nrow(data.matrixNorm.m2[[ts_val]]),ncol(data.matrixNorm.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrixNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_filt2[[ts_val]],nrow(data.matrixNorm.f3[[ts_val]]),ncol(data.matrixNorm.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_stat2[[ts_val]],nrow(data.matrixNorm.s3[[ts_val]]),ncol(data.matrixNorm.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_dge2[[ts_val]],nrow(DE_S_2[[ts_val]]),ncol(DE_S_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_S[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_S[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_S[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsea_gobp[[ts_val]],nrow(GOresultBP_S[[ts_val]]),ncol(GOresultBP_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsea_gomf[[ts_val]],nrow(GOresultMF_S[[ts_val]]),ncol(GOresultMF_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsea_gocc[[ts_val]],nrow(GOresultCC_S[[ts_val]]),ncol(GOresultCC_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsea_kegg[[ts_val]],nrow(KEGGresult_S[[ts_val]]),ncol(KEGGresult_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsta_gobp[[ts_val]],nrow(GOtable.outBP_S[[ts_val]]),ncol(GOtable.outBP_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsta_gomf[[ts_val]],nrow(GOtable.outMF_S[[ts_val]]),ncol(GOtable.outMF_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsta_gocc[[ts_val]],nrow(GOtable.outCC_S[[ts_val]]),ncol(GOtable.outCC_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_gsta_kegg[[ts_val]],nrow(KEGGtable.out_S[[ts_val]]),ncol(KEGGtable.out_S[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_S[[ts_val]]
					plot(graph_gsea_goBP_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bbplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_S[[ts_val]]),ncol(legend_gsea_goBP_S[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_S[[ts_val]]
					plot(graph_gsea_goMF_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_S[[ts_val]]),ncol(legend_gsea_goMF_S[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_S[[ts_val]]
					plot(graph_gsea_goCC_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_S[[ts_val]]),ncol(legend_gsea_goCC_S[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_S[[ts_val]],fillcolor=nodefill_gsea_KEGG_S[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_S[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_S[[ts_val]]),ncol(legend_gsea_KEGG_S[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_S[[ts_val]]
					plot(graph_gsta_goBP_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bbplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_S[[ts_val]]),ncol(legend_gsta_goBP_S[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_S[[ts_val]]
					plot(graph_gsta_goMF_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_S[[ts_val]]),ncol(legend_gsta_goMF_S[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_S[[ts_val]]
					plot(graph_gsta_goCC_S[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_S[[ts_val]]),ncol(legend_gsta_goCC_S[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_S[[ts_val]],fillcolor=nodefill_gsta_KEGG_S[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_S[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_smt_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_S[[ts_val]]),ncol(legend_gsta_KEGG_S[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_smt_genes[[ts_val]],nrow(genes_S[[ts_val]]),ncol(genes_S[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_S[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_S[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_S[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
		if(!is.na(ts_val) && mapapa_name==paste("Online-Data",ts_val,sep="_"))
		{
			if(betabeti_name=="...Normalization")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_norm2[[ts_val]],nrow(data.matrix_onlineNorm.m2[[ts_val]]),ncol(data.matrix_onlineNorm.m2[[ts_val]]))
			} else
			if(betabeti_name=="...QC")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_onlineNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...Filtering")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_filt2[[ts_val]],nrow(data.matrix_onlineNorm.f3[[ts_val]]),ncol(data.matrix_onlineNorm.f3[[ts_val]]))
			} else
			if(betabeti_name=="...Statistical_Significant")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_stat2[[ts_val]],nrow(data.matrix_onlineNorm.s3[[ts_val]]),ncol(data.matrix_onlineNorm.s3[[ts_val]]))
			} else
			if(betabeti_name=="...DGE")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_dge2[[ts_val]],nrow(DE_O_2[[ts_val]]),ncol(DE_O_2[[ts_val]]))
			} else
			if(betabeti_name=="...PCA")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(pca_O[[ts_val]],main="Principle Component Analysis")},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Clustering")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){plot(sample.clust_O[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...Classification")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					heatmap(clas_O[[ts_val]],margins=c(7,7),Rowv=NA,Colv=NA,cexCol=0.8,cexRow=0.8,col=heatcol)
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else
			if(betabeti_name=="...GSEA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsea_gobp[[ts_val]],nrow(GOresultBP_O[[ts_val]]),ncol(GOresultBP_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsea_gomf[[ts_val]],nrow(GOresultMF_O[[ts_val]]),ncol(GOresultMF_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsea_gocc[[ts_val]],nrow(GOresultCC_O[[ts_val]]),ncol(GOresultCC_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsea_kegg[[ts_val]],nrow(KEGGresult_O[[ts_val]]),ncol(KEGGresult_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOBP")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsta_gobp[[ts_val]],nrow(GOtable.outBP_O[[ts_val]]),ncol(GOtable.outBP_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOMF")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsta_gomf[[ts_val]],nrow(GOtable.outMF_O[[ts_val]]),ncol(GOtable.outMF_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_GOCC")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsta_gocc[[ts_val]],nrow(GOtable.outCC_O[[ts_val]]),ncol(GOtable.outCC_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSTA_KEGG")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_gsta_kegg[[ts_val]],nrow(KEGGtable.out_O[[ts_val]]),ncol(KEGGtable.out_O[[ts_val]]))
			} else
			if(betabeti_name=="...GSEA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goBP_O[[ts_val]]
					plot(graph_gsea_goBP_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bbplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsea_gobp[[ts_val]],nrow(legend_gsea_goBP_O[[ts_val]]),ncol(legend_gsea_goBP_O[[ts_val]]),"GSEA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goMF_O[[ts_val]]
					plot(graph_gsea_goMF_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsea_gomf[[ts_val]],nrow(legend_gsea_goMF_O[[ts_val]]),ncol(legend_gsea_goMF_O[[ts_val]]),"GSEA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsea_goCC_O[[ts_val]]
					plot(graph_gsea_goCC_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_O[[ts_val]]),ncol(legend_gsea_goCC_O[[ts_val]]),"GSEA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSEA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsea_KEGG_O[[ts_val]],fillcolor=nodefill_gsea_KEGG_O[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsea_KEGG_O[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsea_kegg[[ts_val]],nrow(legend_gsea_KEGG_O[[ts_val]]),ncol(legend_gsea_KEGG_O[[ts_val]]),"GSEA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOBP_Graph")
			{
				bpplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goBP_O[[ts_val]]
					plot(graph_gsta_goBP_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,bpplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsta_gobp[[ts_val]],nrow(legend_gsta_goBP_O[[ts_val]]),ncol(legend_gsta_goBP_O[[ts_val]]),"GSTA_GOBP_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOMF_Graph")
			{
				mfplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goMF_O[[ts_val]]
					plot(graph_gsta_goMF_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,mfplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsta_gomf[[ts_val]],nrow(legend_gsta_goMF_O[[ts_val]]),ncol(legend_gsta_goMF_O[[ts_val]]),"GSTA_GOMF_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_GOCC_Graph")
			{
				ccplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nAttrs<-list()
					nAttrs$fillcolor<-nodefill_gsta_goCC_O[[ts_val]]
					plot(graph_gsta_goCC_O[[ts_val]],nodeAttrs=nAttrs)
				}
				ma_img2<<-tkrreplot(ma_img,ccplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsta_gocc[[ts_val]],nrow(legend_gsta_goCC_O[[ts_val]]),ncol(legend_gsta_goCC_O[[ts_val]]),"GSTA_GOCC_Graph_Legend")
			} else
			if(betabeti_name=="...GSTA_KEGG_Graph")
			{
				keggplot<-function(){
					try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
					try(tkwm.withdraw(w_legend),silent=TRUE)
					nA<-makeNodeAttrs(graph_gsta_KEGG_O[[ts_val]],fillcolor=nodefill_gsta_KEGG_O[[ts_val]],width=10,height=1.2)
					oldpar<-par(no.readonly=TRUE)
					on.exit(par(oldpar))
					par(mar=c(3,5,0,5),mgp=c(0,0,0))
					layout(mat=matrix(c(rep(1,8),2),ncol=1,byrow=TRUE))
					plot(graph_gsta_KEGG_O[[ts_val]],"dot",nodeAttrs=nA)
					ar<-18
					cols<-colorRampPalette(c("Green", "Red"))(ar)
					image(as.matrix(seq(1,ar)),col=cols,yaxt="n",xaxt="n")
					mtext("down-regulation",side=1,at=0,line=1)
					mtext("up-regulation",side=1,at=1,line=1)
				}
				ma_img2<<-tkrreplot(ma_img,keggplot,hscale=1.7,vscale=1.07)
				tkpack(ma_img)
				legend_tab(ta_onl_legend_gsta_kegg[[ts_val]],nrow(legend_gsta_KEGG_O[[ts_val]]),ncol(legend_gsta_KEGG_O[[ts_val]]),"GSTA_KEGG_Graph_Legend")
			} else
			if(betabeti_name=="...Gene_Symbol")
			{
				try(tkpack.forget(ma_img),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				tabfun(ta_onl_genes[[ts_val]],nrow(genes_O[[ts_val]]),ncol(genes_O[[ts_val]]))
			} else 
			if(betabeti_name=="...Coexprs_Network")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph_O[[ts_val]],nodeAttrs=makeNodeAttrs(myGraph_O[[ts_val]],fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else 
			if(betabeti_name=="...SSE")
			{
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_O[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				tkpack(ma_img)
			} else {}
		}
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
display()
