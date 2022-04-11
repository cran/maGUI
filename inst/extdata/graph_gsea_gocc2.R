graph_gsea_gocc<-function(h,...){

	graph_gsea_goCC_params<-function(h,h2,h3,h4,h5,h6,h7,...){
		ts_name=h;ts_val=h2;gocc<-h3;p_v_1<-h4;p_v_2<-h5;children=h6;r_choice=h7;
pb_ma<-tkProgressBar(title="Graph GSEA GO Cellular Component...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
		p1<-rownames(gocc)[gocc$Pvalue>=p_v_1]
		p2<-rownames(gocc)[gocc$Pvalue<=p_v_2]
		p<-p2[match(p1,p2)]
		GO.vec<-p[!is.na(p)==TRUE]
		g<-GOGraph(GO.vec,GOCCPARENTS)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
		g<-removeNode("all",g)
		nodes<-buildNodeList(g)
		focusnode<-sapply(nodes,name) %in% GO.vec
		names(focusnode)<-names(nodes)
		nodefill<-ifelse(focusnode,"yellow","white")
		nAttrs<-list()
		nAttrs$fillcolor<-nodefill
		graph_gsea_goCC<-plot(g,nodeAttrs=nAttrs)
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
		terms<-getGOTerm(nodes(g))
		legend<-data.frame(terms)
		g_colrs<-data.frame(nodefill)
		legend_gsea_goCC<-data.frame(g_colrs,legend)
		colnames(legend_gsea_goCC)<-c("Node_Color","CC")
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
		ta_legend_gsea_gocc<-tclArray()
		for(i in 0:dim(legend_gsea_goCC)[1]){ for(j in 0:dim(legend_gsea_goCC)[2]){ 
			if(i==0){ ta_legend_gsea_gocc[[i,j]]<-colnames(legend_gsea_goCC)[j] } else { 
			if(j==0){ ta_legend_gsea_gocc[[i,j]]<-rownames(legend_gsea_goCC)[i] } else {
			tem<-legend_gsea_goCC[i,j]
			ta_legend_gsea_gocc[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
			} }
		} }
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
		if(ts_name=="Affymetrix"){
			graph_gsea_goCC_Affy[[ts_val]]<<-g
			nodefill_gsea_goCC_Affy[[ts_val]]<<-nodefill
			legend_gsea_goCC_Affy[[ts_val]]<<-legend_gsea_goCC
			ta_affy_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_Affy[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_Affy[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_affy_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Affy[[ts_val]]),ncol(legend_gsea_goCC_Affy[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Agilent-OneColor"){
			graph_gsea_goCC_Ag1[[ts_val]]<<-g
			nodefill_gsea_goCC_Ag1[[ts_val]]<<-nodefill
			legend_gsea_goCC_Ag1[[ts_val]]<<-legend_gsea_goCC
			ta_ag1_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_Ag1[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_Ag1[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_ag1_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Ag1[[ts_val]]),ncol(legend_gsea_goCC_Ag1[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Agilent-TwoColor"){
			graph_gsea_goCC_Ag2[[ts_val]]<<-g
			nodefill_gsea_goCC_Ag2[[ts_val]]<<-nodefill
			legend_gsea_goCC_Ag2[[ts_val]]<<-legend_gsea_goCC
			ta_ag2_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_Ag2[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_Ag2[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_ag2_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Ag2[[ts_val]]),ncol(legend_gsea_goCC_Ag2[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Illumina-Beadarray"){
			graph_gsea_goCC_Il_B[[ts_val]]<<-g
			nodefill_gsea_goCC_Il_B[[ts_val]]<<-nodefill
			legend_gsea_goCC_Il_B[[ts_val]]<<-legend_gsea_goCC
			ta_ilb_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_Il_B[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_Il_B[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_ilb_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Il_B[[ts_val]]),ncol(legend_gsea_goCC_Il_B[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Illumina-Lumi"){
			graph_gsea_goCC_Il_L[[ts_val]]<<-g
			nodefill_gsea_goCC_Il_L[[ts_val]]<<-nodefill
			legend_gsea_goCC_Il_L[[ts_val]]<<-legend_gsea_goCC
			ta_ill_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_Il_L[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_Il_L[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_ill_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_Il_L[[ts_val]]),ncol(legend_gsea_goCC_Il_L[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Nimblegen"){
			graph_gsea_goCC_N[[ts_val]]<<-g
			nodefill_gsea_goCC_N[[ts_val]]<<-nodefill
			legend_gsea_goCC_N[[ts_val]]<<-legend_gsea_goCC
			ta_nbl_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_N[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_N[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_nbl_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_N[[ts_val]]),ncol(legend_gsea_goCC_N[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Series-Matrix"){
			graph_gsea_goCC_S[[ts_val]]<<-g
			nodefill_gsea_goCC_S[[ts_val]]<<-nodefill
			legend_gsea_goCC_S[[ts_val]]<<-legend_gsea_goCC
			ta_smt_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_S[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_S[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_smt_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_S[[ts_val]]),ncol(legend_gsea_goCC_S[[ts_val]]),"GSEA_GOCC_Graph_Legend")
		}
		if(ts_name=="Online-Data"){
			graph_gsea_goCC_O[[ts_val]]<<-g
			nodefill_gsea_goCC_O[[ts_val]]<<-nodefill
			legend_gsea_goCC_O[[ts_val]]<<-legend_gsea_goCC
			ta_onl_legend_gsea_gocc[[ts_val]]<<-ta_legend_gsea_gocc
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){plot(graph_gsea_goCC_O[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			}else{
				ma_img2<<-tkrreplot(ma_img,function(...){plot(graph_gsea_goCC_O[[ts_val]],nodeAttrs=nAttrs)},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...GSEA_GOCC_Graph"))
			}
			legend_tab(ta_onl_legend_gsea_gocc[[ts_val]],nrow(legend_gsea_goCC_O[[ts_val]]),ncol(legend_gsea_goCC_O[[ts_val]]),"GSEA_GOCC_Graph_Legend")
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
		tkwm.title(w_imp,"Select P-values")
		tkwm.resizable(w_imp,0,0)
		frame_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame_imp,text='Lower limit'),side="left",padx=30)
		tkpack(ttklabel(frame_imp,text='Upper limit'),side="left",padx=30)
		tkpack(frame_imp,expand=TRUE,fill="both",side="top")
		frame2_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		p_list<-c(1,0.5,0.1,0.05,0.01,0.001,0.0001,0.00001,0.0000000001,0)
		p_value1<-tclVar(p_list[5])
		p_comb1<-ttkcombobox(frame2_imp,values=p_list,textvariable=p_value1,state="normal",justify="left",width=10)
		tkpack(p_comb1,side="left",padx=25)
		p_value2<-tclVar(p_list[4])
		p_comb2<-ttkcombobox(frame2_imp,values=p_list,textvariable=p_value2,state="normal",justify="left",width=10)
		tkpack(p_comb2,side="left",padx=12)
		tkpack(frame2_imp,expand=TRUE,fill="both",side="top")
		frame3_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		r_but<-ttkbutton(frame3_imp,text="Ok")
		tkpack(r_but,side="right",padx=12)
		tkpack(frame3_imp,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but,command=function(){
			tkwm.withdraw(w_imp)
			p_v_1<-as.numeric(tclvalue(p_value1))
			p_v_2<-as.numeric(tclvalue(p_value2))
			if(ts_name=="Affymetrix"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_Affy[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Agilent-OneColor"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_Ag1[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Agilent-TwoColor"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_Ag2[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Illumina-Beadarray"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_Il_B[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Illumina-Lumi"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_Il_L[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Nimblegen"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_N[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Series-Matrix"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_S[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
			if(ts_name=="Online-Data"){
				graph_gsea_goCC_params(ts_name,ts_val,GOresultCC_O[[ts_val]],p_v_1,p_v_2,children,r_choice)
			}
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
if(!requireNamespace("GO.db",quietly=TRUE)){BiocManager::install("GO.db",ask=FALSE,update=FALSE);library(GO.db)} else {library(GO.db)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("GOstats",quietly=TRUE)){BiocManager::install("GOstats",ask=FALSE,update=FALSE);library(GOstats)} else {library(GOstats)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("graph",quietly=TRUE)){BiocManager::install("graph",ask=FALSE,update=FALSE);library(graph)} else {library(graph)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("Rgraphviz",quietly=TRUE)){install.packages("Rgraphviz",ask=FALSE,update=FALSE);library(Rgraphviz)} else {library(Rgraphviz)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("tcltk",quietly=TRUE)){install.packages("tcltk",ask=FALSE,update=FALSE);library(tcltk)} else {library(tcltk)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
graph_gsea_gocc()
