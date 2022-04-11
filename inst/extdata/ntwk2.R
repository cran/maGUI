ntwk<-function(h,...){
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
		if(ts_name=="Affymetrix"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-dat2Affy.m2[[ts_val]][rownames(DE_Affy_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=87,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_Affy[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Agilent-OneColor"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-datAgOne2.m2[[ts_val]][rownames(DE_Ag1_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_Ag1[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Agilent-TwoColor"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-datAgTwo2.m2[[ts_val]][rownames(DE_Ag2_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_Ag2[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Illumina-Beadarray"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-datIllBA2.m2[[ts_val]][rownames(DE_Il_B_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_Il_B[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Illumina-Lumi"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-lumi_NQ.m2[[ts_val]][rownames(DE_Il_L_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_Il_L[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Nimblegen"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-data.matrix_Nimblegen2.m2[[ts_val]][rownames(DE_N_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_N[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Series-Matrix"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-data.matrixNorm.m2[[ts_val]][rownames(DE_S_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_S[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
		if(ts_name=="Online-Data"){
pb_ma<-tkProgressBar(title="Co-expression Network...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.1,title=NULL,label="")
			myData_Sel<-data.matrix_onlineNorm.m2[[ts_val]][rownames(DE_O_2[[ts_val]]),]
			myData_Sel<-t(myData_Sel)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			myMat<-adjacency(myData_Sel,type="signed")
			adjMat<-myMat
			adjMat[abs(adjMat)>0.70]<-1
			adjMat[abs(adjMat)<=0.70]<-0
			diag(adjMat)<-0	
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
			myGraph<-as(adjMat,"graphNEL")
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					plot(myGraph,nodeAttrs=makeNodeAttrs(myGraph,fontsize=18,fillcolor="grey"))
				},hscale=1.7,vscale=1.07)
			}
			myGraph_O[[ts_val]]<<-myGraph
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Coexprs_Network"))
			}
		}
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("WGCNA",quietly=TRUE)){install.packages("WGCNA");library(WGCNA)} else {library(WGCNA)}
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
if(!requireNamespace("Rgraphviz",quietly=TRUE)){install.packages("Rgraphviz",ask=FALSE,update=FALSE);library(Rgraphviz)} else {library(Rgraphviz)}
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("graph",quietly=TRUE)){BiocManager::install("graph",ask=FALSE,update=FALSE);library(graph)} else {library(graph)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
ntwk()
