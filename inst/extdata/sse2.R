sse<-function(h,...){
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
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(dat2Affy.m2[[ts_val]])
			ssize_Affy[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_Affy[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Affy[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Agilent-OneColor"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(datAgOne2.m2[[ts_val]])
			ssize_Ag1[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_Ag1[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Ag1[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Agilent-TwoColor"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(datAgTwo2.m2[[ts_val]])
			ssize_Ag2[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_Ag2[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Ag2[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Illumina-Beadarray"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(datIllBA2.m2[[ts_val]])
			ssize_Il_B[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_Il_B[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Il_B[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Illumina-Lumi"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(lumi_NQ.m2[[ts_val]])
			ssize_Il_L[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_Il_L[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_Il_L[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Nimblegen"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(data.matrix_Nimblegen2.m2[[ts_val]])
			ssize_N[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_N[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_N[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Series-Matrix"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(data.matrixNorm.m2[[ts_val]])
			ssize_S[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_S[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_S[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
		if(ts_name=="Online-Data"){
pb_ma<-tkProgressBar(title="Sample Size Estimation...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
			sds<-rowSds(data.matrix_onlineNorm.m2[[ts_val]])
			ssize_O[[ts_val]]<<-ssize(sd=sds,delta=log2(2),sig.level=0.05,power=0.8)
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){
					ssize.plot(ssize_O[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){
					ssize.plot(ssize_O[[ts_val]],xlim=c(0,20),main=paste("Sample size to detect 2-fold change",sep=""))
				},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...SSE"))
			}
		}
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("ssize",quietly=TRUE)){install.packages("ssize");library(ssize)} else {library(ssize)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("genefilter",quietly=TRUE)){BiocManager::install("genefilter",ask=FALSE,update=FALSE);library(genefilter)} else {library(genefilter)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
sse()
