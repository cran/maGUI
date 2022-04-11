qlchk<-function(h,...){
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
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(dat2Affy.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)			
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(dat2Affy.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				tkpack(ma_frame2)			
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Agilent-OneColor"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datAgOne2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgOne2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Agilent-TwoColor"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datAgTwo2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datAgTwo2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Illumina-Beadarray"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(datIllBA2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(datIllBA2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Illumina-Lumi"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(lumi_NQ.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(lumi_NQ.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Nimblegen"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(data.matrix_Nimblegen2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_Nimblegen2.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Series-Matrix"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(data.matrixNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrixNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
		if(ts_name=="Online-Data"){
			if(is.null(ma_img)==TRUE){
				try(tkpack.forget(tableData,scrX,scrY),silent=TRUE)
				try(tkwm.withdraw(w_legend),silent=TRUE)
				ma_img<<-tkrplot(getToolkitWidget(ma_frame2),function(...){boxplot(data.matrix_onlineNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
				ma_img2<<-ma_img
				tkpack(ma_img2)
			} else {
				ma_img2<<-tkrreplot(ma_img,function(...){boxplot(data.matrix_onlineNorm.m2[[ts_val]])},hscale=1.7,vscale=1.07)
			}
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...QC"))
			}
		}
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
qlchk()
