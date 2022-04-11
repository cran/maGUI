prep_log<-function(h,...){
	r_choice="._."
	treeview<-treeview
	children<-as.character(tcl(treeview,"children",""))
	w_analz<-tktoplevel()
	w_analz_frame<-ttkframe(w_analz,padding=c(3,3,50,20),borderwidth=1,relief="groove")
	tkpack(w_analz_frame,expand=TRUE,fill="both",side="left")
	tree_prnts_smo=NULL
	for(i in 1:length(tree_prnts)){
		tree_slc<-strsplit(tree_prnts[i],"_")
		if(tree_slc[[1]][1]=="Series-Matrix" || tree_slc[[1]][1]=="Online-Data"){
			tree_prnts_smo<-c(tree_prnts_smo,tree_prnts[i])
		}
	}
	var<-tclVar(tree_prnts_smo[1])
	sapply(tree_prnts_smo,function(i){
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
	
		if(ts_name=="Series-Matrix"){
pb_ma<-tkProgressBar(title="Log Transformation...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			data.matrix_series<-l_data.matrix_series[[ts_val]]
			Med<-median(data.matrix_series,na.rm=T)
			if(Med>16)
			{
				data.matrixLog<-log2(data.matrix_series)
			} else {
				data.matrixLog<-data.matrix_series
			}
			l_data.matrixLog[[ts_val]]<<-data.matrixLog
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
			
		}
		if(ts_name=="Online-Data"){
pb_ma<-tkProgressBar(title="Log Transformation...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
			data.matrix_online<-l_data.matrix_online[[ts_val]]
			Med<-median(data.matrix_online,na.rm=T)
			if(Med>16)
			{
				data.matrix_onlineLog<-log2(data.matrix_online)
			} else {
				data.matrix_onlineLog<-data.matrix_online
			}
			l_data.matrix_onlineLog[[ts_val]]<<-data.matrix_onlineLog
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
				
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
prep_log()
