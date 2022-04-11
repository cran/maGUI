prep_knn<-function(h,...){
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
		tkdestroy(w_analz)
		tree_slc<-strsplit(r_choice[1],"_")
		ts_name<-tree_slc[[1]][1]
		ts_val<-as.numeric(tree_slc[[1]][2])
	
		if(ts_name=="Series-Matrix"){
			na.length=0
			data.matrixLog<-l_data.matrixLog[[ts_val]]
			if(length(data.matrixLog)!=0){
				na.length<-length(which(is.na(data.matrixLog)==T))
			}
			if(na.length>0){
				w_imp<-tktoplevel(height=500,width=1000)
				tkwm.title(w_imp,"KNN-Imputation")
				tkwm.resizable(w_imp,0,0)
				frame1_imp<-ttkframe(w_imp)
				tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame1_imp,text='KNN'))
				r1_imp<-c(2,5,10,50)
				var1<-tclVar(r1_imp[1])
				sapply(r1_imp,function(i){
					r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
					tkpack(r1_imp_vals,side="top",anchor="w")
				})
				frame2_imp<-ttkframe(w_imp)
				tkpack(frame2_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame2_imp,text='Row Max'),padx=25)
				r2_imp<-c(0.3,0.5,0.8)
				var2<-tclVar(r2_imp[2])
				sapply(r2_imp,function(i){
					r2_imp_vals<-ttkradiobutton(frame2_imp,variable=var2,text=i,value=i)
					tkpack(r2_imp_vals,side="top",anchor="w",padx=25)
				})
				frame3_imp<-ttkframe(w_imp)
				tkpack(frame3_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame3_imp,text='Col Max'))
				r3_imp<-c(0.3,0.5,0.8)
				var3<-tclVar(r3_imp[3])
				sapply(r3_imp,function(i){
					r3_imp_vals<-ttkradiobutton(frame3_imp,variable=var3,text=i,value=i)
					tkpack(r3_imp_vals,side="top",anchor="w")
				})
				tkpack(ttklabel(frame2_imp,text=''))
				tkpack(ttklabel(frame2_imp,text=''))
				r_but<-ttkbutton(frame2_imp,text="Ok")
				tkpack(r_but)
				tkconfigure(r_but,command=function(){
					tkdestroy(w_imp)
pb_ma<-tkProgressBar(title="KNN Imputation...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
					l_data.matrixImp[[ts_val]]<<-impute.knn(data.matrixLog,k=as.numeric(tclvalue(var1)),rowmax=as.numeric(tclvalue(var2)),colmax=as.numeric(tclvalue(var3)))$data
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
				})
			} else {
				l_data.matrixImp[[ts_val]]<<-data.matrixLog
			}			
		}
		if(ts_name=="Online-Data"){
			na.length=0
			data.matrix_onlineLog<-l_data.matrix_onlineLog[[ts_val]]
			if(length(data.matrix_onlineLog)!=0){
				na.length<-length(which(is.na(data.matrix_onlineLog)==T))
			}
			if(na.length>0){
				w_imp<-tktoplevel()
				tkwm.title(w_imp,"KNN-Imputation")
				tkwm.resizable(w_imp,0,0)
				frame1_imp<-ttkframe(w_imp)
				tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame1_imp,text='KNN'))
				r1_imp<-c(2,5,10,50)
				var1<-tclVar(r1_imp[1])
				sapply(r1_imp,function(i){
					r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
					tkpack(r1_imp_vals,side="top",anchor="w")
				})
				frame2_imp<-ttkframe(w_imp)
				tkpack(frame2_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame2_imp,text='Row Max'),padx=25)
				r2_imp<-c(0.3,0.5,0.8)
				var2<-tclVar(r2_imp[2])
				sapply(r2_imp,function(i){
					r2_imp_vals<-ttkradiobutton(frame2_imp,variable=var2,text=i,value=i)
					tkpack(r2_imp_vals,side="top",anchor="w",padx=25)
				})
				frame3_imp<-ttkframe(w_imp)
				tkpack(frame3_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame3_imp,text='Col Max'))
				r3_imp<-c(0.3,0.5,0.8)
				var3<-tclVar(r3_imp[3])
				sapply(r3_imp,function(i){
					r3_imp_vals<-ttkradiobutton(frame3_imp,variable=var3,text=i,value=i)
					tkpack(r3_imp_vals,side="top",anchor="w")
				})
#				tkpack(frame3_imp,expand=TRUE,fill="both",side="left")
				tkpack(ttklabel(frame2_imp,text=''))
				tkpack(ttklabel(frame2_imp,text=''))
				r_but<-ttkbutton(frame2_imp,text="Ok")
				tkpack(r_but)
				tkconfigure(r_but,command=function(){
					tkwm.withdraw(w_imp)
pb_ma<-tkProgressBar(title="KNN Imputation...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
					l_data.matrix_onlineImp[[ts_val]]<<-impute.knn(data.matrix_onlineLog,k=as.numeric(tclvalue(var1)),rowmax=as.numeric(tclvalue(var2)),colmax=as.numeric(tclvalue(var3)))$data
setTkProgressBar(pb_ma,1,title="Completed...",label="")				
close(pb_ma)
				})
			} else {
				l_data.matrix_onlineImp[[ts_val]]<<-data.matrix_onlineLog
			}
		}
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.4,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("BiocManager",quietly=TRUE)){install.packages("BiocManager");library(BiocManager)} else {library(BiocManager)}
setTkProgressBar(pb_ma,0.8,title=NULL,label="")
if(!requireNamespace("impute",quietly=TRUE)){BiocManager::install("impute",ask=FALSE,update=FALSE);library(impute)} else {library(impute)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
prep_knn()
