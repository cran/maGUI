normlz<-function(h,...){
	r_choice="._."
	treeview<-treeview
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
		tree_slc<-strsplit(r_choice,"_")
		ts_name<-tree_slc[[1]][1]
		ts_val<-as.numeric(tree_slc[[1]][2])

		if(ts_name=="Affymetrix"){
			try(({foldAffy<-folder_Affy[[ts_val]];}),silent=TRUE)
			setwd(foldAffy)
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			dat2Affy<-justRMA()
			dat2Affy.m<-exprs(dat2Affy)
			dat2Affy.m<-as.matrix(dat2Affy.m)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
			try({
				if(length(rownames(dat2Affy.m))==0){
					rownames(dat2Affy.m)<-1:length(dat2Affy.m[,1])
				} else {
					rownames(dat2Affy.m)<-as.character(rownames(dat2Affy.m))
				}
			},silent=TRUE)
			ta_affy_norm=NULL
			ta_affy_norm<-tclArray()
			for(i in 0:dim(dat2Affy.m)[1]){ for(j in 0:dim(dat2Affy.m)[2]){ 
				if(i==0){ ta_affy_norm[[i,j]]<-colnames(dat2Affy.m)[j] } else { 
				if(j==0){ ta_affy_norm[[i,j]]<-rownames(dat2Affy.m)[i] } else {
					tem<-dat2Affy.m[i,j]
					ta_affy_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
				} }
			} }
			dat2Affy.m2[[ts_val]]<<-dat2Affy.m
			ta_affy_norm2[[ts_val]]<<-ta_affy_norm
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
			}			
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				
		}
		if(ts_name=="Agilent-OneColor"){
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("none","scale","quantile","cyclicloess")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({foldAg1<-folder_Ag1[[ts_val]];datAgOne<-l_datAgOne[[ts_val]]; }),silent=TRUE)
				setwd(foldAg1)	
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				try(datAgOne2<-backgroundCorrect(datAgOne,"normexp"),silent=TRUE)
				if(length(datAgOne2)==0)datAgOne2<-datAgOne
				datAgOne2<-normalizeBetweenArrays(datAgOne2$R,method=tclvalue(var1))
				datAgOne2<-log(datAgOne2)
				datAgOne2.m<-datAgOne2
				datAgOne2.m<-as.matrix(datAgOne2.m)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(datAgOne2.m))==0){
						rownames(datAgOne2.m)<-1:length(datAgOne2.m[,1])
					} else {
						rownames(datAgOne2.m)<-as.character(rownames(datAgOne2.m))
					}
				},silent=TRUE)
				ta_ag1_norm<-tclArray()
				for(i in 0:dim(datAgOne2.m)[1]){ for(j in 0:dim(datAgOne2.m)[2]){ 
					if(i==0){ ta_ag1_norm[[i,j]]<-colnames(datAgOne2.m)[j] } else { 
					if(j==0){ ta_ag1_norm[[i,j]]<-rownames(datAgOne2.m)[i] } else {
						tem<-datAgOne2.m[i,j]
						ta_ag1_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }		
				datAgOne2.m2[[ts_val]]<<-datAgOne2.m
				ta_ag1_norm2[[ts_val]]<<-ta_ag1_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
		}
		if(ts_name=="Agilent-TwoColor"){
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("none","scale","quantile","Aquantile","Gquantile","Rquantile","Tquantile","cyclicloess")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({foldAg2<-folder_Ag2[[ts_val]];datAgTwo<-l_datAgTwo[[ts_val]]; }),silent=TRUE)
				setwd(foldAg2)
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				try(datAgTwo2<-backgroundCorrect(datAgTwo,"normexp"),silent=TRUE)
				if(length(datAgTwo2)==0)datAgTwo2<-datAgTwo
				datAgTwo2<-normalizeBetweenArrays(datAgTwo2$R,method=tclvalue(var1))
				datAgTwo2<-log(datAgTwo2)
				datAgTwo2.m<-datAgTwo2
				datAgTwo2.m<-as.matrix(datAgTwo2.m)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(datAgTwo2.m))==0){
						rownames(datAgTwo2.m)<-1:length(datAgTwo2.m[,1])
					} else {
						rownames(datAgTwo2.m)<-as.character(rownames(datAgTwo2.m))
					}
				},silent=TRUE)
				ta_ag2_norm=NULL
				ta_ag2_norm<-tclArray()
				for(i in 0:dim(datAgTwo2.m)[1]){ for(j in 0:dim(datAgTwo2.m)[2]){ 
					if(i==0){ ta_ag2_norm[[i,j]]<-colnames(datAgTwo2.m)[j] } else { 
					if(j==0){ ta_ag2_norm[[i,j]]<-rownames(datAgTwo2.m)[i] } else {
					tem<-datAgTwo2.m[i,j]
					ta_ag2_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				datAgTwo2.m2[[ts_val]]<<-datAgTwo2.m				
				ta_ag2_norm2[[ts_val]]<<-ta_ag2_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
		}
		if(ts_name=="Illumina-Beadarray"){
			if (!requireNamespace("beadarray", quietly = TRUE)){
				install.packages("beadarray")
			} else {
				library(beadarray)
			}
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("quantile","qspline","vsn","rankInvariant","median","none","neqc","rsn")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({foldIl_B<-folder_Il_B[[ts_val]];datIllBA<-l_datIllBA[[ts_val]];}),silent=TRUE)
				setwd(foldIl_B)
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				datIllBA2<-normaliseIllumina(datIllBA,method=tclvalue(var1))
				datIllBA2.m<-exprs(datIllBA2)
				datIllBA2.m1<-log2(datIllBA2.m)
				datIllBA2.m1<-as.matrix(datIllBA2.m1)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(datIllBA2.m1))==0){
						rownames(datIllBA2.m1)<-1:length(datIllBA2.m1[,1])
					} else {
						rownames(datIllBA2.m1)<-as.character(rownames(datIllBA2.m1))
					}
				},silent=TRUE)
				ta_ilb_norm<-tclArray()
				for(i in 0:dim(datIllBA2.m1)[1]){ for(j in 0:dim(datIllBA2.m1)[2]){ 
					if(i==0){ ta_ilb_norm[[i,j]]<-colnames(datIllBA2.m1)[j] } else { 
					if(j==0){ ta_ilb_norm[[i,j]]<-rownames(datIllBA2.m1)[i] } else {
					tem<-datIllBA2.m1[i,j]
					ta_ilb_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				datIllBA2.m2[[ts_val]]<<-datIllBA2.m1				
				ta_ilb_norm2[[ts_val]]<<-ta_ilb_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
		}
		if(ts_name=="Illumina-Lumi"){
			if (!requireNamespace("lumi", quietly = TRUE)){install.packages("lumi")} else {	library(lumi) }
			try(({foldIl_L<-folder_Il_L[[ts_val]];lumi_data<-l_lumi_data[[ts_val]];}),silent=TRUE)
			setwd(foldIl_L)	
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
			lumi_NQ<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
			summary(lumi_NQ,'QC')
			lumi_NQ.m<-lumiExpresso(lumi_data,QC.evaluation=TRUE)
			lumi_NQ.m<-as.matrix(lumi_NQ.m)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(lumi_NQ.m))==0){
						rownames(lumi_NQ.m)<-1:length(lumi_NQ.m[,1])
					} else {
						rownames(lumi_NQ.m)<-as.character(rownames(lumi_NQ.m))
					}
				},silent=TRUE)
			ta_ill_norm<-tclArray()
			for(i in 0:dim(lumi_NQ.m)[1]){ for(j in 0:dim(lumi_NQ.m)[2]){ 
				if(i==0){ ta_ill_norm[[i,j]]<-colnames(lumi_NQ.m)[j] } else { 
				if(j==0){ ta_ill_norm[[i,j]]<-rownames(lumi_NQ.m)[i] } else {
				tem<-lumi_NQ.m[i,j]
				ta_ill_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
				} }
			} }
			lumi_NQ.m2[[ts_val]]<<-lumi_NQ.m
			ta_ill_norm2[[ts_val]]<<-ta_ill_norm
			for(i in 1:length(children)){
				x<-as.character(tcl(treeview,"item",children[i],"-values"))
				if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
			}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				
		}
		if(ts_name=="Nimblegen"){
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("none","scale","quantile","cyclicloess")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({foldN<-folder_N[[ts_val]];data.matrix_Nimblegen<-l_data.matrix_Nimblegen[[ts_val]]; }),silent=TRUE)
				setwd(foldN)	
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				data.matrix_Nimblegen2<-normalizeBetweenArrays(data.matrix_Nimblegen,tclvalue(var1))
				data.matrix_Nimblegen2.m<-normalizeBetweenArrays(data.matrix_Nimblegen2,method=tclvalue(var1))
				data.matrix_Nimblegen2.m<-as.matrix(data.matrix_Nimblegen2.m)
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(data.matrix_Nimblegen2.m))==0){
						rownames(data.matrix_Nimblegen2.m)<-1:length(data.matrix_Nimblegen2.m[,1])
					} else {
						rownames(data.matrix_Nimblegen2.m)<-as.character(rownames(data.matrix_Nimblegen2.m))
					}
				},silent=TRUE)
				ta_nbl_norm<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.m)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.m)[2]){ 
					if(i==0){ ta_nbl_norm[[i,j]]<-colnames(data.matrix_Nimblegen2.m)[j] } else { 
					if(j==0){ ta_nbl_norm[[i,j]]<-rownames(data.matrix_Nimblegen2.m)[i] } else {
					tem<-data.matrix_Nimblegen2.m[i,j]
					ta_nbl_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrix_Nimblegen2.m2[[ts_val]]<<-data.matrix_Nimblegen2.m
				ta_nbl_norm2[[ts_val]]<<-ta_nbl_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
		}
		if(ts_name=="Series-Matrix"){
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("none","scale","quantile","cyclicloess")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({foldS<-folder_S[[ts_val]];data.matrixImp<-l_data.matrixImp[[ts_val]];}),silent=TRUE)
				setwd(foldS)	
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				data.matrixNorm<-normalizeBetweenArrays(data.matrixImp,method=tclvalue(var1))
				tmp<-aggregate(data.matrixNorm,list(rownames(data.matrixNorm)),median)
				rownames(tmp)<-tmp[,1]
				data.matrixNorm.m<<-as.matrix(tmp[,-1])
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(data.matrixNorm.m))==0){
						rownames(data.matrixNorm.m)<-1:length(data.matrixNorm.m[,1])
					} else {
						rownames(data.matrixNorm.m)<-as.character(rownames(data.matrixNorm.m))
					}
				},silent=TRUE)
				ta_smt_norm<-tclArray()
				for(i in 0:dim(data.matrixNorm.m)[1]){ for(j in 0:dim(data.matrixNorm.m)[2]){ 
					if(i==0){ ta_smt_norm[[i,j]]<-colnames(data.matrixNorm.m)[j] } else { 
					if(j==0){ ta_smt_norm[[i,j]]<-rownames(data.matrixNorm.m)[i] } else {
						tem<-data.matrixNorm.m[i,j]
						ta_smt_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrixNorm.m2[[ts_val]]<<-data.matrixNorm.m
				
				ta_smt_norm2[[ts_val]]<<-ta_smt_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
		}
		if(ts_name=="Online-Data"){
			w_imp<-tktoplevel()
			tkwm.title(w_imp,"Methods")
			tkwm.resizable(w_imp,0,0)
			frame1_imp<-ttkframe(w_imp,padding=c(3,3,50,20),borderwidth=1,relief="groove")
			tkpack(frame1_imp,expand=TRUE,fill="both",side="left")
			r1_imp<-c("none","scale","quantile","cyclicloess")
			var1<-tclVar(r1_imp[1])
			sapply(r1_imp,function(i){
				r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
				tkpack(r1_imp_vals,side="top",anchor="w")
			})
			tkpack(ttklabel(frame1_imp,text=''))
			tkpack(ttklabel(frame1_imp,text=''))
			r_but<-ttkbutton(frame1_imp,text="Ok")
			tkpack(r_but,pady=3)
			tkconfigure(r_but,command=function(){
				tkwm.withdraw(w_imp)
				try(({data.matrix_onlineImp<-data.matrix_onlineImp[[ts_val]];}),silent=TRUE)
				setwd(folder_O[[ts_val]])
pb_ma<-tkProgressBar(title="Normalization...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.3,title=NULL,label="")
				data.matrix_onlineNorm<-normalizeBetweenArrays(data.matrix_onlineImp,method=tclvalue(var1))
				tmp<-aggregate(data.matrix_onlineNorm,list(rownames(data.matrix_onlineNorm)),median)
				data.matrix_onlineNorm.m<<-as.matrix(tmp[,-1])
setTkProgressBar(pb_ma,0.7,title=NULL,label="")
				try({
					if(length(rownames(data.matrix_onlineNorm.m))==0){
						rownames(data.matrix_onlineNorm.m)<-1:length(data.matrix_onlineNorm.m[,1])
					} else {
						rownames(data.matrix_onlineNorm.m)<-as.character(rownames(data.matrix_onlineNorm.m))
					}
				},silent=TRUE)
				ta_onl_norm<-tclArray()
				for(i in 0:dim(data.matrix_onlineNorm.m)[1]){ for(j in 0:dim(data.matrix_onlineNorm.m)[2]){ 
					if(i==0){ ta_onl_norm[[i,j]]<-colnames(data.matrix_onlineNorm.m)[j] } else { 
					if(j==0){ ta_onl_norm[[i,j]]<-rownames(data.matrix_onlineNorm.m)[i] } else {
						tem<-data.matrix_onlineNorm.m[i,j]
						ta_onl_norm[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrix_onlineNorm.m2[[ts_val]]<<-data.matrix_onlineNorm.m
				ta_onl_norm2[[ts_val]]<<-ta_onl_norm
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Normalization"))
				}
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
			
			})
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
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
if(!requireNamespace("affy",quietly=TRUE)){BiocManager::install("affy",ask=FALSE,update=FALSE);library(affy)} else {library(affy)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
normlz()
