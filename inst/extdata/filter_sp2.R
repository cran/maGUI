filter_sp<-function(h,...){
	sp_filt<-function(h,h2,h3,h4,...){
		ts_name<-h;ts_val<-h2;children<-h3;r_choice=h4
		w_imp2<-tktoplevel()
		tkwm.title(w_imp2,"Groups")
		tkwm.resizable(w_imp2,0,0)
		frame1_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12))
		tkpack(ttklabel(frame1_imp2,text="Controls"),side="left",padx=60)
		tkpack(ttklabel(frame1_imp2,text="Tests"),side="left",padx=100)
		tkpack(frame1_imp2,expand=TRUE,fill="both",side="top")
		frame2_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12),height=40,width=50,borderwidth=1,relief="groove")
		tkpack(frame2_imp2,expand=TRUE,fill="both",side="left")
		frame3_imp2<-ttkframe(w_imp2,padding=c(3,3,12,12),height=40,width=50,borderwidth=1,relief="groove")
		tkpack(frame3_imp2,expand=TRUE,fill="both",side="right")		
		controls<-NULL;tests<-NULL;add_i=NULL;
		frame2_1_imp2<-ttkframe(frame2_imp2,padding=c(3,3,12,12),height=30,width=50,borderwidth=1,relief="groove")
		tkpack(frame2_1_imp2,expand=TRUE,fill="both",side="top")
		c1<-ttkentry(frame2_1_imp2,textvariable=tclVar(""),width=25)
		tkpack(c1,side="left",anchor="e")
		controls<-c(controls,c1)
		t1<-ttkentry(frame2_1_imp2,textvariable=tclVar(""),width=25)
		tkpack(t1,side="right",anchor="e")
		tests<-c(tests,t1)
		add_i=1
		r_but_add<-ttkbutton(frame3_imp2,text="Add")
		tkpack(r_but_add,side="top",anchor="e",pady=5)
		r_but_ok<-ttkbutton(frame3_imp2,text="Ok")
		tkpack(r_but_ok,side="top",anchor="e",pady=5)
		tkpack(frame3_imp2,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but_add,command=function(){
			frame2_2_imp2<-ttkframe(frame2_imp2,padding=c(3,3,12,12),height=30,width=50,borderwidth=1,relief="groove")
			tkpack(frame2_2_imp2,expand=TRUE,fill="both",side="top")
			c2<-ttkentry(frame2_2_imp2,textvariable=tclVar(""),width=25)
			tkpack(c2,side="left",anchor="e")
			controls<<-c(controls,c2)
			t2<-ttkentry(frame2_2_imp2,textvariable=tclVar(""),width=25)
			tkpack(t2,side="right",anchor="e")
			tests<<-c(tests,t2)
			add_i<<-add_i+2
		})
		tkconfigure(r_but_ok,command=function(){
			i=1
			groups<<-c()
			samp<<-c()
			while(i<=add_i){
				xf<-strsplit(tclvalue(tkget(as.tclObj(controls[i]$ID))),",");
				yf<-strsplit(tclvalue(tkget(as.tclObj(tests[i]$ID))),",");
				if(i>2){
					groups<<-c(groups,rep(paste("T",i,sep=""),length(yf[[1]])))
					samp<<-c(samp,yf[[1]])
				} else {
					groups<<-c(groups,rep(paste("C",i,sep=""),length(xf[[1]])),rep(paste("T",i,sep=""),length(yf[[1]])))
					samp<<-c(samp,xf[[1]],yf[[1]])
				}
				i=i+2
			}	
			tkwm.withdraw(w_imp2)
			groups<-as.factor(groups);
			design<-model.matrix(~groups);
			if(ts_name=="Affymetrix"){
				dat2Affy.m<-dat2Affy.m2[[ts_val]];fit<-lmFit(dat2Affy.m[,samp],design);print(fit)
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(dat2Affy.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_Affy[[ts_val]]<<-fit2
				dat2Affy.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_Affy[[ts_val]]<<-samp
				ps1<-as.data.frame(dat2Affy.f)
				rownames(ps1)<-rownames(fit2)
				dat2Affy.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_affy_filt<-tclArray()
				for(i in 0:dim(dat2Affy.f2)[1]){ for(j in 0:dim(dat2Affy.f2)[2]){ 
					if(i==0){ ta_affy_filt[[i,j]]<-colnames(dat2Affy.f2)[j] } else { 
					if(j==0){ ta_affy_filt[[i,j]]<-rownames(dat2Affy.f2)[i] } else {
					tem<-dat2Affy.f2[i,j]
					ta_affy_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				dat2Affy.f3[[ts_val]]<<-dat2Affy.f2
				ta_affy_filt2[[ts_val]]<<-ta_affy_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(dat2Affy.f)!=0){
					err<-try({
						ttx<-topTable(dat2Affy.f,number=nrow(dat2Affy.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(dat2Affy.m)==TRUE)]
					dat2Affy.s<-dat2Affy.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(dat2Affy.f,number=nrow(dat2Affy.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(dat2Affy.m)==TRUE)])
						dat2Affy.s<-dat2Affy.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				dat2Affy.s2<-dat2Affy.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_affy_stat<-tclArray()
				for(i in 0:dim(dat2Affy.s2)[1]){ for(j in 0:dim(dat2Affy.s2)[2]){ 
					if(i==0){ ta_affy_stat[[i,j]]<-colnames(dat2Affy.s2)[j] } else { 
					if(j==0){ ta_affy_stat[[i,j]]<-rownames(dat2Affy.s2)[i] } else {
					tem<-dat2Affy.s2[i,j]
					ta_affy_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				dat2Affy.s3[[ts_val]]<<-dat2Affy.s2			
				ta_affy_stat2[[ts_val]]<<-ta_affy_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Agilent-OneColor"){
				datAgOne2.m<-datAgOne2.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(datAgOne2.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_Ag1[[ts_val]]<<-fit2
				datAgOne2.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_Ag1[[ts_val]]<<-samp
				ps1<-as.data.frame(datAgOne2.f)	
				rownames(ps1)<-rownames(fit2)
				datAgOne2.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_ag1_filt<-tclArray()
				for(i in 0:dim(datAgOne2.f2)[1]){ for(j in 0:dim(datAgOne2.f2)[2]){ 
					if(i==0){ ta_ag1_filt[[i,j]]<-colnames(datAgOne2.f2)[j] } else { 
					if(j==0){ ta_ag1_filt[[i,j]]<-rownames(datAgOne2.f2)[i] } else {
						tem<-datAgOne2.f2[i,j]
					ta_ag1_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				datAgOne2.f3[[ts_val]]<<-datAgOne2.f2
				ta_ag1_filt2[[ts_val]]<<-ta_ag1_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
						if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(datAgOne2.f)!=0){
					err<-try({
						ttx<-topTable(datAgOne2.f,number=nrow(datAgOne2.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datAgOne2.m)==TRUE)]
						datAgOne2.s<-datAgOne2.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datAgOne2.f,number=nrow(datAgOne2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datAgOne2.m)==TRUE)])
						datAgOne2.s<-datAgOne2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datAgOne2.s2<-datAgOne2.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_ag1_stat<-tclArray()
				for(i in 0:dim(datAgOne2.s2)[1]){ for(j in 0:dim(datAgOne2.s2)[2]){ 
					if(i==0){ ta_ag1_stat[[i,j]]<-colnames(datAgOne2.s2)[j] } else { 
					if(j==0){ ta_ag1_stat[[i,j]]<-rownames(datAgOne2.s2)[i] } else {
					tem<-datAgOne2.s2[i,j]
					ta_ag1_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				datAgOne2.s3[[ts_val]]<<-datAgOne2.s2			
				ta_ag1_stat2[[ts_val]]<<-ta_ag1_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}	
			}
			if(ts_name=="Agilent-TwoColor"){
				datAgTwo2.m<-datAgTwo2.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(datAgTwo2.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_Ag2[[ts_val]]<<-fit2
				datAgTwo2.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_Ag2[[ts_val]]<<-samp
				ps1<-as.data.frame(datAgTwo2.f)	
				rownames(ps1)<-rownames(fit2)
				datAgTwo2.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_ag2_filt<-tclArray()
				for(i in 0:dim(datAgTwo2.f2)[1]){ for(j in 0:dim(datAgTwo2.f2)[2]){ 
					if(i==0){ ta_ag2_filt[[i,j]]<-colnames(datAgTwo2.f2)[j] } else { 
					if(j==0){ ta_ag2_filt[[i,j]]<-rownames(datAgTwo2.f2)[i] } else {
					tem<-datAgTwo2.f2[i,j]
					ta_ag2_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				datAgTwo2.f3[[ts_val]]<<-datAgTwo2.f2
				ta_ag2_filt2[[ts_val]]<<-ta_ag2_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(datAgTwo2.f)!=0){
					err<-try({
							ttx<-topTable(datAgTwo2.f,number=nrow(datAgTwo2.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datAgTwo2.m)==TRUE)]
							datAgTwo2.s<-datAgTwo2.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datAgTwo2.f,number=nrow(datAgTwo2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datAgTwo2.m)==TRUE)])
						datAgTwo2.s<-datAgTwo2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datAgTwo2.s2<-datAgTwo2.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_ag2_stat<-tclArray()
				for(i in 0:dim(datAgTwo2.s2)[1]){ for(j in 0:dim(datAgTwo2.s2)[2]){ 
					if(i==0){ ta_ag2_stat[[i,j]]<-colnames(datAgTwo2.s2)[j] } else { 
					if(j==0){ ta_ag2_stat[[i,j]]<-rownames(datAgTwo2.s2)[i] } else {
					tem<-datAgTwo2.s2[i,j]
					ta_ag2_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				datAgTwo2.s3[[ts_val]]<<-datAgTwo2.s2			
				ta_ag2_stat2[[ts_val]]<<-ta_ag2_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Illumina-Beadarray"){
				datIllBA2.m<-datIllBA2.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(datIllBA2.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_Il_B[[ts_val]]<<-fit2
				datIllBA2.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_Il_B[[ts_val]]<<-samp
				datIllBA2.f<-fit2
				ps1<-as.data.frame(datIllBA2.f)	
				rownames(ps1)<-rownames(fit2)
				datIllBA2.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_ilb_filt<-tclArray()
				for(i in 0:dim(datIllBA2.f2)[1]){ for(j in 0:dim(datIllBA2.f2)[2]){ 
					if(i==0){ ta_ilb_filt[[i,j]]<-colnames(datIllBA2.f2)[j] } else { 
					if(j==0){ ta_ilb_filt[[i,j]]<-rownames(datIllBA2.f2)[i] } else {
					tem<-datIllBA2.f2[i,j]
					ta_ilb_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				datIllBA2.f3[[ts_val]]<<-datIllBA2.f2
				ta_ilb_filt2[[ts_val]]<<-ta_ilb_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(datIllBA2.f)!=0){
					err<-try({
							ttx<-topTable(datIllBA2.f,number=nrow(datIllBA2.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(datIllBA2.m1)==TRUE)]
							datIllBA2.s<-datIllBA2.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(datIllBA2.f,number=nrow(datIllBA2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(datIllBA2.m1)==TRUE)])
						datIllBA2.s<-datIllBA2.m1[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				datIllBA2.s2<-datIllBA2.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_ilb_stat<-tclArray()
				for(i in 0:dim(datIllBA2.s2)[1]){ for(j in 0:dim(datIllBA2.s2)[2]){ 
					if(i==0){ ta_ilb_stat[[i,j]]<-colnames(datIllBA2.s2)[j] } else { 
					if(j==0){ ta_ilb_stat[[i,j]]<-rownames(datIllBA2.s2)[i] } else {
					tem<-datIllBA2.s2[i,j]
					ta_ilb_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				datIllBA2.s3[[ts_val]]<<-datIllBA2.s2			
				ta_ilb_stat2[[ts_val]]<<-ta_ilb_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Illumina-Lumi"){
				lumi_NQ.m<-lumi_NQ.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(lumi_NQ.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_Il_L[[ts_val]]<<-fit2
				lumi_NQ.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_Il_L[[ts_val]]<<-samp
				ps1<-as.data.frame(lumi_NQ.f)	
				rownames(ps1)<-rownames(fit2)
				lumi_NQ.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_ill_filt<-tclArray()
				for(i in 0:dim(lumi_NQ.f2)[1]){ for(j in 0:dim(lumi_NQ.f2)[2]){ 
					if(i==0){ ta_ill_filt[[i,j]]<-colnames(lumi_NQ.f2)[j] } else { 
					if(j==0){ ta_ill_filt[[i,j]]<-rownames(lumi_NQ.f2)[i] } else {
					tem<-lumi_NQ.f2[i,j]
					ta_ill_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				lumi_NQ.f3[[ts_val]]<<-lumi_NQ.f2
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				ta_ill_filt2[[ts_val]]<<-ta_ill_filt
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(lumi_NQ.f)!=0){
					err<-try({
							ttx<-topTable(lumi_NQ.f,number=nrow(lumi_NQ.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(lumi_NQ.m)==TRUE)]
							lumi_NQ.s<-lumi_NQ.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(lumi_NQ.f,number=nrow(lumi_NQ.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(lumi_NQ.m)==TRUE)])
						lumi_NQ.s<-lumi_NQ.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				lumi_NQ.s2<-lumi_NQ.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_ill_stat<-tclArray()
				for(i in 0:dim(lumi_NQ.s2)[1]){ for(j in 0:dim(lumi_NQ.s2)[2]){ 
					if(i==0){ ta_ill_stat[[i,j]]<-colnames(lumi_NQ.s2)[j] } else { 
					if(j==0){ ta_ill_stat[[i,j]]<-rownames(lumi_NQ.s2)[i] } else {
					tem<-lumi_NQ.s2[i,j]
					ta_ill_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				lumi_NQ.s3[[ts_val]]<<-lumi_NQ.s2			
				ta_ill_stat2[[ts_val]]<<-ta_ill_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Nimblegen"){
				data.matrix_Nimblegen2.m<-data.matrix_Nimblegen2.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(data.matrix_Nimblegen2.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_N[[ts_val]]<<-fit2
				data.matrix_Nimblegen2.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_N[[ts_val]]<<-samp
				ps1<-as.data.frame(data.matrix_Nimblegen2.f)	
				rownames(ps1)<-rownames(fit2)
				data.matrix_Nimblegen2.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_nbl_filt<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.f2)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.f2)[2]){ 
					if(i==0){ ta_nbl_filt[[i,j]]<-colnames(data.matrix_Nimblegen2.f2)[j] } else { 
					if(j==0){ ta_nbl_filt[[i,j]]<-rownames(data.matrix_Nimblegen2.f2)[i] } else {
					tem<-data.matrix_Nimblegen2.f2[i,j]
					ta_nbl_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrix_Nimblegen2.f3[[ts_val]]<<-data.matrix_Nimblegen2.f2
				ta_nbl_filt2[[ts_val]]<<-ta_nbl_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(data.matrix_Nimblegen2.f)!=0){
					err<-try({
							ttx<-topTable(data.matrix_Nimblegen2.f,number=nrow(data.matrix_Nimblegen2.m))
							rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_Nimblegen2.m)==TRUE)]
							data.matrix_Nimblegen2.s<-data.matrix_Nimblegen2.m[rn,]
						},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(data.matrix_Nimblegen2.f,number=nrow(data.matrix_Nimblegen2.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_Nimblegen2.m)==TRUE)])
						data.matrix_Nimblegen2.s<-data.matrix_Nimblegen2.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				data.matrix_Nimblegen2.s2<-data.matrix_Nimblegen2.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_nbl_stat<-tclArray()
				for(i in 0:dim(data.matrix_Nimblegen2.s2)[1]){ for(j in 0:dim(data.matrix_Nimblegen2.s2)[2]){ 
					if(i==0){ ta_nbl_stat[[i,j]]<-colnames(data.matrix_Nimblegen2.s2)[j] } else { 
					if(j==0){ ta_nbl_stat[[i,j]]<-rownames(data.matrix_Nimblegen2.s2)[i] } else {
					tem<-data.matrix_Nimblegen2.s2[i,j]
					ta_nbl_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				data.matrix_Nimblegen2.s3[[ts_val]]<<-data.matrix_Nimblegen2.s2			
				ta_nbl_stat2[[ts_val]]<<-ta_nbl_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Series-Matrix"){
				data.matrixNorm.m<-data.matrixNorm.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(data.matrixNorm.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_S[[ts_val]]<<-fit2
				data.matrixNorm.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_S[[ts_val]]<<-samp
				ps1<-as.data.frame(data.matrixNorm.f)	
				rownames(ps1)<-rownames(fit2)
				data.matrixNorm.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_smt_filt<-tclArray()
				for(i in 0:dim(data.matrixNorm.f2)[1]){ for(j in 0:dim(data.matrixNorm.f2)[2]){ 
					if(i==0){ ta_smt_filt[[i,j]]<-colnames(data.matrixNorm.f2)[j] } else { 
					if(j==0){ ta_smt_filt[[i,j]]<-rownames(data.matrixNorm.f2)[i] } else {
					tem<-data.matrixNorm.f2[i,j]
					ta_smt_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrixNorm.f3[[ts_val]]<<-data.matrixNorm.f2
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				ta_smt_filt2[[ts_val]]<<-ta_smt_filt
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(data.matrixNorm.f)!=0){
					err<-try({
						ttx<-topTable(data.matrixNorm.f,number=nrow(data.matrixNorm.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrixNorm.m)==TRUE)]
						data.matrixNorm.s<-data.matrixNorm.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(data.matrixNorm.f,number=nrow(data.matrixNorm.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrixNorm.m)==TRUE)])
						data.matrixNorm.s<-data.matrixNorm.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				data.matrixNorm.s2<-data.matrixNorm.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_smt_stat<-tclArray()
				for(i in 0:dim(data.matrixNorm.s2)[1]){ for(j in 0:dim(data.matrixNorm.s2)[2]){ 
					if(i==0){ ta_smt_stat[[i,j]]<-colnames(data.matrixNorm.s2)[j] } else { 
					if(j==0){ ta_smt_stat[[i,j]]<-rownames(data.matrixNorm.s2)[i] } else {
					tem<-data.matrixNorm.s2[i,j]
					ta_smt_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				data.matrixNorm.s3[[ts_val]]<<-data.matrixNorm.s2			
				ta_smt_stat2[[ts_val]]<<-ta_smt_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}
			if(ts_name=="Online-Data"){
				data.matrix_onlineNorm.m<-data.matrix_onlineNorm.m2[[ts_val]]
pb_ma<-tkProgressBar(title="Filtering...",label="",min=0,max=1,initial=0.1,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
				err<-try(fit<-lmFit(data.matrix_onlineNorm.m[,samp],design),silent=TRUE)
				if(length(grep("incompatible dimensions",err))!=0)tkmessageBox(title="Error Comparison",message="Incompatible dimensions or wrong group comparisons",icon="error",type="ok")
				fit2<-eBayes(fit)
				fit2_O[[ts_val]]<<-fit2
				data.matrix_onlineNorm.f<-fit2
#				groups_Affy<<-groups
#				design_Affy<<-design
				samp_O[[ts_val]]<<-samp
				ps1<-as.data.frame(data.matrix_onlineNorm.f)	
				rownames(ps1)<-rownames(fit2)
				data.matrix_onlineNorm.f2<-ps1
setTkProgressBar(pb_ma,0.5,title=NULL,label="")
				ta_onl_filt<-tclArray()
				for(i in 0:dim(data.matrix_onlineNorm.f2)[1]){ for(j in 0:dim(data.matrix_onlineNorm.f2)[2]){ 
					if(i==0){ ta_onl_filt[[i,j]]<-colnames(data.matrix_onlineNorm.f2)[j] } else { 
					if(j==0){ ta_onl_filt[[i,j]]<-rownames(data.matrix_onlineNorm.f2)[i] } else {
					tem<-data.matrix_onlineNorm.f2[i,j]
					ta_onl_filt[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				data.matrix_onlineNorm.f3[[ts_val]]<<-data.matrix_onlineNorm.f2
				ta_onl_filt2[[ts_val]]<<-ta_onl_filt
setTkProgressBar(pb_ma,0.6,title="Statistical_Significant...",label="")
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Filtering"))
				}
				if(length(data.matrix_onlineNorm.f)!=0){
					err<-try({
						ttx<-topTable(data.matrix_onlineNorm.f,number=nrow(data.matrix_onlineNorm.m))
						rn<-rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_onlineNorm.m)==TRUE)]
						data.matrix_onlineNorm.s<-data.matrix_onlineNorm.m[rn,]
					},silent=TRUE)
					if(length(grep("Error",err))!=0)
					{
						ttx<-topTable(data.matrix_onlineNorm.f,number=nrow(data.matrix_onlineNorm.m))
						rn<-as.numeric(rownames(ttx)[which(rownames(ttx)%in%rownames(data.matrix_onlineNorm.m)==TRUE)])
						data.matrix_onlineNorm.s<-data.matrix_onlineNorm.m[rownames(ttx)[ttx$P.Value<=0.01],]
					}
				}
				data.matrix_onlineNorm.s2<-data.matrix_onlineNorm.s
setTkProgressBar(pb_ma,0.9,title=NULL,label="")
				ta_onl_stat<-tclArray()
				for(i in 0:dim(data.matrix_onlineNorm.s2)[1]){ for(j in 0:dim(data.matrix_onlineNorm.s2)[2]){ 
					if(i==0){ ta_onl_stat[[i,j]]<-colnames(data.matrix_onlineNorm.s2)[j] } else { 
					if(j==0){ ta_onl_stat[[i,j]]<-rownames(data.matrix_onlineNorm.s2)[i] } else {
					tem<-data.matrix_onlineNorm.s2[i,j]
					ta_onl_stat[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }			
				data.matrix_onlineNorm.s3[[ts_val]]<<-data.matrix_onlineNorm.s2			
				ta_onl_stat2[[ts_val]]<<-ta_onl_stat
setTkProgressBar(pb_ma,1,title="Completed...",label="")	
close(pb_ma)
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...Statistical_Significant"))
				}
			}	
		})
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
		sp_filt(ts_name,ts_val,children,r_choice)
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
filter_sp()
