dge<-function(h,...){
	dge_params<-function(h,h2,h3,h4,...){
		ts_name<-h;ts_val<-h2;children<-h3;r_choice=h4;
		w_imp<-tktoplevel()
		tkwm.title(w_imp,"Select your data")
		tkwm.resizable(w_imp,0,0)
		frame1_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(frame1_imp,expand=TRUE,fill="both",side="top")
		r1_imp<-c("Top genes","Gene lists")
		var1<-tclVar(r1_imp[1])
		sapply(r1_imp,function(i){
			r1_imp_vals<-ttkradiobutton(frame1_imp,variable=var1,text=i,value=i)
			tkpack(r1_imp_vals,side="left",anchor="w",padx=25)
		})
		tkpack(ttklabel(frame1_imp,text=''),side="bottom",pady=5)
		frame2_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame2_imp,text='Gene names'),side="left",padx=12)
		gn_list<-tclVar("")
		tkpack(ttkentry(frame2_imp,textvariable=gn_list,width=30),side="left",anchor="e")
		tkpack(frame2_imp,expand=TRUE,fill="both",side="top")
		frame3_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame3_imp,text='Number'),side="left",padx=12)
		numb_dge<-tclVar("all")
		tkpack(ttkentry(frame3_imp,textvariable=numb_dge,width=10),side="left",anchor="e",padx=30)
		tkpack(frame3_imp,expand=TRUE,fill="both",side="top")
		frame4_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame4_imp,text='logFC'),side="left",padx=12)
		logFC_dge<-tclVar("2")
		tkpack(ttkentry(frame4_imp,textvariable=logFC_dge,width=10),side="left",anchor="e",padx=45)
		tkpack(frame4_imp,expand=TRUE,fill="both",side="top")
		frame5_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame5_imp,text='P-value'),side="left",padx=12)
		p_list<-c(10,1,0.5,0.1,0.05,0.01,0.001,0.0001)
		p_value<-tclVar(p_list[2])
		p_comb<-ttkcombobox(frame5_imp,values=p_list,textvariable=p_value,state="normal",justify="left",width=7)
		tkpack(p_comb,side="left",padx=33)
		tkpack(frame5_imp,expand=TRUE,fill="both",side="top")
		frame6_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame6_imp,text='Adjustment'),side="left",padx=12)
		adjust_list<-c("BH","BY","holm","hochberg","hommel","bonferroni","fdr","none")
		adjust_value<-tclVar(adjust_list[1])
		adjust_comb<-ttkcombobox(frame6_imp,values=adjust_list,textvariable=adjust_value,state="normal",justify="left",width=12)
		tkpack(adjust_comb,side="left",padx=6)
		tkpack(frame6_imp,expand=TRUE,fill="both",side="top")
		frame7_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		tkpack(ttklabel(frame7_imp,text='Sort by'),side="left",padx=12)
		sort_list<-c("B","F","P","p","t","logFC","AveExpr","none")
		sort_value<-tclVar(sort_list[4])
		sort_comb<-ttkcombobox(frame7_imp,values=sort_list,textvariable=sort_value,state="normal",justify="left",width=12)
		tkpack(sort_comb,side="left",padx=35)
		tkpack(frame7_imp,expand=TRUE,fill="both",side="top")
		frame6_imp<-ttkframe(w_imp,padding=c(3,3,12,12))
		r_but<-ttkbutton(frame6_imp,text="Ok")
		tkpack(r_but,side="right",padx=12)
		tkpack(frame6_imp,expand=TRUE,fill="both",side="top")
		tkconfigure(r_but,command=function(){
			tkwm.withdraw(w_imp)
			if(ts_name=="Affymetrix"){
				dat2Affy.f<-fit2_Affy[[ts_val]]
				DE_AffyAll[[ts_val]]<<-topTable(dat2Affy.f,
				number=nrow(dat2Affy.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_Affy<-DE_AffyAll[[ts_val]]
					} else {
						DE_Affy<-topTable(dat2Affy.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_AffyAll[[ts_val]]))
					colnames(new_dge)<-names(DE_AffyAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_AffyAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_Affy<-new_dge
				}
				ta_affy_dge<-tclArray()
				for(i in 0:dim(DE_Affy)[1]){ for(j in 0:dim(DE_Affy)[2]){ 
					if(i==0){ ta_affy_dge[[i,j]]<-colnames(DE_Affy)[j] } else { 
					if(j==0){ ta_affy_dge[[i,j]]<-rownames(DE_Affy)[i] } else {
					tem<-DE_Affy[i,j]
					ta_affy_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_Affy_2[[ts_val]]<<-DE_Affy
				ta_affy_dge2[[ts_val]]<<-ta_affy_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1]){p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))}
				}
			}
			if(ts_name=="Agilent-OneColor"){
				datAgOne2.f<-fit2_Ag1[[ts_val]]
				DE_Ag1All[[ts_val]]<<-topTable(datAgOne2.f,
				number=nrow(datAgOne2.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_Ag1<-DE_Ag1All[[ts_val]]
					} else {
						DE_Ag1<-topTable(datAgOne2.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_Ag1All[[ts_val]]))
					colnames(new_dge)<-names(DE_Ag1All[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_Ag1All[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_Ag1<-new_dge
				}
				ta_ag1_dge<-tclArray()
				for(i in 0:dim(DE_Ag1)[1]){ for(j in 0:dim(DE_Ag1)[2]){ 
					if(i==0){ ta_ag1_dge[[i,j]]<-colnames(DE_Ag1)[j] } else { 
					if(j==0){ ta_ag1_dge[[i,j]]<-rownames(DE_Ag1)[i] } else {
					tem<-DE_Ag1[i,j]
					ta_ag1_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_Ag1_2[[ts_val]]<<-DE_Ag1
				ta_ag1_dge2[[ts_val]]<<-ta_ag1_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Agilent-TwoColor"){
				datAgTwo2.f<-fit2_Ag2[[ts_val]]
				DE_Ag2All[[ts_val]]<<-topTable(datAgTwo2.f,
				number=nrow(datAgTwo2.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_Ag2<-DE_Ag2All[[ts_val]]
					} else {
						DE_Ag2<-topTable(datAgTwo2.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_Ag2All[[ts_val]]))
					colnames(new_dge)<-names(DE_Ag2All[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_Ag2All[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_Ag2<-new_dge
				}
				ta_ag2_dge<-tclArray()
				for(i in 0:dim(DE_Ag2)[1]){ for(j in 0:dim(DE_Ag2)[2]){ 
					if(i==0){ ta_ag2_dge[[i,j]]<-colnames(DE_Ag2)[j] } else { 
					if(j==0){ ta_ag2_dge[[i,j]]<-rownames(DE_Ag2)[i] } else {
					tem<-DE_Ag2[i,j]
					ta_ag2_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_Ag2_2[[ts_val]]<<-DE_Ag2
				ta_ag2_dge2[[ts_val]]<<-ta_ag2_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Illumina-Beadarray"){
				datIllBA2.f<-fit2_Il_B[[ts_val]]
				DE_Il_BAll[[ts_val]]<<-topTable(datIllBA2.f,
				number=nrow(datIllBA2.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_Il_B<-DE_Il_BAll[[ts_val]]
					} else {
						DE_Il_B<-topTable(datIllBA2.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_Il_BAll[[ts_val]]))
					colnames(new_dge)<-names(DE_Il_BAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_Il_BAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_Il_B<-new_dge
				}
				ta_ilb_dge<-tclArray()
				for(i in 0:dim(DE_Il_B)[1]){ for(j in 0:dim(DE_Il_B)[2]){ 
					if(i==0){ ta_ilb_dge[[i,j]]<-colnames(DE_Il_B)[j] } else { 
					if(j==0){ ta_ilb_dge[[i,j]]<-rownames(DE_Il_B)[i] } else {
					tem<-DE_Il_B[i,j]
					ta_ilb_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_Il_B_2[[ts_val]]<<-DE_Il_B
				ta_ilb_dge2[[ts_val]]<<-ta_ilb_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Illumina-Lumi"){
				lumi_NQ.f<-fit2_Il_L[[ts_val]]
				DE_Il_LAll[[ts_val]]<<-topTable(lumi_NQ.f,
				number=nrow(lumi_NQ.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_Il_L<-DE_Il_LAll[[ts_val]]
					} else {
						DE_Il_L<-topTable(lumi_NQ.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_Il_LAll[[ts_val]]))
					colnames(new_dge)<-names(DE_Il_LAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_Il_LAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_Il_L<-new_dge
				}
				ta_ill_dge<-tclArray()
				for(i in 0:dim(DE_Il_L)[1]){ for(j in 0:dim(DE_Il_L)[2]){ 
					if(i==0){ ta_ill_dge[[i,j]]<-colnames(DE_Il_L)[j] } else { 
					if(j==0){ ta_ill_dge[[i,j]]<-rownames(DE_Il_L)[i] } else {
					tem<-DE_Il_L[i,j]
					ta_ill_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_Il_L_2[[ts_val]]<<-DE_Il_L
				ta_ill_dge2[[ts_val]]<<-ta_ill_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Nimblegen"){
				data.matrix_Nimblegen2.f<-fit2_N[[ts_val]]
				DE_NAll[[ts_val]]<<-topTable(data.matrix_Nimblegen2.f,
				number=nrow(data.matrix_Nimblegen2.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_N<-DE_NAll[[ts_val]]
					} else {
						DE_N<-topTable(data.matrix_Nimblegen2.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_NAll[[ts_val]]))
					colnames(new_dge)<-names(DE_NAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_NAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_N<-new_dge
				}
				ta_nbl_dge<-tclArray()
				for(i in 0:dim(DE_N)[1]){ for(j in 0:dim(DE_N)[2]){ 
					if(i==0){ ta_nbl_dge[[i,j]]<-colnames(DE_N)[j] } else { 
					if(j==0){ ta_nbl_dge[[i,j]]<-rownames(DE_N)[i] } else {
					tem<-DE_N[i,j]
					ta_nbl_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_N_2[[ts_val]]<<-DE_N
				ta_nbl_dge2[[ts_val]]<<-ta_nbl_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Series-Matrix"){
				data.matrixNorm.f<-fit2_S[[ts_val]]
				DE_SAll[[ts_val]]<<-topTable(data.matrixNorm.f,
				number=nrow(data.matrixNorm.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_S<-DE_SAll[[ts_val]]
					} else {
						DE_S<-topTable(data.matrixNorm.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_SAll[[ts_val]]))
					colnames(new_dge)<-names(DE_SAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_SAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_S<-new_dge
				}
				ta_smt_dge<-tclArray()
				for(i in 0:dim(DE_S)[1]){ for(j in 0:dim(DE_S)[2]){ 
					if(i==0){ ta_smt_dge[[i,j]]<-colnames(DE_S)[j] } else { 
					if(j==0){ ta_smt_dge[[i,j]]<-rownames(DE_S)[i] } else {
					tem<-DE_S[i,j]
					ta_smt_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_S_2[[ts_val]]<<-DE_S
				ta_smt_dge2[[ts_val]]<<-ta_smt_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
				}
			}
			if(ts_name=="Online-Data"){
				data.matrix_onlineNorm.f<-fit2_O[[ts_val]]
				DE_OAll[[ts_val]]<<-topTable(data.matrix_onlineNorm.f,
				number=nrow(data.matrix_onlineNorm.f),
				lfc=as.numeric(tclvalue(logFC_dge)),
				p.value=as.numeric(tclvalue(p_value)),
				adjust.method=tclvalue(adjust_value),
				sort.by=tclvalue(sort_value))
				if(tclvalue(var1)=="Top genes"){
					if(tclvalue(numb_dge)=="all"){
						DE_O<-DE_OAll[[ts_val]]
					} else {
						DE_O<-topTable(data.matrix_onlineNorm.f,
						number=as.numeric(tclvalue(numb_dge)),
						lfc=as.numeric(tclvalue(logFC_dge)),
						p.value=as.numeric(tclvalue(p_value)),
						adjust.method=tclvalue(adjust_value),
						sort.by=tclvalue(sort_value))
					}
				}
				if(tclvalue(var1)=="Gene lists"){
					xf<-strsplit(tclvalue(gn_list),",")
					new_dge<-matrix(,ncol=ncol(DE_OAll[[ts_val]]))
					colnames(new_dge)<-names(DE_OAll[[ts_val]])
					for(i in 1:length(unlist(xf))){
						new_dge<-rbind(new_dge,DE_OAll[[ts_val]][as.character(xf[[1]][i]),])
					}
					new_dge<-new_dge[-1,]
					DE_O<-new_dge
				}
				ta_onl_dge<-tclArray()
				for(i in 0:dim(DE_O)[1]){ for(j in 0:dim(DE_O)[2]){ 
					if(i==0){ ta_onl_dge[[i,j]]<-colnames(DE_O)[j] } else { 
					if(j==0){ ta_onl_dge[[i,j]]<-rownames(DE_O)[i] } else {
					tem<-DE_O[i,j]
					ta_onl_dge[[i,j]]<-ifelse(is.na(tem),".",ifelse(is.numeric(tem),round(tem,digits=5),as.character(tem)))
					} }
				} }
				DE_O_2[[ts_val]]<<-DE_O
				ta_onl_dge2[[ts_val]]<<-ta_onl_dge
				for(i in 1:length(children)){
					x<-as.character(tcl(treeview,"item",children[i],"-values"))
					if(x==r_choice[1])p_enter_class<-tcl(treeview,"insert",children[i],"end",values=as.tclObj("...DGE"))
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
		dge_params(ts_name,ts_val,children,r_choice)
	})
}
pb_ma<-tkProgressBar(title="Loading packages...",label="",min=0,max=1,initial=0,width=500)
setTkProgressBar(pb_ma,0.2,title=NULL,label="")
if(!requireNamespace("gWidgets2tcltk",quietly=TRUE)){install.packages("gWidgets2tcltk");library(gWidgets2tcltk)} else {library(gWidgets2tcltk)}
setTkProgressBar(pb_ma,0.6,title=NULL,label="")
if(!requireNamespace("tkrplot",quietly=TRUE)){install.packages("tkrplot");library(tkrplot)} else {library(tkrplot)}
setTkProgressBar(pb_ma,1,title="Done...",label="")				
close(pb_ma)
dge()
