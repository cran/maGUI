stat_sign<-function(h,...){

	try(({
	dat2Affy.f<-dat2Affy.f;datAgOne2.f<-datAgOne2.f;datAgTwo2.f<-datAgTwo2.f;datIllBA2.f<-datIllBA2.f;
	lumi_NQ.f<-lumi_NQ.f;data.matrix_Nimblegen2.f<-data.matrix_Nimblegen2.f;
	data.matrixNorm.f<-data.matrixNorm.f;data.matrix_onlineNorm.f<-data.matrix_onlineNorm.f;
	use.dat2Affy.m<-use.dat2Affy.m;use.datAgOne2.m<-use.datAgOne2.m;use.datAgTwo2.m<-use.datAgTwo2.m;
	use.datIllBA2.m2<-use.datIllBA2.m2;use.lumi_NQ.m<-use.lumi_NQ.m;
	use.data.matrix_Nimblegen2.m<-use.data.matrix_Nimblegen2.m;
	use.data.matrixNorm.m<-use.data.matrixNorm.m;use.data.matrix_onlineNorm.m<-use.data.matrix_onlineNorm.m;
	groups_Affy<-groups_Affy;groups_Ag1<-groups_Ag1;groups_Ag2<-groups_Ag2;groups_Il_B<-groups_Il_B;groups_Il_L<-groups_Il_L;
	groups_N<-groups_N;groups_S<-groups_S;groups_O<-groups_O;l<-l;tree<-tree;
	}),silent=TRUE)
	aa=0;bb=0;cc=0;dd=0;ee=0;fff=0;gg=0;hh=0;
	try(({
		if(exists("dat2Affy.f"))aa=length(dat2Affy.f)
		if(exists("datAgOne2.f"))bb=length(datAgOne2.f)
		if(exists("datAgTwo2.f"))cc=length(datAgTwo2.f)
		if(exists("datIllBA2.f"))dd=length(datIllBA2.f)
		if(exists("lumi_NQ.f"))ee=length(lumi_NQ.f)
		if(exists("data.matrix_Nimblegen2.f"))fff=length(data.matrix_Nimblegen2.f)
		if(exists("data.matrixNorm.f"))gg=length(data.matrixNorm.f)
		if(exists("data.matrix_onlineNorm.f"))hh=length(data.matrix_onlineNorm.f)
	}),silent=TRUE)

	dat2Affy.s=NULL;datAgOne2.s=NULL;datAgTwo2.s=NULL;datIllBA2.s=NULL;lumi_NQ.s=NULL;data.matrix_Nimblegen2.s=NULL;
	data.matrixNorm.s=NULL;data.matrix_onlineNorm.s=NULL;types=NULL
	ttx=NULL

	rm(dat2Affy.s,datAgOne2.s,datAgTwo2.s,datIllBA2.s,lumi_NQ.s,data.matrix_Nimblegen2.s,
	data.matrixNorm.s,data.matrix_onlineNorm.s,types,
	ttx)

	if(aa!=0){
		err<-try(types<<-length(unique(groups_Affy)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(dat2Affy.f,coef=types,number=nrow(use.dat2Affy.m))
		} else {
			ttx<<-topTable(dat2Affy.f,number=nrow(use.dat2Affy.m))
		}
		rn<-row.names(ttx)[ttx$P.Value<=0.01]
		err<-try(dat2Affy.s<<-use.dat2Affy.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(dat2Affy.s<<-use.dat2Affy.m[as.numeric(rn),],silent=TRUE)
		if(length(dat2Affy.s)!=0){
			visible(g1_1)<-FALSE
			l$Affymetrix$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			}
		display()
	}
	if(bb!=0){
		err<-try(types<<-length(unique(groups_Ag1)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(datAgOne2.f,coef=types,number=nrow(use.datAgOne2.m))
		} else {
			ttx<<-topTable(datAgOne2.f,number=nrow(use.datAgOne2.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(datAgOne2.s<<-use.datAgOne2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(datAgOne2.s<<-use.datAgOne2.m[as.numeric(rn),],silent=TRUE)
		if(length(datAgOne2.s)!=0){
			visible(g1_1)<-FALSE
			l$Agilent_OneColor$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()
	
	}
	if(cc!=0){
		err<-try(types<<-length(unique(groups_Ag2)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(datAgTwo2.f,coef=types,number=nrow(use.datAgTwo2.m))
		} else {
			ttx<<-topTable(datAgTwo2.f,number=nrow(use.datAgTwo2.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(datAgTwo2.s<<-use.datAgTwo2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(datAgTwo2.s<<-use.datAgTwo2.m[as.numeric(rn),],silent=TRUE)
		if(length(datAgTwo2.s)!=0){
			visible(g1_1)<-FALSE
			l$Agilent_TwoColor$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()
	}
	if(dd!=0){
		err<-try(types<<-length(unique(groups_Il_B)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(datIllBA2.f,coef=types,number=nrow(use.datIllBA2.m2))
		} else {
			ttx<<-topTable(datIllBA2.f,number=nrow(use.datIllBA2.m2))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(datIllBA2.s<<-use.datIllBA2.m2[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(datIllBA2.s<<-use.datIllBA2.m2[as.numeric(rn),],silent=TRUE)
		if(length(datIllBA2.s)!=0){
			visible(g1_1)<-FALSE
			l$Illumina_Beadarray$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()	
	}
	if(ee!=0){
		err<-try(types<<-length(unique(groups_Il_L)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(lumi_NQ.f,coef=types,number=nrow(use.lumi_NQ.m))
		} else {
			ttx<<-topTable(lumi_NQ.f,number=nrow(use.lumi_NQ.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(lumi_NQ.s<<-use.lumi_NQ.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(lumi_NQ.s<<-use.lumi_NQ.m[as.numeric(rn),],silent=TRUE)
		if(length(lumi_NQ.s)!=0){
			visible(g1_1)<-FALSE
			l$Illumina_Lumi$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()
	
	}
	if(fff!=0){
		err<-try(types<<-length(unique(groups_N)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(data.matrix_Nimblegen2.f,coef=types,number=nrow(use.data.matrix_Nimblegen2.m))
		} else {
			ttx<<-topTable(data.matrix_Nimblegen2.f,number=nrow(use.data.matrix_Nimblegen2.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(data.matrix_Nimblegen2.s<<-use.data.matrix_Nimblegen2.m[as.numeric(rn),],silent=TRUE)
		if(length(data.matrix_Nimblegen2.s)!=0){
			visible(g1_1)<-FALSE
			l$Nimblegen$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()

	}
	if(gg!=0){
		err<-try(types<<-length(unique(groups_S)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(data.matrixNorm.f,coef=types,number=nrow(use.data.matrixNorm.m))
		} else {
			ttx<<-topTable(data.matrixNorm.f,number=nrow(use.data.matrixNorm.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(data.matrixNorm.s<<-use.data.matrixNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(data.matrixNorm.s<<-use.data.matrixNorm.m[as.numeric(rn),],silent=TRUE)
		if(length(data.matrixNorm.s)!=0){
			visible(g1_1)<-FALSE
			l$Series_Matrix$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()

	}
	if(hh!=0){
		err<-try(types<<-length(unique(groups_O)),silent=TRUE)
		if(length(grep("Error",err))==0)
		{
			ttx<<-topTable(data.matrix_onlineNorm.f,coef=types,number=nrow(use.data.matrix_onlineNorm.m))
		} else {
			ttx<<-topTable(data.matrix_onlineNorm.f,number=nrow(use.data.matrix_onlineNorm.m))
		}
		rn<-rownames(ttx)[ttx$P.Value<=0.01]
		err<-try(data.matrix_onlineNorm.s<<-use.data.matrixNorm.m[rn,],silent=TRUE)
		if(length(grep("Error",err))!=0)try(data.matrix_onlineNorm.s<<-use.data.matrixNorm.m[as.numeric(rn),],silent=TRUE)
		if(length(data.matrix_onlineNorm.s)!=0){
			visible(g1_1)<-FALSE
			l$Online_Data$Stat_Significant<<-list()
			tr<<-gtree(offspring=tree,container=g1_1)
			size(tr)<-c(300,480)
			visible(g1_1)<-TRUE
			
			}
		display()
	}
}
stat_sign()
