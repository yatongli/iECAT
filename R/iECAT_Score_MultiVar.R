iECAT_MultiVar_Score<-function(Z, Y, X, internal.indicator, method, MAF.adjust){
	if(ncol(cbind(Z)) == 1){
		iECAT_SingleVar_Score(Z, Y, X, internal.indicator, method, MAF.adjust)
		stop("Single variant test")
	}
	
	if(length(Z) != length(internal.indicator)) {
		stop("Need indicator of internval vs external sources for all samples!")
	}
	
	#--- Choice of method ---#
	# "Internal", "Naive", "iECAT", "iECATminP"
	if (method=="Internal") {
		if (MAF.adjust==TRUE) {print("MAF adjustment not applicable, default to FALSE")}
		output<- MultiVar_Score_GetKernal_internal_naive(Z, Y, X, internal.indicator)
	} else if (method=="Naive") { #"Naive"
		if (MAF.adjust==TRUE) {print("MAF adjustment not applicable, default to FALSE")}
		internal.indicator=rep(1, length=length(internal.indicator))
		output<- MultiVar_Score_GetKernal_internal_naive(Z, Y, X, internal.indicator)
	} else {
		output<- MultiVar_Score_GetKernal_iECAT(Z, Y, X, internal.indicator, method, MAF.adjust)
	}
	
	
	MAF<- colMeans(Z)/2
	weight<- Beta_Weight(MAF, c(1,25))
	idx.internal<- which(internal.indicator==1)
	Cor.mat<- cor(Z[idx.internal,])
	
	re <- MultiVar_Score_GetPval(output, weight, Cor.mat)
	
	return(re)


}



MultiVar_Score_GetKernal_internal_naive<- function(Z, Y, X, internal.indicator){

	#--- Get index ---#
	idx.internal<- which(internal.indicator==1)

	#--- Create obj ---#
	data<- data.frame(cbind(Y, X))[internal.indicator,]
	G<- Z[idx.internal,]
	
	dat.null<-model.frame(Y~., data=data, na.action = na.pass)
	# there is an issue if only one column exists in dat.null
	if(ncol(dat.null)==1){
    	dat.null$add_one_column = 1
	}

	
	Null.obj<-ScoreTest_NULL_Model(Y~., data=as.data.frame(dat.null) )
	
	output <- MultiVar_Score_Kernel_internal_naive(G, Y, X, Null.obj)
	
	return(output)
}



MultiVar_Score_Kernel_internal_naive<- function(G, Y, X, Null.obj) {
	
	Sw.spa.all <- Var.Sw.spa.all <- pval.Sw.spa.all <- rep(0, ncol(G))
	
	for (jj in 1:ncol(G)) {
		#--- Testing ---#
		A1<- (G[,jj]  -  Null.obj$XXVX_inv %*%  (Null.obj$XV %*% G[,jj]))
		S1<- sum(A1 * Null.obj$res)
		Var1 <- sum(A1 *(G[,jj] * Null.obj$V) )
	
		#--- SPA and ER calibration and update variance ---#
		S1.q <- sum(A1*(Null.obj$res + Null.obj$mu))/sqrt(sum(G[,jj]))
		g1 <- A1/sqrt(sum(G[,jj]))
		mu.internal <- sum(Null.obj$mu * g1)
		var.internal <- sum(Null.obj$mu * (1-Null.obj$mu) *g1^2)
		stat_S1.q <- (S1.q-mu.internal)^2/var.internal
		p_S1.q <- pchisq(stat_S1.q, df=1, lower.tail=F)
		zscore_S1.q <- (S1.q-mu.internal)*sqrt(sum(G[,jj]))
		pnew_S1.q <- SPA_ER_pval(tempdat=as.data.frame(dat.null), G=G[,jj], q=S1.q, stat.qtemp=stat_S1.q, mu=Null.obj$mu, g=g1)
		if (pnew_S1.q>0) {Var.spa.S1 <- zscore_S1.q^2/qchisq(pnew_S1.q, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)} else{Var.spa.S1 <- Var1}
		
		Sw.spa.all[jj] <- S1
		Var.Sw.spa.all[jj] <- Var.spa.S1
		pval.Sw.spa.all[jj] <- pchisq(S1^2 / Var.spa.S1, df=1, lower.tail=FALSE)
	
	} #end of for loop jj
	
	
	output <- list()
	output$Sw.spa.all <- Sw.spa.all
	output$Var.Sw.all <- Var.Sw.all
	output$pval.Sw.all <- pval.Sw.spa.all
	
	return(output)

}


### Utility functions
Beta_Weight<-function(MAF,weights.beta){

	n<-length(MAF)
	weights<-rep(0,n)	
	IDX_0<-which(MAF == 0)
	if(length(IDX_0) == n){
		stop("No polymorphic SNPs")
	} else if( length(IDX_0) == 0){
		weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
	} else {
		weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
	}

	
	#print(length(IDX_0))
	#print(weights[-IDX_0])
	return(weights)
	
}
