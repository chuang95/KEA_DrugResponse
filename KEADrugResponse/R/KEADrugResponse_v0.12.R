#' Drug response prediction
#'
#' R package use for predicting drug response of patients
#'
#' @docType package
#' @name DrugResponse
NULL

#########################################
##			 Drug Response             ##
#########################################
#' predicting drug response of patients
#'
#' Returns a drug response report for each patient.
#'
#' This function predict drug response by using svm based recersive feature elimination.
#'
#' @param patient_inputfile		gene expression data for patient.
#'
#' @param inputfile_formate     input file formate, can be cel (cel file) or exp (gene expression file)
#'
#' @param patient_name     		patient name
#'
#' @param path.to.model     	where the model stored, default to '../data'
#' 
#' @return a matrix of positive or negtive score for each drug
#'
#' @export

DrugResponse.predict <- function(patient_inputfile,inputfile_formate,patient_name,path.to.model='../data',path.output='.'){

    if(!file.exists(patient_inputfile) && !dir.exists(patient_inputfile))
    	stop("cannot load pateint data: ", patient_inputfile)
    	
    cat("read in test sample gene expression data","\n")
    
    if(inputfile_formate =="cel"){
    	
    	test_geneExp=DrugResponse.readcel(patient_inputfile,patient_name,path.output)
        p_index=which(test_geneExp[,3] %in% 1)
        test_geneExp_all=t(test_geneExp[,c(1,2)])
    	test_geneExp=t(test_geneExp[p_index,c(1,2)])
    	test_geneExp=cbind(rbind('prob_id', patient_name), test_geneExp)
        test_geneExp_all=cbind(rbind('prob_id', patient_name), test_geneExp_all)
    	
    }else{
    	test_geneExp_f=file(patient_inputfile)
    	test_geneExp=read.table(test_geneExp_f)
    	test_geneExp=t(test_geneExp)
        test_geneExp_all=test_geneExp
    }
    
    
    cat("load model and cell line expression","\n")
    #path.to.model="../data"
    cellline_geneExp_f=file(paste(path.to.model,"cell_line_expression.csv",sep="/"))
    cellline_geneExp=read.csv(cellline_geneExp_f)
    cellline_score_f=file(paste(path.to.model,"celline_score.csv",sep="/"))
    cellline_score=read.csv(cellline_score_f,header=TRUE)
    #prob_gene_f=file(paste(path.to.model,"prob_gene.csv",sep="/"))
    #prob_gene =read.csv(prob_gene_f)
    
    drug1_f=file(paste(path.to.model,"prob/svm4_Gemcitabine.txt",sep="/"))
    drug1_prob=read.table(drug1_f)
    drug2_f=file(paste(path.to.model,"prob/svm4_Gefitinib.txt",sep="/"))
    drug2_prob= read.table(drug2_f)
    drug3_f=file(paste(path.to.model,"prob/svm4_Cisplatin.txt",sep="/"))
    drug3_prob= read.table(drug3_f)
    drug4_f=file(paste(path.to.model,"prob/svm4_Doxorubicin.txt",sep="/"))
    drug4_prob= read.table(drug4_f)
    drug5_f=file(paste(path.to.model,"prob/svm4_Docetaxel.txt",sep="/"))
    drug5_prob= read.table(drug5_f)
    drug6_f=file(paste(path.to.model,"prob/svm4_Paclitaxel.txt",sep="/"))
    drug6_prob= read.table(drug6_f)
    drug7_f=file(paste(path.to.model,"prob/svm4_Carboplatin.txt",sep="/"))
    drug7_prob= read.table(drug7_f)
   
    drug1wb_f=file(paste(path.to.model,"for_Rtest/Gemcitabine_wbx.txt",sep="/"))
    drug1wb_prob=read.table(drug1wb_f)
    drug2wb_f=file(paste(path.to.model,"for_Rtest/Gefitinib_wbx.txt",sep="/"))
    drug2wb_prob= read.table(drug2wb_f)
    drug3wb_f=file(paste(path.to.model,"for_Rtest/Cisplatin_wbx.txt",sep="/"))
    drug3wb_prob= read.table(drug3wb_f)
    drug4wb_f=file(paste(path.to.model,"for_Rtest/Doxorubicin_wbx.txt",sep="/"))
    drug4wb_prob= read.table(drug4wb_f)
    drug5wb_f=file(paste(path.to.model,"for_Rtest/Docetaxel_wbx.txt",sep="/"))
    drug5wb_prob= read.table(drug5wb_f)
    drug6wb_f=file(paste(path.to.model,"for_Rtest/Paclitaxel_wbx.txt",sep="/"))
    drug6wb_prob= read.table(drug6wb_f)
    drug7wb_f=file(paste(path.to.model,"for_Rtest/Carboplatin_wbx.txt",sep="/"))
    drug7wb_prob= read.table(drug7wb_f)
    
    drug_name=cbind('Gemcitabine','Gefitinib','Cisplatin','Doxorubicin','Docetaxel','Paclitaxel','Carboplatin')

	cat("normalization","\n")
	t_d=dim(test_geneExp)
	c_d=dim(cellline_geneExp)
	cellline_geneExp_sort=t(apply(cellline_geneExp,1,sort))
	cellline_median=apply(cellline_geneExp_sort,2,median)
	delta_mh=c_d[2]-t_d[2]+1;
	test_geneExp_data=matrix(as.numeric(test_geneExp[2:t_d[1],2:t_d[2]]),t_d[1]-1, t_d[2]-1)
	if(delta_mh>0){
    	n_i=floor(delta_mh/2);
    	if(t_d[1]<3){
    		test_geneExp_index=order(test_geneExp_data)
    		test_geneExp_index = test_geneExp_index +n_i;
    		test_normalized= cellline_median[test_geneExp_index]
    		test_normalized=t(test_normalized)
    	}else{
    		test_geneExp_index=t(apply(test_geneExp_data,1,order))
    		test_geneExp_index = test_geneExp_index +n_i;
    		for(i in 1:(t_d[1]-1))
    		{
    			if(i<2){
    				test_normalized= cellline_median[test_geneExp_index[i,]]
    			}else{
        		test_normalized= rbind(test_normalized, cellline_median[test_geneExp_index[i,]])
        		}
    		}
    	}
    	#
    }else{
    	n_j= ceiling(-delta_mh/2);
    	test_normalized=test_geneExp_data;
    	if(t_d[1]<3){
    		test_geneExp_index=order( test_geneExp_data[,n_j:(c_d[2]+n_j-1)] )
        	test_normalized[,n_j:(c_d[2]+n_j-1)]= cellline_median[test_geneExp_index];
    	}else{
    		test_geneExp_index=t(apply(test_geneExp_data[,n_j:(c_d[2]+n_j-1)],1,order))	
    		for(i in 1:(t_d[1]-1))
    		{
        		test_normalized[i,n_j:(c_d[2]+n_j-1)]= cellline_median[test_geneExp_index[i,]];
    		}
    	}
	}
    
    if(inputfile_formate =="cel"){
        t_d_a=dim(test_geneExp_all)
        test_geneExp_data_all=matrix(as.numeric(test_geneExp_all[2:t_d_a[1],2:t_d_a[2]]),t_d_a[1]-1, t_d_a[2]-1)
        mc=mean(cellline_median, na.rm=TRUE)
        sdc=sd(cellline_median,na.rm=TRUE)
        zt= ( test_geneExp_data_all - mean(test_geneExp_data_all, na.rm=TRUE) )/sd(test_geneExp_data_all,na.rm=TRUE)
        test_normalized_all= zt*sdc+mc
    }else{
        test_normalized_all=test_normalized
    }




    cat("calculate score for each patient","\n")
    cell_line_number=59
    patient_in_database=17
    drug_number=7
    n_p=t_d[1]-1  
    t_score=matrix(sample(0, n_p * drug_number, replace = TRUE), n_p, drug_number)# initiate score for each patient
    md=median(apply(cellline_geneExp,2,median))
    add=rep(md,t_d[1]-1)
    test_normalized =cbind(test_normalized,add)
    
    #create a matrix to plot cell line/pateint
    xc_drug=rep(1,cell_line_number)
    for(k in 2:drug_number)
    {
    	xc_drug=rbind(xc_drug,rep(k, cell_line_number))
    }
    
    xp_drug=rep(1,patient_in_database)
    for(k in 2:drug_number)
    {
    	xp_drug=rbind(xp_drug,rep(k,patient_in_database))
    }
    
    for(i in 1:n_p)
    {
    	#drug score
    	cat('Search model selected probes in patient ', patient_name,'\t', 'probe level gene expression', '\n')
		t_score[i,1]=DrugResponse.score(i, drug1wb_prob, drug1_prob, test_geneExp, test_normalized, t_d, drug_name[1], test_geneExp_all, test_normalized_all)
		t_score[i,2]=DrugResponse.score(i, drug2wb_prob, drug2_prob, test_geneExp, test_normalized, t_d, drug_name[2], test_geneExp_all, test_normalized_all)
		t_score[i,3]=DrugResponse.score(i, drug3wb_prob, drug3_prob, test_geneExp, test_normalized, t_d, drug_name[3], test_geneExp_all, test_normalized_all)
		t_score[i,4]=DrugResponse.score(i, drug4wb_prob, drug4_prob, test_geneExp, test_normalized, t_d, drug_name[4], test_geneExp_all, test_normalized_all)
		t_score[i,5]=DrugResponse.score(i, drug5wb_prob, drug5_prob, test_geneExp, test_normalized, t_d, drug_name[5], test_geneExp_all, test_normalized_all)
		t_score[i,6]=DrugResponse.score(i, drug6wb_prob, drug6_prob, test_geneExp, test_normalized, t_d, drug_name[6], test_geneExp_all, test_normalized_all)
		t_score[i,7]=DrugResponse.score(i, drug7wb_prob, drug7_prob, test_geneExp, test_normalized, t_d, drug_name[7], test_geneExp_all, test_normalized_all)

		cat("plot and save result","\n")
		cellline_score=rbind(cellline_score,t_score[i,])
		#print(toString(test_geneExp[i+1,1]))
		pdf(file.path(path.output, paste("Patient_", patient_name[i], "_drug_response_prediction.pdf",sep="")))
		xlimit=max(abs(cellline_score),na.rm = TRUE)
		dotchart(t(cellline_score[(cell_line_number+patient_in_database)+i,]),xlim=c(0-xlimit, xlimit),cex=1.5)
		points(t(cellline_score[1: cell_line_number,]), xc_drug,col="blue",cex=0.2)
		points(t(cellline_score[(cell_line_number+1):(cell_line_number+patient_in_database),]), xp_drug,col="red",cex=0.2)
		points(t(cellline_score[(cell_line_number+patient_in_database)+i,]), c(1:drug_number),cex=1.5)
		lines(c(0,0),c(0,10),pch=22,lty=2,col="grey",lwd=1.5)
		title(main=paste("Patient ", patient_name[i], " drug response prediction",sep=""), font.main=4)
		title(xlab= "Response score", col.lab=rgb(0,0.5,0))
		title(ylab= "Drugs", col.lab=rgb(0,0.5,0))
		dev.off()	
		
		#pdf(paste("Patient_", patient_name[i], "_drug_response_prediction_cellline_distribution.pdf",sep=""))
		#layout(matrix(c(1:drug_number), 4, 2, byrow = TRUE))
		#for(k in drug_number:1)
		#{
		#	p=sum(cellline_score[1:(cell_line_number+patient_in_database),k]>cellline_score[(cell_line_number+patient_in_database)+i,k],na.rm = TRUE)/(cell_line_number+patient_in_database)
		#	hist(cellline_score[1:((cell_line_number+patient_in_database)+i),k],200,xlab=paste("cell line score fraction below ",signif((1-p)*100,3),"%",sep=""),main=toString(colnames(cellline_score)[k]))
		#	abline(v = cellline_score[(cell_line_number+patient_in_database)+i,k], col = "red")
		#}
		#dev.off()
		
    }
    #output result
    d_score=data.frame(t_score)
    row.names(d_score)=patient_name
    colnames(d_score)=drug_name
    write.table(d_score, file = file.path(path.output, paste("Patient_", patient_name[i],"_drug_response_prediction_cellline_distribution.csv",sep="")), quote=F, sep=",",col.names=NA)
    return(t_score)
    cat("finished","\n")
}

#' calculate score
#'
#' Returns a positive or negtive response value for each drug.
#'
#' This function calculate drug response score by using svm model vector.
#'
#' @param test_index     patient index.
#'
#' @param drugwb_prob
#'
#' @param drug_prob
#'
#' @param test_geneExp       
#'
#' @param test_normalized          
#'
#' @param t_d       
#'
#' @param drug          
#'
#' @param test_geneExp_all       
#'
#' @param test_normalized          
#'
#' @param test_normalized_all       
#' 
#' @return a positive or negtive score for each drug
#'
#' @export
DrugResponse.score <- function(test_index, drugwb_prob, drug_prob, test_geneExp, test_normalized, t_d, drug, test_geneExp_all, test_normalized_all){
	    l=length(drugwb_prob)-2
    	drug_index=rep(0, l)
    	match=0
    	org_match=0;
        mis_match=0;
        test_for_drug=matrix(0,1,l)
        cat('Drug: ', drug, '\n')
    	for(j in 1:l)
    	{
    		t_index=which(test_geneExp[1,] %in% drug_prob[j,3])
    		if(length(t_index)==0){ 
                t_index_org=which(test_geneExp_all[1,] %in% drug_prob[j,3])
                if(length(t_index_org)==0){
                    test_for_drug[j]=test_normalized[test_index,t_d[2]]
                    cat('Could not find: ', toString(drug_prob[j,3]),'\n')
                    mis_match=mis_match+1
                }else{
                    test_for_drug[j]=test_normalized_all[test_index,(t_index_org-1)]
                    org_match=org_match+1
                }
                
    		}else{
                test_for_drug[j]=test_normalized[test_index,(t_index-1)]
    			match=match+1
    		}
    	}
    	cat('\n')
    	cat('match: ', match,'\t', 'org match: ', org_match, '\t', 'miss match: ', mis_match,'\n')
		
		
		w=drugwb_prob[1:l]
		w=t(w)
		b=drugwb_prob[l+1]
		t=as.numeric(test_for_drug)
		score=t%*%w + b
		#return(test_for_drug)
		return(-unlist(score))
}

#' read cel file
#'
#' Returns a vector of gene expression value for each patient.
#'
#' This function reads cel file.
#'
#' @param path.data
#'
#' @param patient_name
#' 
#' @return a vector of gene expression
#'
#' @export

DrugResponse.readcel <- function(path.data, patient_name, path.output){
    #read in cell file
    if (dir.exists(path.data)) {
        affy.data = ReadAffy(celfile.path = path.data)
    } else {
        affy.data = ReadAffy(filenames = path.data)
    }
    #mas5 normalize
    eset.mas5 = mas5(affy.data)
    exprSet.nologs = exprs(eset.mas5)
    exprSet = log(exprSet.nologs, 2)
    # Get the actual A/P calls
    data.mas5calls = mas5calls(affy.data)
    data.mas5calls.calls = exprs(data.mas5calls)
    x=cbind(exprSet,data.mas5calls.calls)
    pn = row.names(exprSet)
    
    dm=dim(x)
    zf=NULL
    st=dm[2]/2+1
    ap_call=matrix(0,dm[1],1)
    
    for(i in 1:dm[1])
    {
        flg=FALSE
        for(j in st:dm[2])
        {
            if(x[i,j]=="P")
            {
                flg=TRUE
            }
        }
        t=st-1
        sum=0
        for(k in 1:t)
        {
        	sum=sum+as.numeric(x[i,k])
        }
        
        if(flg==TRUE)
        {
            ap_call[i,1]=1
        }
        zf=rbind( zf, cbind(pn[i], sum/t ) )
  
    } 
    zf=cbind(zf,ap_call)
    # Print the calls as a matrix
    write.table(zf, file= file.path(path.output, paste("Patient_",patient_name,"_tp.txt",sep="")), quote=F, sep="\t", row.names = FALSE)
    return(zf)
}

#' test package
#'
#' test if a package is installed.
#'
#' This function reads cel file.
#'
#' @param x
#'
#' @export

DrugResponse.pkgTest <- function(x)
{
    if (!require(x,character.only = TRUE))
    {
        install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
}
