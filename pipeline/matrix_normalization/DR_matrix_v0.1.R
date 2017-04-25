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
#' Returns a positive or negtive response value for each drug.
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

DrugResponse.predict <- function(patient_inputfile,inputfile_formate,patient_name,path.to.model='../data'){

    if(!file.exists(patient_inputfile))
    	stop("cannot load pateint data: ", patient_inputfile)
    	
    cat("read in test sample gene expression data","\n")
    
    if(inputfile_formate =="cel"){
    	
    	test_geneExp=DrugResponse.readcel(patient_inputfile,patient_name)
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
    
    
    drug_name=cbind('Topotecan','Gemcitabine','Gefitinib','Cisplatin','Doxorubicin','Erlotinib','Taxotere','Taxol','Carboplatin')

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


    test_normalized_out=test_normalized_all
    test_normalized_out[p_index]=test_normalized
    test_normalized_out=t(test_normalized_out)
    colnames(test_normalized_out) = patient_name
    write.table(test_normalized_out, file= paste("Patient_",patient_name,"_tp.txt",sep=""), quote=F, sep="\t", row.names = FALSE)

    cat("finished","\n")
}

DrugResponse.readcel <- function(path.data, patient_name){
    #read in cell file
    affy.data = ReadAffy(celfile.path = path.data)
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
    
    return(zf)
}
