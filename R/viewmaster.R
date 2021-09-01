#' Viewmaster
#' @description ip
#' @param query_cds cds to query
#' @param ref_cds reference cds
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @importFrom Matrix colSums
#' @export


viewmaster <-function(query_cds, 
                      ref_cds, 
                      ref_celldata_col="celltype", 
                      query_celldata_col=NULL, 
                      FUNC=c("naive_bayes", "neural_network", "bagging","softmax_regression", "logistic_regression", "deep_belief_nn"),
                      selected_genes=NULL,
                      train_frac = 0.8,
                      verbose = T){
  FUNC=match.arg(FUNC)
  common_list<-viewmaster::common_features(list(ref_cds, query_cds))

  rcds<-common_list[[1]]
  qcds<-common_list[[2]]
  
  if(is.null(selected_genes)){
    selected_common<-rownames(qcds)
  }
  
  # #no tf_idf
  query_mat<-monocle3::normalized_counts(rcds)
  ref_mat<-monocle3::normalized_counts(rcds[selected_common,])
  query_mat<-monocle3::normalized_counts(qcds[rownames(ref_mat),])

  data<-as.matrix(ref_mat)
  query<-as.matrix(query_mat)

  labf<-as.factor(colData(ref_cds)[[ref_celldata_col]])
  labn<-as.numeric(labf)-1
  labels<-levels(labf)
  laboh<-matrix(model.matrix(~0+labf), ncol = length(labels))
  colnames(laboh)<-NULL
  rownames(data)<-NULL
  colnames(data)<-NULL

  train_idx<-sample(1:dim(data)[2], round(train_frac*dim(data)[2]))
  test_idx<-which(!1:dim(data)[2] %in% train_idx)

  # dim(t(data[train_idx,]))
  # dim(t(data[test_idx,]))
  # dim(laboh[train_idx,])
  # dim(laboh[test_idx,])
  # length(labn[train_idx])
  # length(labn[test_idx])
  # length(labels)
  # dim(query)
  
 switch(FUNC, 
        naive_bayes={FUNC = naive_bayes
          funclabel="naive_bayes_"
          output = "labels"},
        neural_network={FUNC = af_nn
          funclabel="nn_"
          output = "probs"},
        softmax_regression={FUNC = smr
          funclabel="smr_"
          output = "probs"},
        deep_belief_nn={FUNC = af_dbn
          funclabel="dbnn_"
          output = "probs"},
        logisting_regression={FUNC = lr
          funclabel="lr_"
          output = "probs"},
        bagging={FUNC = lr()
          funclabel="bagging_"
          output = "probs"}
        )
        
  if(is.null(query_celldata_col)){
    coldata_label<-paste0(funclabel, "celltype")
  }else{
    coldata_label = query_celldata_col
  }
  

  if(output=="probs"){
    args<-list(data[,train_idx], 
               data[,test_idx], 
               laboh[train_idx,], 
               laboh[test_idx,], 
               length(labels), 
               query, 
               verbose = verbose)
    out<-do.call(FUNC, args)
    colnames(out)<-labels
    colData(query_cds)[[query_celldata_col]]<-colnames(as.data.frame(out))[apply(as.data.frame(out),1,which.max)]
    return(query_cds)
  }
  if(output=="labels"){
    args<-list(data[,train_idx], 
               data[,test_idx], 
               labn[train_idx], 
               labn[test_idx], 
               length(labels), 
               query, 
               verbose = verbose)
    out<-do.call(FUNC, args)
    colData(query_cds)[[coldata_label]]<-labels[out+1]
    return(query_cds)
  }
}
