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
                      FUNC=c("naive_bayes", "neural_network", "bagging","softmax_regression", "logistic_regression", "deep_belief_nn", "perceptron"),
                        selected_genes=NULL,
                      train_frac = 0.8,
                      verbose = T){
  FUNC=match.arg(FUNC)
  common_list<-viewmaster::common_features(list(ref_cds, query_cds))

  rcds<-common_list[[1]]
  qcds<-common_list[[2]]
  
  if(is.null(selected_genes)){
    selected_common<-rownames(qcds)
  }else{
    selected_common<-selected_genes
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
        logistic_regression={FUNC = lr
          funclabel="lr_"
          output = "probs"},
        bagging={FUNC = bagging
          funclabel="bagging_"
          output = "labels"},
        perceptron={FUNC = perceptron
          funclabel="perceptron_"
          output = "probs"},
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
    colData(query_cds)[[coldata_label]]<-colnames(as.data.frame(out))[apply(as.data.frame(out),1,which.max)]
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

#' Common Variant Genes
#' @description Find common variant genes between two cds objects
#' @param cds1 cds 
#' @param cds2 
#' @return a vector of similarly variant genes
#' @export


common_variant_genes <-function(cds1, 
                      cds2,
                      top_n=2000,
                      logmean_ul = 2, 
                      logmean_ll = -6,
                      row_data_column = "gene_short_name",
                      unique_data_column = "id",
                      verbose = T){
  cds1<-calculate_gene_dispersion(cds1)

  cds1<-select_genes(cds1, top_n = top_n, logmean_ul = logmean_ul, logmean_ll = logmean_ll)
  p<-plot_gene_dispersion(cds1)
  print(p)
  qsel<-rowData(cds1)[[row_data_column]][rowData(cds1)[[unique_data_column]] %in% get_selected_genes(cds1)]
  cds2<-calculate_gene_dispersion(cds2)
  p<-plot_gene_dispersion(cds2)
  print(p)
  cds2<-select_genes(cds2, top_n = top_n, logmean_ul = logmean_ul, logmean_ll = logmean_ll)
  p<-plot_gene_dispersion(cds2)
  print(p)
  rsel<-rowData(cds2)[[row_data_column]][rowData(cds2)[[unique_data_column]] %in% get_selected_genes(cds2)]
  selected_common<-intersect(qsel, rsel)
  selected_common
  }




#' Pseudo-singlets
#' @description ip
#' @param se Summarized Experiment
#' @param cds reference cds
#' @return a cell_data_set object or a list of items if unfiltered data is returned (see unfiltered)
#' @importFrom monocle3 new_cell_data_set
#' @importFrom Matrix rowSums
#' @export

pseudo_singlets <- function(se, 
                            sc_cds, 
                            assay_name="logcounts", 
                            logtransformed=T, 
                            selected_genes=NULL,
                            ncells_per_col=50, 
                            threads = detectCores()){
  if(!assay_name %in% names(assays(se))){stop("Assay name not found in Summarized Experiment object;  Run: 'names(assays(se))' to see available assays")}
  mat<-as.matrix(assays(se)[[assay_name]])
  if(logtransformed){
    cmat<-ceiling(exp(mat)-1)
  }else{
    cmat<-mat
  }
  if(is.null(selected_genes)){
    selected_genes<-rownames(se)
  }
  exprs_bulk<-cmat[selected_genes,]
  exprs_sc<-counts(sc_cds)[selected_genes[selected_genes %in% rownames(sc_cds)],]
  depth <- round(sum(rowSums(exprs_sc) / ncol(exprs_sc)))
  nRep <- 5
  n2 <- ceiling(ncells_per_col / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  message(paste0("Simulating ", (ncells_per_col * dim(se)[2]), " single cells"))
  syn_mat <- pbmcapply::pbmclapply(seq_len(ncol(exprs_bulk)), function(x){
    counts <- exprs_bulk[, x]
    counts <- rep(seq_along(as.numeric(counts)), as.numeric(counts))
    simMat <- lapply(seq_len(nRep), function(y){
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depth) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[,-1,drop=FALSE]
      simMat[,1] <- simMat[,1] + (y - 1) * n2
      simMat
    }) %>%  Reduce("rbind", .)
    simMat <- Matrix::sparseMatrix(i = simMat[,2], j = simMat[,1], x = rep(1, nrow(simMat)), dims = c(nrow(exprs_bulk), n2 * nRep))
    colnames(simMat) <- paste0(colnames(exprs_bulk)[x], "#", seq_len(ncol(simMat)))
    rownames(simMat)<-rownames(exprs_bulk)
    simMat}, mc.cores =  threads)
  syn_mat <- Reduce("cbind", syn_mat)
  if(any(is.nan(syn_mat@x))){
    syn_mat@x[is.nan(syn_mat@x)]<-0
    warning("NaN calculated during single cell generation")
  }
  slice<-rep.int(1:nrow(colData(se)), ncells_per_col)
  slice<-slice[order(slice)]
  sim_meta_data<-colData(se)[slice,]
  rownames(sim_meta_data)<-colnames(syn_mat)
  new_cell_data_set(syn_mat, sim_meta_data, DataFrame(row.names = rownames(syn_mat), gene_short_name = rownames(syn_mat), id = rownames(syn_mat)))
}



#' Seurat to Monocle3
#' @description Conver Seurat to Monocle3
#' @param seu Seurat Object
#' @param seu_rd Reduced dimname for seurat ('i.e. UMAP')
#' @param assay_name Name of data slot ('i.e. RNA')
#' @param mon_rd Reduced dimname for monocle3 ('i.e. UMAP')
#' @import Seurat
#' @import monocle3
#' @return a cell_data_set object
#' @export


seurat_to_monocle3 <-function(seu, seu_rd="umap", mon_rd="UMAP", assay_name="RNA"){
  cds<-new_cell_data_set(seu@assays[[assay_name]]@counts, 
                         cell_metadata = seu@meta.data, 
                         gene_metadata = DataFrame(
                           row.names = rownames(seu@assays[[assay_name]]@counts), 
                           id=rownames(seu@assays[[assay_name]]@counts), 
                           gene_short_name=rownames(seu@assays[[assay_name]]@counts)))
  reducedDims(cds)[[mon_rd]]<-seu@reductions[[seu_rd]]@cell.embeddings
  cds
}


