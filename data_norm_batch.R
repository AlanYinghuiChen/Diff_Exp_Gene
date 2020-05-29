#' @title read counts matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()]
#'
#' @return integer matrix
#'
#' @details In any case, first part (separated by `|`) of row names must be
#'   Ensembl transcript id
#'
#' @export


#' @title sample classinfo
#'
#' @param path string.
#'
#' @return string matrix
#'
#' @details column 1 represents sample name, column 2 represents classinfo
#'
#' @export

read_classinfo <- function(path, ...) {
  read.table(path, sep='\t', header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
}

read_matrix <- function(filename){
  read.table(filename, sep='\t', header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
}

write_matrix <- function(mat, filename){
  write.table(mat, filename, sep='\t',quote=FALSE, row.names=TRUE, col.names=TRUE)
}
################################################################################
###############################normalization####################################
################################################################################
suppressPackageStartupMessages(library(edgeR))


library(magrittr)

#' @title martix normalization
#' @examples
#' \donotrun{
#'     norm_mat(
#'         '/path/to/matrix'
#'     )
#' }
normalize <- function(
  mat,
  method,
  top_n = 20, 
  rm_gene_types = NULL,
  ref_genes=NULL
) {
  if (method == 'TMM')       mat <- norm_tmm(mat)
  else if (method == 'RLE')       mat <- norm_rle(mat)
  else if (method == 'CPM')       mat <- norm_cpm_total(mat)
  else if (method == 'CPM_top')   mat <- norm_cpm_top(mat, top_n)
  else if (method == 'null')      mat <- mat
  else stop('unknown normalization method: ', method)
  mat
}


#' @rdname  norm_scater
#'
#' @details `norm_tmm()` performs TMM normalization
#'
#' @examples
#' norm_tmm(sim_mat)
#'
#' @export
#norm_tmm <- function(mat) {
#    print('start normalization using TMM')
#    mat %>% as_SingleCellExperiment() %>%
#    {suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
#    scater::normalise() %>% SingleCellExperiment::normcounts()
#}

norm_tmm <- function(mat) {
  print('start normalization using TMM')
  dl <- edgeR::DGEList(counts=mat)
  dl <- edgeR::calcNormFactors(dl, method='TMM')
  return(edgeR::cpm(dl))
}


#' @rdname  norm_scater
#'
#' @details `norm_rle()` performs RLE normalization
#'
#' @examples
#' norm_rle(sim_mat)
#'
#' @export
#norm_rle <- function(mat) {
#    print('start normalization using RLE')
#    mat %>% as_SingleCellExperiment() %>%
#    {suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
#    scater::normalise() %>% SingleCellExperiment::normcounts()
#}
norm_rle <- function(mat) {
  print('start normalization using RLE')
  dl <- edgeR::DGEList(counts=mat)
  dl <- edgeR::calcNormFactors(dl, method='RLE')
  edgeR::cpm(dl)
}

# norm_cpm ------------------

#' @title CPM normalization by some genes
#'
#' @param mat integer matrix. counts
#' @param row integer or logical. Use which rows (genes) as normalization factor
norm_cpm_impl <- function(mat, row) {
  t(t(mat*1e6) / colSums(mat[row, , drop = F], na.rm = T))
  #edgeR::cpm(mat)
}


#' @title CPM normalization
#'
#' @description CPM normalization using counts sum of _certain_ genes as scaling factor
#'
#' @param mat integer matrix. counts.
#'
#' @details some functions may throw errors
#'
#' @family matrix normalization
#'
#' @name norm_cpm




#' @rdname norm_cpm
#'
#' @details `norm_cpm_total()` uses total genes
#'
#' @examples
#' norm_cpm_total(sim_mat)
#'
#' @export
norm_cpm_total <- function(mat) {
  print('start normalization using CPM')
  row_all <- nrow(mat) %>% seq_len()
  
  norm_cpm_impl(mat, row_all)
}

#' @rdname norm_cpm
#'
#' @param top_n integer scalar. see `norm_cpm_top()` below
#'
#' @details `norm_cpm_top()` uses top 20 genes sorted by counts (assuming `top_n = 20L`)
#'
#' @examples
#' norm_cpm_top(sim_mat, 20L)
#'
#' @export
#norm_cpm_top <- function(mat, top_n) {
#   print(paste('start normalization using top',top_n,'genes as scale factor',sep=' '))
#   if (nrow(mat) < top_n)
#    stop('too few feature for CPM top n normalization')
#    
#    row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
#    {.$ix[seq_len(top_n)]}
#    
#    norm_cpm_impl(mat, -row_top)
#}
norm_cpm_top <- function(mat, top_n) {
  print(paste('start normalization using top',top_n,'genes as scale factor',sep=' '))
  if (nrow(mat) < top_n)
    stop('too few feature for CPM top n normalization')
  
  row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
  {.$ix[seq_len(top_n)]}
  
  top = t(t(mat[row_top,]*1e6) / colSums(mat[row_top, , drop = F], na.rm = T))
  top_down= t(t(mat[setdiff(seq_len(dim(mat)[1]),row_top),]*1e6) / colSums(mat[setdiff(seq_len(dim(mat)[1]),row_top), , drop = F], na.rm = T))
  mat_top <- rbind(top,top_down)
  mat_top[rownames(mat),]
}


################################################################################
#################################batch removal##################################
################################################################################

suppressPackageStartupMessages(library(sva))


remove_batch <- function( mat, method, class_info=NULL, batch_info=NULL, ruv_k=1){
  # only remove batch for samples with batch information
  if(!is.null(batch_info)){
    samples_with_batch <- names(batch_info)[!is.na(batch_info)]
  }else{
    samples_with_batch <- colnames(mat)
  }
  if (method == 'RUV')       mat <- ruvs(mat, class_info=class_info, k=ruv_k)
  else if(method == 'RUVn')       {
    class_info <- rep(1, ncol(mat))
    names(class_info) <- colnames(mat)
    mat <- ruvs(mat, class_info=class_info, k=ruv_k)
  }
  else if (method == 'ComBat')    mat[, samples_with_batch] <- combat(mat[, samples_with_batch], class_info=class_info, batch_info=batch_info)
  else if (method == 'limma')     mat[, samples_with_batch] <- limma(mat[, samples_with_batch], class_info=class_info, batch_info=batch_info)
  else if (method == 'null')      mat <- mat
  else stop("unknown batch effect removal method: ", method)
  mat
}

ruv <- function(
  mat,
  classinfo_path,
  label_column = 2,
  k = 10
){
  suppressMessages(library(RUVSeq))
  
  print('start batch removal using RUVs')
  cIdx <- rownames(mat)
  
  sample_info <- read.table(classinfo_path,sep='\t',header=TRUE,  check.names=FALSE,  stringsAsFactors=FALSE)
  ##rank by mat
  if(unique(is.na(sample_info$sample_id))) 
    stop("sample_id not in file")
  rownames(sample_info) = sample_info$sample_id
  sample_info=sample_info[names(mat),]
  rownames(sample_info) <- c()
  
  names(sample_info)[label_column]="label"
  scIdx <- matrix(-1, ncol = max(table(sample_info$label)), nrow = dim(table(sample_info$label)))
  labellist <- names(table(sample_info$label))
  for(i in c(1:dim(table(sample_info$label)))) {
    tmp <- which(sample_info$label == labellist[i])
    scIdx[i, 1:length(tmp)] <- tmp
  }
  mat <- log(mat+0.001)
  ruv <- RUVs(as.matrix(mat), cIdx, k = k, scIdx = scIdx, isLog = TRUE)
  exp(ruv$normalizedCounts)
}

ruvs <- function(mat, class_info, batch_info=NULL, k = 1){
  if(is.null(class_info)) stop('class_info is needed for RUVs')
  
  message('start batch removal using RUVs')
  suppressMessages(library(RUVSeq))
  
  cIdx <- rownames(mat)
  
  class_sizes <- table(class_info)
  scIdx <- matrix(-1, ncol = max(class_sizes), nrow = dim(class_sizes))
  for(i in c(1:dim(class_sizes))) {
    tmp <- which(class_info == names(class_sizes)[i])
    scIdx[i, 1:length(tmp)] <- tmp
  }
  mat <- log(mat + 0.25)
  seq_ruvs <- RUVs(as.matrix(mat), cIdx, k = k, scIdx = scIdx, isLog = TRUE)
  exp(seq_ruvs$normalizedCounts)
}


combat <- function(mat,class_info=NULL, batch_info=NULL){
  if(is.null(batch_info)) stop('batch_info is needed for ComBat')
  
  message('start batch removal using combat')
  suppressMessages(library(sva))
  #batch_info <-read.table(batchinfo_path,sep='\t',row.names=1,header=T,check.names = FALSE)
  #if (!(dim(mat)[2]==dim(batch_info)[1]))
  #    stop('sample numbers in batch info and expression matrix should be same')
  #batchname <-toString(names(batch_info)[batch_column])
  #batch_info=as.data.frame(batch_info[names(mat),])
  #batch_info <- as.factor(batch_info)
  mod <- model.matrix(~ 1, data = as.factor(batch_info))
  combat <- ComBat(
    dat = log(as.matrix(mat) + 0.25),
    batch = as.factor(batch_info),
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  mat <- exp(combat)
  mat
}

limma <- function(
  mat,
  class_info=NULL,
  batch_info=NULL
){
  if(is.null(batch_info)) stop('batch_info is needed for limma')
  
  print('start batch removal using limma')
  suppressMessages(library(limma))
  
  mat <- removeBatchEffect(log(as.matrix(mat) + 0.25), as.factor(batch_info))
  mat <- exp(mat)
  mat
}

################################################################################
#################################数据计算##################################
################################################################################
data = read_matrix('D:/pythonfile/biology/pico_feature_select/data_0405/count/featurecounts.txt')
sample_ids <- colnames(data)
class_info <- NULL
batch <- NULL
class_info <- read.table('D:/pythonfile/biology/pico_feature_select/row_data/sample_classes.txt', sep='\t', header=TRUE, check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
class_info <- class_info[sample_ids, 1]
batch <- read.table('D:/pythonfile/biology/pico_feature_select/row_data/batch_info.txt', sep='\t', header=TRUE, check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
batch <- batch[sample_ids, 1]
names(batch) <- sample_ids



data_normalize = normalize(data,'TMM')

data_batch = remove_batch(data_normalize,'ComBat',class_info,batch_info=batch)
write_matrix(data_batch,'D:/pythonfile/biology/pico_feature_select/data_0405/count/filter.null.Norm_TMM.Batch_Combat_1.featurecounts.txt')

data_batch = remove_batch(data_normalize,'RUV',class_info,batch_info=batch)
write_matrix(data_batch,'D:/pythonfile/biology/pico_feature_select/data_0405/count/filter.null.Norm_TMM.Batch_RUV_1.featurecounts.txt')

data_batch = remove_batch(data_normalize,'limma',class_info,batch_info=batch)
write_matrix(data_batch,'D:/pythonfile/biology/pico_feature_select/data_0405/count/filter.null.Norm_TMM.Batch_limma_1.featurecounts.txt')


