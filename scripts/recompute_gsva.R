# Optional modern rerun. Requires a complete local gene x sample expression matrix.
# Rscript scripts/recompute_gsva.R /path/to/TCGA.LUNG.sampleMap_HiSeqV2.txt output_dir
# This does NOT claim bit-for-bit equivalence with the undocumented 2023 environment.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Usage: Rscript scripts/recompute_gsva.R expression.tsv output_dir')
if (!requireNamespace('GSVA', quietly=TRUE)) stop('Install GSVA through BiocManager first')
expr <- as.matrix(read.delim(args[1],row.names=1,check.names=FALSE))
storage.mode(expr) <- 'double'
if (anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) stop('Duplicate IDs')
if (any(!is.finite(expr))) stop('Nonfinite expression values')
genes <- readLines('data/source/mir885_target_genes.txt')
sets <- list(mir885_targets=genes, mir885_targets_without_MCM5=setdiff(genes,'MCM5'))
for (f in list.files('archive/gene_sets',pattern='[.]gmt$',full.names=TRUE)) {
  for (line in readLines(f)) {
    fields <- strsplit(line,'\t',fixed=TRUE)[[1]]
    if (length(fields) >= 3) sets[[fields[1]]] <- unique(fields[-c(1,2)])
  }
}
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
coverage <- data.frame(set=names(sets),requested=lengths(sets),matched=vapply(sets,function(s)sum(s %in% rownames(expr)),integer(1)))
write.table(coverage,file=file.path(args[2],'gene_set_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Saved HiSeqV2 values are continuous; verify source normalization before interpreting a rerun.
param <- GSVA::gsvaParam(expr,sets,kcdf='Gaussian',minSize=1,maxSize=Inf)
scores <- GSVA::gsva(param,verbose=FALSE)
write.table(scores,file=file.path(args[2],'recomputed_scores.tsv'),sep='\t',quote=FALSE,col.names=NA)
writeLines(capture.output(sessionInfo()),file.path(args[2],'sessionInfo.txt'))
