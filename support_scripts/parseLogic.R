returnProcessedVariants=function(id,variants){
  variants$mutation_id=paste(id,1:nrow(variants),sep = "-")
  
  # rename column headers in case their differ from NKI GCF columns
  if(all(c("Gene","transcriptid","Cufflinks FPKM value (expression level)") %in% names(variants))){
    setnames(x = variants,
             old = c("Gene","transcriptid","Cufflinks FPKM value (expression level)"),
             new = c("symbol","transcript","gene_FPKM"))
  }
  
  # take subset of columns
  variantssubset=unique(x = subset(x = variants,
                                   select = c("mutation_id","symbol","transcript","peptidecontextnormal","peptidecontexttumor","gene_FPKM")),
                        by = c("peptidecontextnormal","peptidecontexttumor"))
  
  # remove any amino acid sequence after stop codon(s)
  variantssubset$peptidecontextnormal=gsub(pattern = "\\*[A-Z*]*",
                                           replacement = "",
                                           x = variantssubset$peptidecontextnormal)
  
  variantssubset$peptidecontexttumor=gsub(pattern = "\\*[A-Z*]*",
                                          replacement = "",
                                          x = variantssubset$peptidecontexttumor)
  
  # take unique subset where normal context != tumor context
  variantssubset=unique(x = subset(x = variantssubset,
                                   subset = !variantssubset$peptidecontextnormal==variantssubset$peptidecontexttumor),
                        by = c("peptidecontextnormal","peptidecontexttumor"))
  
  # if alpha characters are found in cufflinks data, set to 0 (means failed or low data)
  variantssubset$gene_FPKM[grepl(pattern = "[A-Za-z]",x = variantssubset$gene_FPKM)]=0
  
  return(variantssubset)
}