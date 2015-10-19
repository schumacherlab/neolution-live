buildPeptideList=function(variant,peptidelength){
  # determine how many peptides can be made
  n_normal=nchar(variant$peptidecontextnormal)-(peptidelength-1)
  n_tumor=nchar(variant$peptidecontexttumor)-(peptidelength-1)
  
  # make peptides
  if (n_normal>0){
    peptide=sapply(seq(1,n_normal,1),function(i) substr(x = variant$peptidecontextnormal,start = i, stop = i+(peptidelength-1)))
    c_term_pos=seq(peptidelength,nchar(variant$peptidecontextnormal),1)
    normal=unique(data.table(peptide,c_term_pos), by = "peptide")
    normal[,variant_id := variant$variant_id]
    normal[,gene_symbol := variant$symbol]
    normal[,rna_expression_fpkm := variant$gene_FPKM]
  } else{
    stop()
  }
  
  if (n_tumor>0){
    peptide=sapply(seq(1,n_tumor,1),function(i) substr(x = variant$peptidecontexttumor,start = i, stop = i+(peptidelength-1)))
    c_term_pos=seq(peptidelength,nchar(variant$peptidecontexttumor),1)
    tumor=unique(x = data.table(peptide,c_term_pos), by = "peptide")
    tumor[,variant_id := variant$variant_id]
    tumor[,gene_symbol := variant$symbol]
    tumor[,rna_expression_fpkm := variant$gene_FPKM]
  } else{
    stop()
  }
  
  # select tumor peptides != normal peptides; select corresponding normal peptides
  tumor=tumor[tumor$peptide %ni% normal$peptide]
  normal=normal[match(tumor$c_term_pos,normal$c_term_pos,nomatch = FALSE)]
  
  return(list(normal,tumor))
}