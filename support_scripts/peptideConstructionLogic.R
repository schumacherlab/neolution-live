buildPeptideList=function(variant,peptidelength){
  # determine how many peptides can be made
  n_normal=nchar(variant$peptidecontextnormal)-(peptidelength-1)
  n_tumor=nchar(variant$peptidecontexttumor)-(peptidelength-1)
  
  # make peptides
  if (n_normal>0){
    normal_peptide=sapply(seq(1,n_normal,1),function(i) substr(x = variant$peptidecontextnormal,start = i, stop = i+(peptidelength-1)))
    normal_c_term_pos=seq(peptidelength,nchar(variant$peptidecontextnormal),1)
    normal=unique(data.table(normal_peptide,normal_c_term_pos), by = "normal_peptide")
  } else{
    stop()
  }
  
  if (n_tumor>0){
    tumor_peptide=sapply(seq(1,n_tumor,1),function(i) substr(x = variant$peptidecontexttumor,start = i, stop = i+(peptidelength-1)))
    tumor_c_term_pos=seq(peptidelength,nchar(variant$peptidecontexttumor),1)
    tumor=unique(x = data.table(tumor_peptide,tumor_c_term_pos), by = "tumor_peptide")
  } else{
    stop()
  }
  
  # select tumor peptides != normal peptides; select corresponding normal peptides
  tumor=tumor[tumor$tumor_peptide %ni% normal$normal_peptide]
  normal=normal[match(tumor$tumor_c_term_pos,normal$normal_c_term_pos,nomatch = FALSE)]
  
  return(list(normal,tumor))
}