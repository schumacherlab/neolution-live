buildPeptideList=function(variant,peptidelength){
  # determine how many peptides can be made
  n_normal=nchar(variant$peptidecontextnormal)-(peptidelength-1)
  n_tumor=nchar(variant$peptidecontexttumor)-(peptidelength-1)
  
  # make peptides
  if (n_normal>0){
    normal_peptides=sapply(seq(1,n_normal,1),function(i) substr(x = variant$peptidecontextnormal,start = i, stop = i+(peptidelength-1)))
    c_term_pos=seq(peptidelength,nchar(variant$peptidecontextnormal),1)
    normal=data.table(normal_peptides,c_term_pos)
  } else{
    stop()
  }
  
  if (n_tumor>0){
    tumor_peptides=sapply(seq(1,n_tumor,1),function(i) substr(x = variant$peptidecontexttumor,start = i, stop = i+(peptidelength-1)))
    c_term_pos=seq(peptidelength,nchar(variant$peptidecontexttumor),1)
    tumor=data.table(tumor_peptides,c_term_pos)
  } else{
    stop()
  }
  
  # select tumor peptides != normal peptides; select corresponding normal peptides
  tumor=tumor[tumor$tumor_peptides %ni% normal$normal_peptides]
  normal=normal[match(tumor$c_term_pos,normal$c_term_pos,nomatch = FALSE)]
  
  return(list(normal,tumor))
}