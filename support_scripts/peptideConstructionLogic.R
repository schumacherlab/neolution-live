buildPeptideList=function(sequences,peptidelength){
  if (doSingleSequencePrediction){
    # determine how many peptides can be made
    n_seq=nchar(sequences$sequence)-(peptidelength-1)
    
    # make peptides
    if (n_seq>0){
      peptide=sapply(seq(1,n_seq,1),function(i) substr(x = sequences$sequence,start = i, stop = i+(peptidelength-1)))
      c_term_pos=seq(peptidelength,nchar(sequences$sequence),1)
      table=unique(data.table(peptide,c_term_pos), by = "peptide")
      table[,sequence_id := sequences$sequence_id]
    } else{
      table=data.table()
    }
    
    return(table)
    
  } else {
    # determine how many peptides can be made
    n_normal=nchar(sequences$peptidecontextnormal)-(peptidelength-1)
    n_tumor=nchar(sequences$peptidecontexttumor)-(peptidelength-1)
    
    # make peptides
    if (n_normal>0){
      peptide=sapply(seq(1,n_normal,1),function(i) substr(x = sequences$peptidecontextnormal,start = i, stop = i+(peptidelength-1)))
      c_term_pos=seq(peptidelength,nchar(sequences$peptidecontextnormal),1)
      normal=unique(data.table(peptide,c_term_pos), by = "peptide")
      normal[,variant_id := sequences$variant_id]
      normal[,gene_symbol := sequences$symbol]
      normal[,rna_expression_fpkm := sequences$gene_FPKM]
    } else{
      normal=data.table()
    }
    
    if (n_tumor>0){
      peptide=sapply(seq(1,n_tumor,1),function(i) substr(x = sequences$peptidecontexttumor,start = i, stop = i+(peptidelength-1)))
      c_term_pos=seq(peptidelength,nchar(sequences$peptidecontexttumor),1)
      tumor=unique(x = data.table(peptide,c_term_pos), by = "peptide")
      tumor[,variant_id := sequences$variant_id]
      tumor[,gene_symbol := sequences$symbol]
      tumor[,rna_expression_fpkm := sequences$gene_FPKM]
    } else{
      tumor=data.table()
    }
    
    # select tumor peptides != normal peptides; select corresponding normal peptides
    tumor=tumor[tumor$peptide %ni% normal$peptide]
    normal=normal[match(tumor$c_term_pos,normal$c_term_pos,nomatch = FALSE)]
    
    return(list(normal,tumor))
  }
}