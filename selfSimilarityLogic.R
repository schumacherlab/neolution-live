# load self-epitope lists
#message("********** Loading self-epitope lists")
for(i in 1:length(hlaTypes)){
  selfEpitopes[[i]]=fread(paste(scriptPath,"/selflists_v2.4/","20140610_human_proteome_",
                                hlaTypes[i],"_",
                                xMer,
                                "mer_epitopes.csv",
                                sep=""),
                          header=TRUE,
                          stringsAsFactors=FALSE)
}

# prepare some stuff for use in self-similarity function
scoreMatrix = matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
              ncol=21)
rownames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
                          'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')
colnames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G', 
                          'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')

# self-similarity check functions
matchManySequences=function(single.seq, seq.list, scoreMatrix) {
  all(sapply(seq.list, function(seq) matchSequences(single.seq, seq, scoreMatrix)$keep.in.list))
}

matchSequences=function(seq1, seq2, scoreMatrix, threshold=Inf) {
  # Split the sequences into vectors of amino acids
  peptide.list <- strsplit(c(seq1, seq2), '')
  
  # Lookup the score in 'm' function position (p1, p2, ...) 
  score.vec <- mapply(peptide.list[[1]], peptide.list[[2]], FUN=function(aa1, aa2) scoreMatrix[aa1, aa2], USE.NAMES=FALSE)
  
  # Vector of matches by position
  p.matches <- peptide.list[[1]] == peptide.list[[2]]
  
  # Return many stats
  r <- list(
    #seq1                             = seq1,
    #seq2                             = seq2,
    #change                           = paste(peptide.list[[1]], '->', peptide.list[[2]]),
    p.matches                       = p.matches,
    score.per.p                     = score.vec,
    total.score                     = sum(score.vec[c(3,4,5,6,7,8)]),
    p5.match                        = p.matches[5],
    n.mutations                     = sum(!p.matches[c(3,4,5,6,7,8)]),
    n.mutations.p3.and.p4           = sum(!p.matches[c(3,4)]),
    n.mutations.p6.p7.and.p8        = sum(!p.matches[c(6,7,8)])
  )
  r$keep.in.list <-
    (
      r$n.mutations >= 3                    # the total number of mutations within p3-p8 equals 2
      |
      r$total.score <= 5                    # the total PMBEC score equals 5 or less on p3-p8 
      |
      !r$p5.match                       # p5 is a mismatch
      |
      r$n.mutations >= 2 & r$n.mutations.p3.and.p4 == 2     # the total # of mutations equals 2 and are both located on the left side of p5
      |
      r$n.mutations >= 2 & r$n.mutations.p6.p7.and.p8 >= 2  # the total # of mutations equals 2 and are both located on the right side of p5
    )
  return(r)
}

## perform self-similarity check
performSelfSimilarityCheck=function(epitopeInputNormal,epitopeInputTumor,epitopeInputSelf){
  tumorPeptides=as.character(epitopeInputTumor$PEPTIDE)
  selfPeptides=c(as.character(epitopeInputSelf$PEPTIDE),
                  as.character(epitopeInputNormal$PEPTIDE_NORMAL)) # self-sim testing with complete human proteome peptides (500 nM & 0.5 chop score cutoff) + WT peptides from predictions list (500nM & 0.5 chop score cutoff)
  
  ## test whether peptide is similar to self
  epitopeInputTumor$NON_SELF=mclapply(tumorPeptides,
                                      matchManySequences,
                                      selfPeptides,
                                      scoreMatrix,
                                      mc.cores=(numberOfWorkers))

  epitopeInputTumor=as.data.table(mclapply(epitopeInputTumor,unlist, mc.cores=(numberOfWorkers)))
  
  # filter duplicate rows without taking SAMPLE_ID into account, use SAMPLE_ID for production code
  #epitopeInputTumor=unique(epitopeInputTumor,by=c(colnames(epitopeInputTumor)[-match("SAMPLE_ID",colnames(epitopeInputTumor))]))
  epitopeInputTumor=unique(epitopeInputTumor,by=colnames(epitopeInputTumor))
  return(epitopeInputTumor)
}