# load self-epitope lists
loadSelfEpitopeList = function(path, allele, peptidelength) {
  availableSelfLists = dir(path = path,
                           pattern = paste0(allele, ".*", paste0(peptidelength,"mer")),
                           include.dirs = FALSE,
                           full.names = TRUE)

  if (length(availableSelfLists) == 1) {
    selfEpitopes = unique(x = fread(availableSelfLists,
                                    header = TRUE,
                                    stringsAsFactors = FALSE),
                          by = "peptide")
    setnames(x = selfEpitopes,
             old = names(selfEpitopes),
             new = tolower(names(selfEpitopes)))
    write(x = paste0("Self-epitope list:\t\t\t", basename(availableSelfLists[1]),"\n"),
          file = file.path(runParameters$filepath,
                           'predictions_logs',
                           paste0(paste(runStart,
                                        runParameters$filename_no_ext,
                                        runParameters$allele,
                                        runParameters$peptidelength, sep = '_'),
                                  'mer_runInfo.txt')),
          append = TRUE)
    return(selfEpitopes)
  } else {
    stop(paste0("Zero or more than one self-epitope lists found,
                please make sure only 1 list per HLA & peptide length is present in ", path))
  }
}

# prepare matrix used in self-similarity function
loadSelfSimilarityMatrix = function() {
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
                       ncol = 21)
  rownames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
                            'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')
  colnames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
                            'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')
  return(scoreMatrix)
}

# simple self-similarity check
performSimpleSelfSimilarityCheck = function(epitopes, selfepitopes, scorematrix, normalepitopes = list()) {
  if (length(epitopes) < 1) {
    return(c(""))
  }

  selfepitopes = c(as.character(selfepitopes), # add complete human proteome epitopes
                   as.character(normalepitopes[!is.na(normalepitopes)]))  # add normal epitopes from predictions list, when available

  ## test whether peptide is similar to self
  different_from_self = mclapply(X = epitopes,
                                 FUN = matchManySequencesSimple,
                                 seq.list = selfepitopes,
                                 scorematrix = scorematrix,
                                 mc.cores = numberOfWorkers)

  different_from_self = unlist(different_from_self)

  return(different_from_self)
}

# extended self-similarity check using predicted human proteome epitopes (and normal epitopes)
performExtendedSelfSimilarityCheck = function(epitopes, selfepitopes, scorematrix, normalepitopes = list()) {
  if (length(epitopes) < 1) {
    return(c(""))
  }

  selfepitopes = c(as.character(selfepitopes), # add complete human proteome epitopes
                   as.character(normalepitopes[!is.na(normalepitopes)]))  # add normal epitopes from predictions list, when available

  ## test whether peptide is similar to self
  different_from_self = mclapply(X = epitopes,
                                 FUN = matchManySequencesExtended,
                                 seq.list = selfepitopes,
                                 scorematrix = scorematrix,
                                 mc.cores = numberOfWorkers)

  different_from_self = unlist(different_from_self)

  return(different_from_self)
}

# simple self-similarity check
matchManySequencesSimple = function(single.seq, seq.list, scorematrix) {
  all(sapply(seq.list, function(seq) matchSequencesSimpleCompiled(seq1 = single.seq,
                                                                  seq2 = seq,
                                                                  scorematrix = scorematrix)$keep.in.list))
}

matchSequencesSimple = function(seq1, seq2, scorematrix, threshold = Inf) {
  # Split the sequences into vectors of amino acids
  peptide.list <- strsplit(x = c(seq1, seq2),
                           split = '')

  # Lookup the score in 'm' function position (p1, p2, ...)
  score.vec <- mapply(peptide.list[[1]], peptide.list[[2]],
                      FUN = function(aa1, aa2) scorematrix[aa1, aa2],
                      USE.NAMES = FALSE)

  # Vector of matches by position
  p.matches <- peptide.list[[1]] == peptide.list[[2]]

  # Return many stats.
  r <- list(
    # seq1          = seq1,
    # seq2          = seq2,
    # change        = paste(peptide.list[[1]], '->', peptide.list[[2]]),
    p.matches     = p.matches,
    score.per.p   = score.vec,
    n.mutations   = sum(!p.matches[c(3:(runParameters$peptidelength - 1))]),
    total.score   = sum(score.vec[c(3:(runParameters$peptidelength - 1))])
  )

  r$keep.in.list <-
    (
      r$n.mutations >= 2		# the total number of mutations within the core peptide = 2 or more
      | r$total.score <= (runParameters$peptidelength - 4)  # the total PMBEC score equals 5, 6 or 7 for 9-, 10- or 11-mers, repectively
    )
  return(r)
}

# extended self-similarity check
matchManySequencesExtended = function(single.seq, seq.list, scorematrix) {
  all(sapply(seq.list, function(seq) matchSequencesExtendedCompiled(seq1 = single.seq,
                                                                    seq2 = seq,
                                                                    scorematrix = scorematrix)$keep.in.list))
}

matchSequencesExtended = function(seq1, seq2, scorematrix, threshold = Inf) {
  # Split the sequences into vectors of amino acids
  peptide.list <- strsplit(x = c(seq1, seq2),
                           split = '')

  # Lookup the score in 'm' function position (p1, p2, ...)
  score.vec <- mapply(peptide.list[[1]], peptide.list[[2]],
                      FUN = function(aa1, aa2) scorematrix[aa1, aa2],
                      USE.NAMES = FALSE)

  # Vector of matches by position
  p.matches <- peptide.list[[1]] == peptide.list[[2]]

  # Return many stats
  r <- list(
    # seq1                             = seq1,
    # seq2                             = seq2,
    # change                           = paste(peptide.list[[1]], '->', peptide.list[[2]]),
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
      r$n.mutations >= 3     # the total number of mutations within p3-p8 equals 3 or more
      | r$total.score <= 5   # the total PMBEC score equals 5 or less on p3-p8
      | !r$p5.match          # p5 is a mismatch
      | r$n.mutations >= 2 & r$n.mutations.p3.and.p4 == 2     # the total # of mutations equals 2 and are both located on the left side of p5
      | r$n.mutations >= 2 & r$n.mutations.p6.p7.and.p8 >= 2  # the total # of mutations equals 2 and are both located on the right side of p5
    )
  return(r)
}

# compile self-sim functions
matchSequencesSimpleCompiled = cmpfun(matchSequencesSimple)
matchSequencesExtendedCompiled = cmpfun(matchSequencesExtended)
