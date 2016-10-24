buildPeptideList = function(sequences, peptidelength) {
  if (runParameters$single_sequence) {
    # determine how many peptides can be made
    n_seq = nchar(sequences$sequence) - (peptidelength - 1)

    # make peptides
    if (n_seq > 0) {
      peptide = sapply(seq(from = 1, to = n_seq, by = 1), function(i) substr(x = sequences$sequence,
                                                                             start = i,
                                                                             stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                       to = nchar(sequences$sequence),
                       by = 1)
      table = unique(x = data.table(peptide, c_term_pos),
                     by = "peptide")
      table[, sequence_id := sequences$sequence_id]
    } else {
      table = data.table()
    }

    return(table)

  } else if (runParameters$single_sequence == FALSE & runParameters$structural_variants == FALSE) {
    # determine how many peptides can be made
    n_normal = nchar(sequences$peptidecontextnormal) - (peptidelength - 1)
    n_tumor = nchar(sequences$peptidecontexttumor) - (peptidelength - 1)

    # make peptides
    if (n_normal > 0) {
      peptide = sapply(seq(from = 1, to = n_normal, by = 1), function(i) substr(x = sequences$peptidecontextnormal,
                                                                                start = i,
                                                                                stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                       to = nchar(sequences$peptidecontextnormal),
                       by = 1)
      normal = data.table(peptide, c_term_pos)
      normal = cbind(normal, subset(x = sequences, select = names(sequences) %ni% c("peptidecontextnormal", "peptidecontexttumor")))
    } else {
      normal = data.table()
    }

    if (n_tumor > 0) {
      peptide = sapply(seq(from = 1, to = n_tumor, by = 1), function(i) substr(x = sequences$peptidecontexttumor,
                                                                               start = i,
                                                                               stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                       to = nchar(sequences$peptidecontexttumor),
                       by = 1)
      tumor = unique(x = data.table(peptide, c_term_pos),
                     by = "peptide")
      tumor = cbind(tumor, subset(x = sequences, select = names(sequences) %ni% c("peptidecontextnormal", "peptidecontexttumor")))
    } else {
      tumor = data.table()
    }

    # select tumor peptides != normal peptides; select corresponding normal peptides (NOTE: in case of ins or dels corresponding normal peptide will likely be wrong)
    tumor = tumor[tumor$peptide %ni% normal$peptide]
    normal = normal[match(x = tumor$c_term_pos, table = normal$c_term_pos, nomatch = FALSE)]

    return(list(normal, tumor))
  } else if (runParameters$structural_variants) {
    # determine how many peptides can be made
    n_normal_a = nchar(sequences$a_full_aa_seq) - (peptidelength - 1)
    n_normal_b = nchar(sequences$b_full_aa_seq) - (peptidelength - 1)
    n_tumor = nchar(sequences$fusion_aa_sequence) - (peptidelength - 1)

    # make peptides
    if (n_normal_a > 0) {
      peptide = sapply(seq(from = 1, to = n_normal_a, by = 1), function(i) substr(x = sequences$a_full_aa_seq,
                                                                                start = i,
                                                                                stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                         to = nchar(sequences$a_full_aa_seq),
                         by = 1)
      normal_a = data.table(peptide, c_term_pos)
      normal_a = cbind(normal_a, subset(x = sequences, select = names(sequences) %ni% c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence")))
    } else {
      normal_a = data.table()
    }

    if (n_normal_b > 0) {
      peptide = sapply(seq(from = 1, to = n_normal_b, by = 1), function(i) substr(x = sequences$b_full_aa_seq,
                                                                                  start = i,
                                                                                  stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                       to = nchar(sequences$b_full_aa_seq),
                       by = 1)
      normal_b = data.table(peptide, c_term_pos)
      normal_b = cbind(normal_b, subset(x = sequences, select = names(sequences) %ni% c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence")))
    } else {
      normal_b = data.table()
    }

    if (n_tumor > 0) {
      peptide = sapply(seq(from = 1, to = n_tumor, by = 1), function(i) substr(x = sequences$fusion_aa_sequence,
                                                                               start = i,
                                                                               stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
                       to = nchar(sequences$fusion_aa_sequence),
                       by = 1)
      tumor = unique(x = data.table(peptide, c_term_pos),
                     by = "peptide")
      tumor = cbind(tumor, subset(x = sequences, select = names(sequences) %ni% c("a_full_aa_seq", "b_full_aa_seq", "fusion_aa_sequence")))
    } else {
      tumor = data.table()
    }

    # select tumor peptides != normal peptides; select corresponding normal peptides (NOTE: in case of ins or dels corresponding normal peptide will likely be wrong)
    tumor = tumor[tumor$peptide %ni% c(normal_a$peptide, normal_b$peptide)]

    return(list(list(normal_a, normal_b), list(tumor)))
  }
}

findVariantsContributingToEpitope = function(predicted_variants, all_variants) {
  if (nrow(predicted_variants) > 0) {
    contributing_variant_info = lapply(seq(1, nrow(predicted_variants), 1),
                                       function(x) {
                                         transcript_variants = subset(x = all_variants,
                                                                      subset = transcript_id == predicted_variants[x, ]$transcript_id)

                                         epitope_variants = unique(subset(x = transcript_variants,
                                                                          subset = sapply(seq(1, nrow(transcript_variants), 1),
                                                                                          function(y) {
                                                                                            any(c((predicted_variants[x, ]$c_term_pos - runParameters$peptidelength) : predicted_variants[x, ]$c_term_pos)[-1] >= transcript_variants$aa_pos_tumor_start[y]) &
                                                                                              any(c((predicted_variants[x, ]$c_term_pos - runParameters$peptidelength) : predicted_variants[x, ]$c_term_pos)[-1] <= transcript_variants$aa_pos_tumor_stop[y])
                                                                                          })
                                         ))

                                         # epitope_variants[, aa_pos_germline := runParameters$peptidelength - (predicted_variants[x, ]$c_term_pos - aa_pos_germline)]

                                         epitope_variants[, aa_pos_tumor_start := ifelse(test = aa_pos_tumor_start <= predicted_variants[x, ]$c_term_pos - runParameters$peptidelength + 1,
                                                                                         yes = 1,
                                                                                         no = runParameters$peptidelength - (predicted_variants[x, ]$c_term_pos - aa_pos_tumor_start))]
                                         epitope_variants[, aa_pos_tumor_stop := ifelse(test = predicted_variants[x, ]$c_term_pos <= aa_pos_tumor_stop,
                                                                                        yes = runParameters$peptidelength,
                                                                                        no = runParameters$peptidelength - (predicted_variants[x, ]$c_term_pos - aa_pos_tumor_stop))]
                                         setkey(x = epitope_variants, aa_pos_tumor_start)

                                         contributing_variants = paste(epitope_variants$variant_id,
                                                                       epitope_variants$variant_classification,
                                                                       sep = " @ ",
                                                                       collapse = "!")

                                         # contributing_aa_pos_germline = paste(epitope_variants$aa_pos_germline,
                                         #                                 collapse = ";")
                                         contributing_aa_pos_tumor = paste(ifelse(test = epitope_variants$aa_pos_tumor_start == epitope_variants$aa_pos_tumor_stop,
                                                                                  yes = epitope_variants$aa_pos_tumor_start,
                                                                                  no = paste(epitope_variants$aa_pos_tumor_start,
                                                                                             epitope_variants$aa_pos_tumor_stop,
                                                                                             sep = '-')),
                                                                           collapse = ";")

                                         return(data.table(contributing_variants = contributing_variants,
                                                           # contributing_aa_pos_germline = contributing_aa_pos_germline,
                                                           contributing_aa_pos_tumor = contributing_aa_pos_tumor))
                                       })
    return(contributing_variant_info)
  } else {
    return(list(data.table(contributing_variants = NA,
                           # contributing_aa_pos_germline = NA,
                           contributing_aa_pos_tumor = NA)))
  }
}
