gen_pep_table <- function(N_start_residues, context, peptidelength, sequences) {
  if (N_start_residues <= 0) {
    return(data.table)
  }

  dtf <- data.table(
    peptide = vapply(seq(1, N_start_residues, by=1),
      function(i) substr(x = context,
        start = i, stop = i + (peptidelength - 1)), character(1)))
  dtf[, 'c_term_pos' :=
    seq(from = peptidelength, to = nchar(context))]

  dtf <- cbind(dtf,
    sequences[, setdiff(colnames(sequences),
      c('peptidecontextnormal', 'peptidecontexttumor')), with = F])

  return(dtf)
}


buildPeptideList <- function(sequences, peptidelength, runParameters) {
  if (runParameters$single_sequence) {
    ## Determine how many peptides can be made
    N_residues <- stringr::str_length(sequences$sequence)
    N_seq <- N_residues - peptidelength + 1

    if (N_seq > 0) {
      peptide <- sapply(seq(from=1, to=N_seq, by=1),
        function(i) substr(sequences$sequence, i, i + (peptidelength - 1)))
      c_term_pos <- seq(from=peptidelength, to=N_residues, by=1)
      ## 2019-01-11 19:28 Noticed a problem here, unique was being run on
      ## peptide only
      dtf <- data.table(peptide, c_term_pos)
      dtf[, sequence_id := sequences$sequence_id]
    } else {
      dtf <- data.table()
    }
    return(dtf)
  } else if (runParameters$structural_variants == FALSE) {
    ## Determine how many peptides can be made
    N_residues_normal <- stringr::str_length(sequences$peptidecontextnormal)
    N_residues_tumor <- stringr::str_length(sequences$peptidecontexttumor)
    N_normal <- N_residues_normal - (peptidelength - 1)
    N_tumor <- N_residues_tumor - (peptidelength - 1)
    normal <- gen_pep_table(N_normal, sequences$peptidecontextnormal,
      peptidelength, sequences)
    tumor <- gen_pep_table(N_tumor, sequences$peptidecontexttumor,
      peptidelength, sequences)

    ## select tumor peptides != normal peptides; select corresponding normal
    ## peptides (NOTE: in case of ins or dels corresponding normal peptide will
    ## likely be wrong)
    idxs <- tumor$peptide %ni% normal$peptide
    normal_idxs <-
      match(x = tumor$c_term_pos, table = normal$c_term_pos, nomatch = FALSE)
    tumor <- tumor[idxs]
    normal <- normal[normal_idxs]
    return(list(normal, tumor))
  } else if (runParameters$structural_variants) {
    # determine how many peptides can be made
    n_normal_a <- nchar(sequences$a_full_aa_seq) - (peptidelength - 1)
    n_normal_b <- nchar(sequences$b_full_aa_seq) - (peptidelength - 1)
    n_tumor <- nchar(sequences$fusion_aa_sequence) - (peptidelength - 1)

    # make peptides
    if (n_normal_a > 0) {
      peptide = sapply(seq(from = 1, to = n_normal_a, by = 1),
        function(i) substr(x = sequences$a_full_aa_seq,
          start = i, stop = i + (peptidelength - 1)))
      c_term_pos = seq(from = peptidelength,
        to = nchar(sequences$a_full_aa_seq), by = 1)
      normal_a = data.table(peptide, c_term_pos)
      normal_a = cbind(normal_a, subset(x = sequences,
          select = names(sequences) %ni%
            c('a_full_aa_seq', 'b_full_aa_seq', 'fusion_aa_sequence'))) } else {
      normal_a = data.table()
    }

    if (n_normal_b > 0) {
      peptide <- sapply(seq(from = 1, to = n_normal_b, by = 1),
        function(i) substr(x = sequences$b_full_aa_seq,
          start = i, stop = i + (peptidelength - 1)))
      c_term_pos <- seq(from = peptidelength,
        to = nchar(sequences$b_full_aa_seq),
        by = 1)
      normal_b <- data.table(peptide, c_term_pos)
      normal_b <- cbind(normal_b, subset(x = sequences,
          select = names(sequences) %ni%
            c('a_full_aa_seq', 'b_full_aa_seq', 'fusion_aa_sequence')))
    } else {
      normal_b = data.table()
    }

    if (n_tumor > 0) {
      peptide <- sapply(seq(from = 1, to = n_tumor, by = 1),
        function(i) substr(x = sequences$fusion_aa_sequence,
          start = i, stop = i + (peptidelength - 1)))
      c_term_pos <- seq(from = peptidelength,
        to = nchar(sequences$fusion_aa_sequence),
        by = 1)
      dumor <- data.table(peptide, c_term_pos)
      # tumor <- unique(tumor, by = 'peptide')
      tumor <- cbind(tumor,
        subset(x = sequences, select = names(sequences) %ni%
          c('a_full_aa_seq', 'b_full_aa_seq', 'fusion_aa_sequence')))
    } else {
      tumor <- data.table()
    }

    # select tumor peptides != normal peptides
    tumor = tumor[tumor$peptide %ni% c(normal_a$peptide, normal_b$peptide)]

    return(list(list(normal_a, normal_b), list(tumor)))
  }
}

findVariantsContributingToEpitope <- function(predicted_variants, all_variants,
  runParameters) {
  if (runParameters$verbose) message('Determining contributing variants')

  if (nrow(predicted_variants) > 0) {
    contributing_variant_info = lapply(seq(1, nrow(predicted_variants), 1),
      function(x) {
        transcript_variants = subset(x = all_variants,
          subset = transcript_id == predicted_variants[x, ]$transcript_id)

        epitope_variants = unique(subset(x = transcript_variants,
            subset = sapply(seq(1, nrow(transcript_variants), 1),
              function(y) {
                any(c((predicted_variants[x, ]$c_term_pos - runParameters$peptidelength) :
                    predicted_variants[x, ]$c_term_pos)[-1] >=
                  transcript_variants$aa_pos_tumor_start[y]) &
                  any(c((predicted_variants[x, ]$c_term_pos - runParameters$peptidelength) :
                      predicted_variants[x, ]$c_term_pos)[-1] <=
                    transcript_variants$aa_pos_tumor_stop[y])
              })
            ))


        epitope_variants[, aa_pos_tumor_start :=
          ifelse(test = aa_pos_tumor_start <= predicted_variants[x, ]$c_term_pos -
            runParameters$peptidelength + 1,
          yes = 1,
          no = runParameters$peptidelength - (predicted_variants[x, ]$c_term_pos -
            aa_pos_tumor_start))]
        epitope_variants[, aa_pos_tumor_stop :=
          ifelse(test = predicted_variants[x, ]$c_term_pos <= aa_pos_tumor_stop,
            yes = runParameters$peptidelength,
            no = runParameters$peptidelength -
              (predicted_variants[x, ]$c_term_pos - aa_pos_tumor_stop))]
        setkey(x = epitope_variants, aa_pos_tumor_start)

        contributing_variants <-
          paste(epitope_variants$variant_id,
            epitope_variants$variant_classification, sep = ' @ ', collapse = '!')

        contributing_aa_pos_tumor <- paste(ifelse(test =
            epitope_variants$aa_pos_tumor_start == epitope_variants$aa_pos_tumor_stop,
          yes = epitope_variants$aa_pos_tumor_start,
          no = paste(epitope_variants$aa_pos_tumor_start,
            epitope_variants$aa_pos_tumor_stop, sep = '-')),
        collapse = ';')

        if ('rna_alt_expression' %in% names(epitope_variants)) {
          contributing_variants_alt_expression = paste(epitope_variants$rna_alt_expression,
            collapse = '!')

          data <- data.table(contributing_variants = contributing_variants,
            contributing_aa_pos_tumor = contributing_aa_pos_tumor,
            contributing_variants_alt_expression = contributing_variants_alt_expression)
        } else {
          data <- data.table(contributing_variants = contributing_variants,
            contributing_aa_pos_tumor = contributing_aa_pos_tumor)
        }

        return(data)
      })
    return(contributing_variant_info)
  } else {
    if ('rna_alt_expression' %in% names(all_variants)) {
      data <- data.table(contributing_variants = NA,
                        # contributing_aa_pos_germline = NA,
                        contributing_aa_pos_tumor = NA,
                        contributing_variants_alt_expression = NA)
    } else {
      data <- data.table(contributing_variants = NA,
                        # contributing_aa_pos_germline = NA,
                        contributing_aa_pos_tumor = NA)
    }
    return(list(data))
  }
}
