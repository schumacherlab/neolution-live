readFastaFile = function(file) {
  # read data
  lines = readLines(file)

  if (length(lines) < 1) {
    stop("No input data found")
  }

  # find lines indices with name
  ind = which(substr(lines, 1L, 1L) == ">")

  # find how many sequences in data
  nseq = length(ind)
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }

  # find start(s) and end(s) of sequence data
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))

  # get sequences
  sequences = sapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))

  # get names of sequences
  seqnames = sapply(seq_len(nseq), function(i) gsub(pattern = "^> *", replacement = "", x = lines[ind[i]]))

  # put in a table
  table = data.table(sequence_id = seqnames, sequence = sequences)
  table$sequence = gsub(pattern = "\\*[A-Z*]*",
                        replacement = "",
                        x = table$sequence)

  return(table)
}

processVariants = function(sid, variants) {
  if (runParameters$verbose) message('Processing variants')

  # determine if there is any RNA expression data
  if (any(c("rna_expression", "Cufflinks FPKM value (expression level)", "gene_FPKM", "FPKM", 'tpm', "hugo_expression", "entrez_expression") %in% names(variants))) {
    ###
    # RNA EXPRESSION DATA PRESENT - determine source of data and take relevant subset
    ###
    if (all(c('variant_id', 'chromosome', 'start_position', 'end_position', 'variant_strand', 'ref_allele', 'alt_allele', 'dna_ref_read_count', 'dna_alt_read_count', 'dna_vaf',
              'rna_ref_read_count', 'rna_alt_read_count', 'rna_vaf', 'gene_id', 'transcript_id', 'transcript_strand', 'hugo_symbol', 'variant_classification',
              'transcript_remark', 'transcript_extension', 'codon_ref', 'codon_germline', 'codon_tumor', 'aa_ref', 'aa_germline', 'aa_tumor', 'aa_pos_ref',
              'aa_pos_germline', 'aa_pos_tumor_start', 'aa_pos_tumor_stop', 'protein_seq_ref', 'protein_seq_germline', 'protein_seq_tumor') %in% names(variants))) {
      # dealing with NKI varcontext (adapted from foreign antigen space project) data: rename column headers, take subset
      if (any(is.na(variants[, variant_id])) | any(variants[, variant_id] == '.') | all(variants[, variant_id] == variants[, variant_id][1])) {
        variants[, variant_id := ifelse(test = variant_id == '.' | is.na(variant_id) | variant_id == variants[, variant_id][1],
                                        yes = paste(sid, 1:nrow(variants), sep = '-'),
                                        no = paste(sid, variant_id, sep = '-'))]
      }
      variants[, chromosome := as.character(chromosome)]
      variants[, c('codon_ref', 'aa_ref', 'aa_pos_ref', 'protein_seq_ref') := NULL]

      setnames(x = variants,
               old = c('protein_seq_germline', 'protein_seq_tumor', grep('fpkm|tpm', names(variants), ignore.case = T, value = T)),
               new = c('peptidecontextnormal', 'peptidecontexttumor', 'rna_expression'))

      # don't take SNP lines along (variants are already applied in lines with tumor-specific variants)
      variantssubset = unique(x = subset(x = variants[!grepl(pattern = regexPatterns$gs_identifier, x = variant_id, perl = TRUE)],
                                         select = names(variants) %ni% c('dna_ref_read_count', 'dna_alt_read_count', 'dna_total_read_count', 'dna_vaf',
                                                                         'rna_ref_read_count', 'rna_alt_read_count', 'rna_total_read_count', 'rna_vaf', 'rna_alt_expression')),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c("variant_id", "chromosome", "start_position", "end_position", "strand", "ref_allele", "alt_allele", "ref_read_count", "alt_read_count", "vaf",
              "aa_pos_ref", "aa_pos_tumor_start", "aa_pos_tumor_stop", "protein_seq_ref", "protein_seq_tumor") %in% names(variants))) {
      # dealing with antigenic space input, proceed with all columns
      variants[, chromosome := as.character(chromosome)]
      setnames(x = variants,
               old = c('strand', "protein_seq_germline", "protein_seq_tumor"),
               new = c('transcript_strand', "peptidecontextnormal", "peptidecontexttumor"))

      variantssubset = unique(x = subset(x = variants,
                                         select = names(variants) %ni% c("codon_ref", "aa_ref", "aa_pos_ref", "protein_seq_ref", 'ref_read_count', 'alt_read_count', 'vaf')),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c("Gene", "transcriptid", "peptidecontextnormal", "peptidecontexttumor", "Cufflinks FPKM value (expression level)") %in% names(variants))) {
      # dealing with Sanger data: rename column headers, take subset
      variants[, variant_id := paste(id, 1:nrow(variants), sep = "-")]
      setnames(x = variants,
               old = c("Gene", "transcriptid", "Cufflinks FPKM value (expression level)"),
               new = c("gene_symbol", "transcript_id", "rna_expression"))

      variantssubset = unique(x = subset(x = variants,
                                         select = c("variant_id", "gene_symbol", "transcript_id", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c("symbol", "gene", "transcript", "peptidecontextnormal", "peptidecontexttumor", "gene_FPKM") %in% names(variants))) {
      # dealing with NKI kitchensink data: rename column headers, take subset
      variants[, variant_id := paste(sid, 1:nrow(variants), sep = "-")]
      setnames(x = variants,
               old = c("symbol", "gene", "transcript", "gene_FPKM"),
               new = c("gene_symbol", "gene_id","transcript_id", "rna_expression"))

      variantssubset = unique(x = subset(x = variants,
                                         select = c("variant_id", "gene_symbol", "gene_id", "transcript_id", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else {
      stop("Input format not recognized")
    }
  } else {
    ###
    # RNA EXPRESSION DATA ABSENT - determine source of data and take relevant subset
    ###
    if (all(c('variant_id', 'chromosome', 'start_position', 'end_position', 'variant_strand', 'ref_allele', 'alt_allele', 'dna_ref_read_count', 'dna_alt_read_count', 'dna_vaf',
              'rna_ref_read_count', 'rna_alt_read_count', 'rna_vaf', 'gene_id', 'transcript_id', 'transcript_strand', 'hugo_symbol', 'variant_classification',
              'transcript_remark', 'transcript_extension', 'codon_ref', 'codon_germline', 'codon_tumor', 'aa_ref', 'aa_germline', 'aa_tumor', 'aa_pos_ref',
              'aa_pos_germline', 'aa_pos_tumor_start', 'aa_pos_tumor_stop', 'protein_seq_ref', 'protein_seq_germline', 'protein_seq_tumor') %in% names(variants))) {
      # dealing with NKI new varcontext (adapted from foreign antigen space project) data: rename column headers, take subset
      if (any(is.na(variants[, variant_id])) | any(variants[, variant_id] == '.') | all(variants[, variant_id] == variants[, variant_id][1])) {
        variants[, variant_id := ifelse(test = variant_id == '.' | is.na(variant_id) | variant_id == variants[, variant_id][1],
                                        yes = paste(sid, 1:nrow(variants), sep = '-'),
                                        no = paste(sid, variant_id, sep = '-'))]
      }
      variants[, chromosome := as.character(chromosome)]
      variants[, c('codon_ref', 'aa_ref', 'aa_pos_ref', 'protein_seq_ref') := NULL]
      variants[, rna_expression := NA]

      setnames(x = variants,
               old = c('protein_seq_germline', 'protein_seq_tumor'),
               new = c('peptidecontextnormal', 'peptidecontexttumor'))

      # don't take SNP lines along (even though they shouldn't result in peptides for prediction)
      variantssubset = unique(x = subset(x = variants[!grepl(pattern = regexPatterns$gs_identifier, x = variant_id, perl = TRUE)],
                                         select = c('dna_ref_read_count', 'dna_alt_read_count', 'dna_vaf',
                                                    'rna_ref_read_count', 'rna_alt_read_count', 'rna_total_read_count', 'rna_vaf', 'rna_alt_expression')),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c("Gene", "transcriptid", "peptidecontextnormal", "peptidecontexttumor") %in% names(variants))) {
      # dealing with Sanger data: rename column headers, add "no data" for rna_expression & take subset
      variants[, variant_id := paste(sid, 1:nrow(variants), sep = "-")]
      setnames(x = variants,
               old = c("Gene", "transcriptid"),
               new = c("gene_symbol", "transcript_id"))

      variants[, rna_expression := NA]

      variantssubset = unique(x = subset(x = variants,
                                         select = c("variant_id", "gene_symbol", "transcript_id", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c("symbol", "gene", "transcript", "peptidecontextnormal", "peptidecontexttumor") %in% names(variants))) {
      # dealing with NKI kitchensink data: rename column headers, add "no data" for rna_expression & take subset
      variants[, variant_id := paste(sid, 1:nrow(variants), sep = "-")]
      setnames(x = variants,
               old = c("symbol", "gene", "transcript"),
               new = c("gene_symbol", "gene_id", "transcript_id"))

      variants[, rna_expression := NA]

      variantssubset = unique(x = subset(x = variants,
                                         select = c("variant_id", "gene_symbol", "gene_id", "transcript_id", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else if (all(c('a_full_aa_seq', 'b_full_aa_seq', 'fusion_aa_sequence') %in% colnames(variants))) {
      # dealing with structural variants
      variants[, a_chr := as.character(a_chr)]
      variants[, b_chr := as.character(b_chr)]
      variants[, rna_expression := NA]

      variantssubset = subset(x = variants,
                              select = names(variants) %ni% c("a_transcript_strand", "b_transcript_strand", "a_partner_nt_seq", "b_partner_nt_seq", "fusion_nt_sequence"))

    } else {
      stop("Input format not recognized")
    }
  }

  if (runParameters$structural_variants) {
    variantssubset$fusion_aa_sequence = gsub(pattern = "\\*[A-Z*]*",
                                             replacement = "",
                                             x = variantssubset$fusion_aa_sequence)
  } else {
    # remove stop codon and any amino acid sequence after (if present)
    variantssubset$peptidecontextnormal = gsub(pattern = "\\*[A-Z*]*",
                                               replacement = "",
                                               x = variantssubset$peptidecontextnormal)

    variantssubset$peptidecontexttumor = gsub(pattern = "\\*[A-Z*]*",
                                              replacement = "",
                                              x = variantssubset$peptidecontexttumor)

    # take unique subset where normal context != tumor context
    variantssubset = unique(x = subset(x = variantssubset,
                                       subset = variantssubset$peptidecontextnormal != variantssubset$peptidecontexttumor),
                            by = c("peptidecontextnormal", "peptidecontexttumor"))

    # clean (Sanger) expression data, set "Low confidence = 0" to 0 and set rest with alpha characters to NA (e.g. "Status: 'FAILED'|'LOW DATA'")
    if ("rna_expression" %in% names(variantssubset)) {
      variantssubset$rna_expression[grepl(pattern = "Low confidence = 0",
                                          fixed = TRUE,
                                          x = variantssubset$rna_expression)] = 0

      variantssubset$rna_expression[grepl(pattern = "[A-Za-z]+",
                                          x = variantssubset$rna_expression)] = NA

      # convert rna_expression column to numeric, as otherwise filters won't work properly
      variantssubset[, rna_expression := as.numeric(variantssubset$rna_expression)]
    }
  }
  return(variantssubset)
}
