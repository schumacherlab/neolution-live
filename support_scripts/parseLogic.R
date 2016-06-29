readFastaFile = function(file) {
  # read data
  lines = readLines(file)

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
  # determine if there is any RNA expression data
  if (any(c("rna_expression", "Cufflinks FPKM value (expression level)", "gene_FPKM", "FPKM", "hugo_expression", "entrez_expression") %in% names(variants))) {
    ###
    # RNA EXPRESSION DATA PRESENT - determine source of data and take relevant subset
    ###
    if (all(c("donor_id", "mut_id", "chromosome", "start_position", "end_position", "strand", "ref_allele", "mut_allele", "vaf", "snp6_amp", "protein_pos_ref", "protein_pos_alt_start", "protein_pos_alt_stop", "protein_seq_ref", "protein_seq_alt") %in% names(variants))) {
      # dealing with antigenic space input, proceed with all columns
      variants[, chromosome := as.character(chromosome)]
      setnames(x = variants,
               old = c("mut_id", "protein_pos_ref", "protein_pos_alt_start", "protein_pos_alt_stop", "protein_seq_ref", "protein_seq_alt"),
               new = c("variant_id", "aa_pos_ref", "aa_pos_tumor_start", "aa_pos_tumor_stop", "peptidecontextnormal", "peptidecontexttumor"))

      variantssubset = unique(x = variants,
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
    } else if (all(c("variant_id", "chromosome", "start_position", "end_position", "ref_allele", "alt_allele", "gene_id", "transcript_id", "gene_symbol", "variant_classification", "transcript_remark", "transcript_extension", "nmd_status",
    								 "codon_ref", "codon_germline", "codon_tumor", "aa_ref", "aa_germline", "aa_tumor", "aa_pos_ref", "aa_pos_germline", "aa_pos_tumor_start", "aa_pos_tumor_stop", "protein_seq_ref", "protein_seq_germline", "protein_seq_tumor", "FPKM") %in% names(variants))) {
    	# dealing with NKI varcontext (adapted from foreign antigen space project) data: rename column headers, take subset
      if (any(is.na(variants[, variant_id])) | any(variants[, variant_id] == ".") | all(variants[, variant_id] == variants[, variant_id][1])) {
        variants[, variant_id := paste(sid, 1:nrow(variants), sep = "-")]
      }
    	variants[, chromosome := as.character(chromosome)]

    	setnames(x = variants,
    					 old = c('protein_seq_germline', 'protein_seq_tumor', 'FPKM'),
    					 new = c('peptidecontextnormal', 'peptidecontexttumor', 'rna_expression'))

    	variantssubset = unique(x = subset(x = variants,
    																		 select = c("variant_id", "gene_symbol", "gene_id", "transcript_id", "variant_classification", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
    													by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else {
      stop("Input format not recognized")
    }
  } else {
    ###
    # RNA EXPRESSION DATA ABSENT - determine source of data and take relevant subset
    ###
    if (all(c("donor_id", "mut_id", "chromosome", "start_position", "end_position", "strand", "ref_allele", "mut_allele", "vaf", "snp6_amp") %in% names(variants))) {
      # dealing with antigenic space input, proceed with all columns, add "no data" for rna_expression
      variants[, chromosome := as.character(chromosome)]
      setnames(x = variants,
               old = c("mut_id"),
               new = c("variant_id"))

      variantssubset = unique(x = variants,
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
    } else if (all(c("variant_id", "chromosome", "start_position", "end_position", "ref_allele", "alt_allele", "gene_id", "transcript_id", "gene_symbol", "variant_classification", "transcript_remark", "transcript_extension", "nmd_status",
                     "codon_ref", "codon_germline", "codon_tumor", "aa_ref", "aa_germline", "aa_tumor", "aa_pos_ref", "aa_pos_germline", "aa_pos_tumor_start", "aa_pos_tumor_stop", "protein_seq_ref", "protein_seq_germline", "protein_seq_tumor") %in% names(variants))) {
      # dealing with NKI new varcontext (adapted from foreign antigen space project) data: rename column headers, take subset
      if (any(is.na(variants[, variant_id])) | any(variants[, variant_id] == ".") | all(variants[, variant_id] == variants[, variant_id][1])) {
        variants[, variant_id := paste(sid, 1:nrow(variants), sep = "-")]
      }
      variants[, chromosome := as.character(chromosome)]

      setnames(x = variants,
               old = c('protein_seq_germline', 'protein_seq_tumor'),
               new = c('peptidecontextnormal', 'peptidecontexttumor'))

      variants[, rna_expression := NA]

      variantssubset = unique(x = subset(x = variants,
                                         select = c("variant_id", "gene_symbol", "gene_id", "transcript_id", "variant_classification", "peptidecontextnormal", "peptidecontexttumor", "rna_expression")),
                              by = c("peptidecontextnormal", "peptidecontexttumor"))
    } else {
      stop("Input format not recognized")
    }
  }

  # remove stop codon and any amino acid sequence after (if present)
  variantssubset$peptidecontextnormal = gsub(pattern = "\\*[A-Z*]*",
                                             replacement = "",
                                             x = variantssubset$peptidecontextnormal)

  variantssubset$peptidecontexttumor = gsub(pattern = "\\*[A-Z*]*",
                                            replacement = "",
                                            x = variantssubset$peptidecontexttumor)

  # take unique subset where normal context != tumor context
  variantssubset = unique(x = subset(x = variantssubset,
                                     subset = !variantssubset$peptidecontextnormal == variantssubset$peptidecontexttumor),
                          by = c("peptidecontextnormal", "peptidecontexttumor"))

  # clean (Sanger) expression data, set "Low confidence = 0" to 0 and set rest with alpha characters to NA (e.g. "Status: 'FAILED'|'LOW DATA'")
  if ("rna_expression" %in% names(variantssubset)) {
    variantssubset$rna_expression[grepl(pattern = "Low confidence = 0",
                                        fixed = TRUE,
                                        x = variantssubset$rna_expression)] = 0

    variantssubset$rna_expression[grepl(pattern = "[A-Za-z]",
                                        x = variantssubset$rna_expression)] = NA

    # convert rna_expression column to numeric, as otherwise filters won't work properly
    variantssubset[, rna_expression := as.numeric(variantssubset$rna_expression)]
  }

  return(variantssubset)
}
