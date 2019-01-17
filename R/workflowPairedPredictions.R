performPairedSequencePredictions <- function(runParameters, unique_cols) {
  if (runParameters$verbose) {
    message('Starting paired sequence workflow')
  }

  ## import data & clean up
  variantInput <- fread(
    input = file.path(runParameters$filepath, runParameters$filename),
    na.strings = c('NA', '', '-')) %>% unique

  sampleId <- gsub(
      ## remove some seq facility-specific info
      pattern = paste('_[ATCG]{7,8}|_S1_LMG_TR1_001|_L[0-9]{3}_t',
      ## substitute various suffixes
      paste0('[-_]tumor|[-_]kitchensink|[-_]mdup|[-_]ra|[-_]bq|[-_]table',
        '|[-_]complete|[-_]merged|[-_]varcontext|[-_]rna'),
      '|\\.txt.*$|\\.vcf.*$|\\.[ct]sv.*', # substitute various extensions
      '^[0-9]{8}_|^[0-9]{8}[-_][0-9]{4}_', # substitute various prefixes
      sep = '|'),
    replacement = '', x = runParameters$filename_no_ext, perl = T) %>% toupper

  variantInfo <- processVariants(sid = sampleId, variants = variantInput,
    runParameters = runParameters)

  ## write input to disk
  if (runParameters$copyinput) {
    write.table(x = variantInput,
      file = file.path(runParameters$filepath, 'predictions_input',
        paste0(runParameters$filename_no_ext, '_input.tsv')),
      quote = FALSE,
      sep = '\t',
      row.names = FALSE)
  }

  ## prepare vectors with colnames and colclasses for making empty tables, in
  ## case needed
  columnNamesEmptyTable <- c(
      names(variantInfo) %>%
      setdiff(c('peptidecontextnormal', 'peptidecontexttumor')),
    'c_term_pos', 'hla_allele', 'xmer', 'tumor_peptide', 'tumor_c_term_aa',
    paste0('tumor_', runParameters$allele, 'affinity'),
    paste0('tumor_', runParameters$allele, 'percentile_rank'),
    'tumor_processing_score', 'normal_peptide', 'normal_c_term_aa',
    paste0('normal_', runParameters$allele, 'affinity'),
    paste0('normal_', runParameters$allele, 'percentile_rank'),
    'normal_processing_score')

  columnClassesEmptyTable <- c(
    sapply(variantInfo, class) %>%
      .[names(.) %ni% c('peptidecontextnormal', 'peptidecontexttumor')] %>%
      setNames(NULL),
    'numeric', 'character', 'numeric', 'character', 'character', 'numeric',
    'numeric', 'numeric', 'character', 'character', 'numeric', 'numeric',
    'numeric')

  ## write empty output if no variants found
  if (nrow(variantInfo) == 0) {
    writePredictionsToDisk(
      runParameters = runParameters,
      dtf = emptyTableWithColumnNamesAndColumnClasses(
        colnames = columnNamesEmptyTable,
        colclasses = columnClassesEmptyTable),
      unique_cols = unique_cols,
      filepath = runParameters$filepath,
      filename = runParameters$filename_no_ext,
      allele = runParameters$allele,
      peptidelength = runParameters$peptidelength)
    return(NULL)
  }

  ## load required data for self-similarity check
  if (runParameters$simple_selfsim || runParameters$extended_selfsim) {
    if (runParameters$use_selflist) {
      selfEpitopes <- loadSelfEpitopeList(
        path = runParameters$selfEpitopeListPath,
        allele = runParameters$allele,
        peptidelength = runParameters$peptidelength)
    } else {
      selfEpitopes <- emptyTableWithColumnNamesAndColumnClasses(
        colnames = c('sequence_id', 'hla_allele', 'xmer', 'peptide', 'c_term_aa',
          'c_term_pos', paste0(runParameters$allele,'affinity'),
          'processing_score'),
        colclasses = c('character', 'character', 'numeric',
          'character', 'character', 'numeric',
          'numeric', 'numeric'))
    }
    scoreMatrix <- loadSelfSimilarityMatrix()
  }

  ## For each variant line, list tumor peptides which are different from
  ## normal (and corresponding normal peptides) and run netChop on entire
  ## protein sequences
  epitopePredictionsAll <- plyr::llply(1:nrow(variantInfo), function(i) {
  # epitopePredictions <- lapply(1:nrow(variantInfo), function(i) {
    maartenutils::mymessage(
      sprintf('Processing peptides related to variant %d/%d', 
        i, nrow(variantInfo)), 
      instance = 'epitopePredictions')

    peptideList <- buildPeptideList(sequences = variantInfo[i, ],
      peptidelength = runParameters$peptidelength, 
      runParameters = runParameters)

    ## If no tumor peptides found, stop early
    if (is.null(peptideList) || nrow(peptideList[[2]]) < 1) {
      return(NULL)
    }

    ## Do predictions for all peptides
    normalAndTumorPredictions <- plyr::llply(1:2, function(k) {
    # normalAndTumorPredictions <- lapply(1:2, function(k) {
      l_peptidestretch <- 
        ifelse(k == 1, 'peptidecontextnormal', 'peptidecontexttumor')

      dtf <- performParallelPredictions(peptides = peptideList[[k]],
        runParameters = runParameters,
        peptidestretch = variantInfo[i, get(l_peptidestretch)],
        allele = runParameters$allele,
        peptidelength = runParameters$peptidelength)

      ## Rename columns
      setnames(x = dtf,
        old = c('peptide', paste0(runParameters$allele, 'affinity'),
          paste0(runParameters$allele, 'percentile_rank'), 'c_term_aa',
          'processing_score'),
        new = c('normal_peptide',
          paste0('normal_', runParameters$allele, 'affinity'),
          paste0('normal_', runParameters$allele, 'percentile_rank'),
          'normal_c_term_aa', 'normal_processing_score'))
      if (k == 2) {
        setnames(dtf, gsub('normal', 'tumor', colnames(dtf)))
      }
      return(dtf)
    })

    ## Merge tumor and normal predictions
    if (nrow(normalAndTumorPredictions[[2]]) > 0) {
      mergedPredictions <- merge(
        normalAndTumorPredictions[[1]],
        normalAndTumorPredictions[[2]],
        by = setdiff(c('c_term_pos', 'xmer', 'hla_allele', names(variantInfo)), 
          c('peptidecontextnormal', 'peptidecontexttumor')),
        all.y = TRUE, all.x = F)
    } else {
      mergedPredictions <- NULL
    }
    return(mergedPredictions)
  }, .parallel = (runParameters$ncores > 1)) %>% rbindlist(use.names = TRUE)

  if (maartenutils::null_dat(epitopePredictionsAll)) {
    return(emptyTableWithColumnNamesAndColumnClasses(
        colnames = columnNamesEmptyTable, colclasses = columnClassesEmptyTable))
  }

  ## Determine which variants contributed to the formation of predicted epitopes
  if (all(c('variant_id', 'aa_pos_germline', 'aa_pos_tumor_start',
        'aa_pos_tumor_stop', 'variant_classification') %in%
      names(variantInput))) {
    allContributingVariantsInfo <- 
      findVariantsContributingToEpitope(
        predicted_variants = epitopePredictionsAll, 
        all_variants = variantInput,
        runParameters = runParameters) %>% rbindlist(use.names = TRUE)

    epitopePredictionsAll[, contributing_variants :=
      allContributingVariantsInfo[, contributing_variants]]
    epitopePredictionsAll[, aa_peptide_pos_tumor :=
      allContributingVariantsInfo[, contributing_aa_pos_tumor]]

    if ('rna_alt_expression' %in% names(variantInput)) {
      epitopePredictionsAll[, contributing_variants_alt_expression :=
        allContributingVariantsInfo[, contributing_variants_alt_expression]]
    }

    ## If contributing variants are determined, remove variant_id &
    ## variant_classification columns (not relevant anymore)
    columnsToRetain <- setdiff(colnames(epitopePredictionsAll),
      c('variant_id', 'chromosome', 'start_position', 'end_position',
        'variant_strand', 'ref_allele', 'alt_allele', 'variant_classification',
        'codon_germline', 'codon_tumor', 'aa_germline', 'aa_tumor',
        'aa_pos_germline', 'aa_pos_tumor_start', 'aa_pos_tumor_stop'))

    epitopePredictionsAll <- 
      epitopePredictionsAll[, columnsToRetain, with = FALSE]
  }

  ## Model epitope cell surface presentation using random forest model trained
  ## on mass spec data
  if (runParameters$use_rfModel &&
    is.numeric(epitopePredictionsAll$rna_expression)) {

    ## apply random forest model to all predictions
    model_rna <- get(load(runParameters$randomForestModelPath))

    setnames(epitopePredictionsAll,
      old = c(paste0('tumor_', runParameters$allele, 'percentile_rank'),
        'tumor_processing_score', 'rna_expression'),
      new = c('affMIN', 'chop', 'rna'))

    rfModelPrediction <- predict(model_rna, epitopePredictionsAll,
      type = 'prob', na.action = 'na.pass')[, 'yes']

    if (nrow(epitopePredictionsAll) != length(rfModelPrediction))
      stop('Random forest output length unequal to input')

    epitopePredictionsAll[, model_prediction := rfModelPrediction]

    setnames(x = epitopePredictionsAll,
      old = c('affMIN', 'chop', 'rna'),
      new = c(paste0('tumor_', runParameters$allele, 'percentile_rank'),
        'tumor_processing_score', 'rna_expression'))
  }

  col_names <- c('hugo_symbol',
    setdiff(names(variantInfo), c('hugo_symbol', 'rna_expression',
        'peptidecontextnormal', 'peptidecontexttumor')),
    'rna_expression', 'c_term_pos', 'hla_allele', 'xmer',
    'tumor_peptide', 'tumor_c_term_aa',
    paste0('tumor_', runParameters$allele, 'affinity'),
    paste0('tumor_', runParameters$allele, 'percentile_rank'),
    'tumor_processing_score', 'normal_peptide', 'normal_c_term_aa',
    paste0('normal_', runParameters$allele, 'affinity'),
    paste0('normal_', runParameters$allele, 'percentile_rank'),
    'normal_processing_score', 'contributing_variants',
    'aa_peptide_pos_tumor', 'contributing_variants_alt_expression',
    'model_prediction', 'different_from_self')

  setcolorder(epitopePredictionsAll,
    intersect(col_names, colnames(epitopePredictionsAll)))

  if (runParameters$verbose) {
    maartenutils::mymessage('Writing predictions to disk')
  }

  ## write predictions to disk
  writePredictionsToDisk(dtf = epitopePredictionsAll,
    filepath = runParameters$filepath,
    filename = runParameters$filename_no_ext,
    allele = runParameters$allele,
    unique_cols = unique_cols,
    peptidelength = runParameters$peptidelength)
}
