#### Lookup parameters and other configuration options can be set here ####
data_path = '/DATA/users/l.fanchi'

runOptions = list(general = list(numberOfWorkers = round(detectCores()/6, 0), # set number of cores available for parallel processing
                                 temporaryDirectoryPath = file.path(data_path, 'scratch')), # set path of temporary directory (set chmod 1777 on tempdir)

                  predictors = list(netMHCpan = file.path(data_path, paste0('predictors/netMHCpan-', runParameters$panversion, '/netMHCpan')),
                                    netChop = file.path(data_path, 'predictors/netchop-3.1/bin/netChop')),

                  resources = list(selfEpitopeListPath = file.path(data_path, 'resources/neolution_selflists'), # set path of self-epitope lists
                                   randomForestModelPath = './resources/rf_model.RData') # set path to random forest model
)

regexPatterns = list(file_extension = '\\.[^.]+$', # match file extension (everything after last dot, inclusive)
                     snp_identifier = '[gr]s\\d+', # for matching SNPs
                     gs_identifier = 'gs\\d+', # for matching snps not found in dbSNP, keep boundless (no '^' or '$')
                     rs_identifier = 'rs\\d+', # for matching snps found in dbSNP, keep boundless (no '^' or '$')
                     cosmic_identifier = 'COSM\\d+', # for matching variants found in COSMIC coding muts database, keep boundless
                     seqdata_prefix = '(^.+[ATCG]{7,8})(.+)', # for isolating GCF prefix
                     allele_exclusion = 'C[0-9]{4}') # for excluding particular alleles from analysis
