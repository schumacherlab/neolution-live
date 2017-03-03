#### Lookup parameters and other configuration options can be set here ####

runOptions = list(general = list(numberOfWorkers = detectCores()/4, # set number of cores available for parallel processing
																 temporaryDirectoryPath = '~/scratch'), # set path of temporary directory (set chmod 1777 on tempdir)

									predictors = list(netMHCpan = paste0('/home/NFS/users/l.fanchi/predictors/netMHCpan-', runParameters$panversion, '/netMHCpan'),
																		netChop = '/home/NFS/users/l.fanchi/predictors/netchop-3.1/bin/netChop'),

									resources = list(selfEpitopeListPath = './selflists', # set path of self-epitope lists
																	 randomForestModelPath = './resources/rf_model.RData') # set path to random forest model
									)
