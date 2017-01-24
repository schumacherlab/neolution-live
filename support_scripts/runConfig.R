#### Lookup parameters and other configuration options can be set here ####

# set number of cores available for parallel processing
numberOfWorkers = detectCores()/4

# set predictor paths
predictorPaths = data.table(netMHCpan = paste0("/home/NFS/users/l.fanchi/predictors/netMHCpan-", runParameters$panversion, "/netMHCpan"),
                            netChop = "/home/NFS/users/l.fanchi/predictors/netchop-3.1/bin/netChop")

# set path of self-epitope lists
selfEpitopeListPath = "./selflists"

# set path of temporary directory (set chmod 1777 on tempdir)
temporaryDirectoryPath = "~/scratch"
