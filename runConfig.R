#### Lookup parameters and other configuration options can be set here ####

# get number of cores available for parallel processing
numberOfWorkers=5

# set min. affinity (in nM), min. chop score & min. RNA expression
affinityLimit=500
chopLimit=0.5
expressionLimit=0

# set peptide length to be analysed (only 9mers supported at the moment)
xMer=9

# set hla types to be analysed and location of file containing RNA expression data
hlaTypes=c("A0101","A0201","A0301","B0702","B0801")
rnaExpressionFile="20150723_RNAexpression_likelihood_TCGA_IlluminaHiSeq.csv"

# set predictor paths
netMHCpath="/home/NKI/l.fanchi/netMHC-3.0/netMHC-3.0"
netMHCpanpath="/home/NKI/l.fanchi/netMHCpan-2.4/netMHCpan"
netChoppath="/home/NKI/l.fanchi/netchop-3.1/bin/netChop"

# set mysql configuration
sqlhost = "medoid"
sqluser = "l.fanchi"
sqlpass = "MpRi1RKd"
sqldbname = "SchumiDB"

sqlprogressdb = "otacon"