#### Lookup parameters and other configuration options can be set here ####

# set number of cores available for parallel processing
numberOfWorkers=10

# set min. affinity (in nM), min. chop score & min. RNA expression
affinityLimit=500
chopLimit=0.5
expressionLimit=0

# set peptide length to be analysed (only 9mers supported at the moment)
peptideLength=9

# set hla types to be analysed 
hlaTypes=c("A0101",
           "A0201",
           "A0301",
           "B0702",
           "B0801")

# set location of file containing RNA expression data
rnaExpressionPath="./rna_expr_data/20150723_RNAexpression_likelihood_TCGA_IlluminaHiSeq.csv"

# set location and base filename of self-peptide lists
# self-peptide list filenames should follow convention: 'selfPeptideListBaseName'_'hlaTypes[x]'_'peptideLength'mer_epitopes.csv
selfPeptideListPath="./selflists_v2.4"
selfPeptideListBaseName="20140610_human_proteome"

# set predictor paths
predictorPaths=data.table(netMHCpan="/home/NKI/l.fanchi/netMHCpan-2.4/netMHCpan",
                          netChop="/home/NKI/l.fanchi/netchop-3.1/bin/netChop")

# set mysql configuration
sqlhost = "medoid"
sqluser = "l.fanchi"
sqlpass = "MpRi1RKd"
sqldbname = "SchumiDB"

sqlprogressuser = "otacon"
sqlprogresspass = "DSmxaoxA"
sqlprogressdbname = "otacon"