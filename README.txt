Foreign Antigenic Space database query script
Pipeline for the prediction of neo-antigens from variant calls

Main script (queryDatabase.R) takes takes an input file (specified as a commandline argument) and queries the Schumi FAS database with the entries it parses. For additional required commandline arguments, keep reading.

Required commandline arguments (see below for extended info):  

1. input file path + filename (e.g. /input_data/AML.csv)  
2. tissue type (e.g. LAML)  
3. lookup progress (e.g. 123)  

Extended information:

@1 input file should be a csv file containing the variants, containing the following columns:

|SAMPLE_ID|Chromosome_ID|ChromosomeLocation|MutantCode|
|---------|-------------|------------------|----------|
|  .....  |    .....    |       .....      |  ......  |

@2 tissue type depends on RNA expression data used; for RNAlikelihood score the following can be used:

|Cancer type 					|Tissue type code	|Cancer type 			|Tissue type code	|
|-------------------------------|-------------------|-----------------------|-------------------|
|ALL (acute lymphoid leukemia)	|DLBC				|Lung Adeno				|LUAD				|
|AML (acute myeloid leukemia)	|LAML				|Lung Small Cell		|LUSMC				|
|Bladder						|BLCA				|Lung Squamous Cell		|LUSC				|
|Breast							|BRCA				|Lymphoma B-cell		|DLBC				|
|Cervix							|DLBC				|Medulloblastoma		|GBM				|
|CLL (chronic lymphoid leukemia)|CESC				|Melanoma				|SKCM				|
|Colorectum						|COAD				|Myeloma				|DLBC				|
|Esophageal						|STAD				|Neuroblastoma			|ACC				|
|Glioblastoma					|GBM				|Ovary					|OV					|
|Glioma (Low grade)				|LGG				|Pancreas				|PAAD				|
|Head and Neck					|HNSC				|Pilocytic Astrocytoma	|LGG				|
|Kidney Chromophobe				|KICH				|Prostate				|PRAD				|
|Kidney Clear Cell				|KIRC				|Stomach				|STAD				|
|Kidney Papillary				|KIRP				|Thyroid				|THCA				|
|Liver							|LIHC				|Uterus					|UCEC				|

@3 lookup progress specifies where to start in input file, starts at passed index +1; is automatically obtained from FAS database when lookups are controlled through masterLookupController.R