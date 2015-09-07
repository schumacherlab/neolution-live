#### **Foreign Antigenic Space (FAS) database prediction pipeline**  
*Pipeline for the prediction of neo-antigens from (Alexandrov) variant calls through FASdb lookups*

---

Main script (queryDatabase.R) takes takes an input file (specified as a commandline argument) and queries contents of the Schumi FAS database with the entries it parses.  
**IMPORTANT:** For additional information regarding drawbacks & the required commandline arguments, read the **Extended information** segment below.

**Usage example:**  
`Rscript ~/fasdb/queryDatabase.R /input_data/AML.csv LAML 123`

The preceding call can be run in the Terminal and will start neo-antigen predictions from line 124 for AML.csv, using RNAlikelihood scores specified in the LAML column of the RNAlikelihood file.  

Important to note is that a MySQL server should be up-and-running before starting the script (see extended info for tables).  

**Required commandline arguments:**  

1. input file path + filename (e.g. /input_data/AML.csv)  
2. tissue type (e.g. LAML)  
3. lookup progress (e.g. 123)  

---
#### **Extended information:**  

**Important drawbacks of this version of the prediction pipeline are:**  

Context generation is performed *on-the-fly*, by querying the contents of the FAS database. This has the following limitations:  

1. in-dels are not taken into account, as nucleotide information for transcripts is not present in the database
2. the contents of the FAS database are generated from Ensembl version 58 (2010) and are therefore, arguably, outdated

**Explanation of commandline arguments:**

@1 **input file** should be a csv file containing the variants, containing the following columns:

|SAMPLE_ID|ChromosomeID|ChromosomeLocation|MutantCode|
|---------|-------------|------------------|----------|
|  .....  |    .....    |       .....      |  ......  |

@2 **tissue type** depends on RNA expression data used; for RNAlikelihood score the following tissue type codes can be used:

| Cancer type | Tissue type code | Cancer type | Tissue type code |
|-----------------------------|----|-------------|------|
|ALL (acute lymphoid leukemia)|DLBC|Lung Adeno|LUAD|
|AML (acute myeloid leukemia)|LAML|Lung Small Cell|LUSMC|
|Bladder|BLCA|Lung Squamous Cell|LUSC|
|Breast|BRCA|Lymphoma B-cell|DLBC|
|Cervix|DLBC|Medulloblastoma|GBM|
|CLL (chronic lymphoid leukemia)|CESC|Melanoma|SKCM|
|Colorectum|COAD|Myeloma|DLBC|
|Esophageal|STAD|Neuroblastoma|ACC|
|Glioblastoma|GBM|Ovary|OV|
|Glioma (Low grade)|LGG|Pancreas|PAAD|
|Head and Neck|HNSC|Pilocytic Astrocytoma|LGG|
|Kidney Chromophobe|KICH|Prostate|PRAD|
|Kidney Clear Cell|KIRC|Stomach|STAD|
|Kidney Papillary|KIRP|Thyroid|THCA|
|Liver|LIHC|Uterus|UCEC|

@3 **lookup progress** specifies where to start in input file, starts at passed index +1; is automatically obtained from FAS database when lookups are controlled through masterLookupController.R

**MySQL database definition:**

The MySQL database has been filled with information obtained from the group of Michael Stratton at the Sanger Institute. They generated the following data, by making all possible mutations in the coding sequence of all canonical transcripts contained in Ensembl version 58.

| Table | Contents | Headers |
|-------|----------|---------|
|codonToAAC| Conversion table of codon to amino acid residue |'codon' ; 'AA' |
|ENSGtoStrand| Transcript orientation for each ENSG# | 'ENSG' ; 'Strand' |
|idInfo| For every ENSG#, canonical ENST#, CCDS# and Gene Symbol | 'geneID' ; 'ENSG' ; 'CCDS' ; 'ENST' |
|lookupProgress| Index of last processed line per dataset | 'dataset' ; 'progress' ; 'total' |
|perChromosomeLocationENSG| Genomic coordinates and transcript-/codon-/protein positions for all ENSG# | 'chromosomeID' ; 'chromosomeLocation' ; 'ENSG' ; 'transcriptPosition' ; 'codonPosition' ; 'proteinPosition' |
|perGenePepProcessingScores| Proteasomal processing scores for all peptides per ENSG# | 'ENSG' ; 'peptide' ; 'peptideStart' ; 'processingScore' |
|perPepAffinityScores_NEW| MHC binding affinity scores for all generated peptides | 'peptide' ; 'A0101affinity' ; 'A0201affinity' ; 'A0301affinity' ; 'B0702affinity' ; 'B0801affinity' |
|perProteinAAC| Wildtype codon and amino acid at all protein positions per ENSG# | 'ENSG' ; 'proteinPosition' ; 'codon' ; 'AA' |
|pmbecScores| PMBEC scores for all possible amino acid substitutions | 'AA1' ; 'AA2' ; 'PMBECscore' |