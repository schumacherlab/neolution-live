### **Neolution pipeline**  
*Pipeline for the prediction of neo-antigens*

---

Two branches are currently under development:
1. fasdb-based - uses csv input and epitopes in fasdb for initial lookups (any missing data is generated on-the-fly)
2. live - uses vcf input and performs all predictions live (to be used e.g. for patient predictions)

**IMPORTANT:** For additional information regarding the required commandline arguments, read the **Extended information** segment below.

**Usage example:**  
`Rscript ~/fasdb/queryDatabase.R /input_data/AML.csv LAML 123`

The preceding call can be run in the Terminal and will start neo-antigen predictions from line 124 for AML.csv, using RNAlikelihood scores specified in the LAML column of the RNAlikelihood file.  

Important to note is that a MySQL server should be up-and-running before starting the script (see extended info for tables).  

**Required commandline arguments:**  

1. input file path + filename (e.g. /input_data/AML.csv)  
2. tissue type (e.g. LAML)  
3. lookup progress (e.g. 123)  

---

#### **Explanation of commandline arguments:**

@1 **input file** should be a vcf file containing the variants

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

@3 **lookup progress** specifies where to start in input file, starts at passed index +1; is automatically obtained from MySQL server when lookups are controlled through masterLookupController.R

---

#### **MySQL database definition:**

The MySQL database has been filled with information obtained from the group of Michael Stratton at the Sanger Institute. They generated the following data by making all possible mutations in the coding sequence of all canonical transcripts contained in Ensembl version 58.

| Table | Contents | Headers |
|-------|----------|---------|
|lookupProgress| Index of last processed line per dataset | 'dataset' ; 'progress' ; 'total' |
|perPepAffinityScores_NEW| MHC binding affinity scores for all generated peptides | 'peptide' ; 'A0101affinity' ; 'A0201affinity' ; 'A0301affinity' ; 'B0702affinity' ; 'B0801affinity' |