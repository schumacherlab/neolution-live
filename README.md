### **Neolution-live pipeline**  
*Pipeline for the prediction of neo-antigens*

---

Two branches are currently under development: 

1. fasdb-based - performs initial affinity lookups in the fasdb (any missing data is generated on-the-fly)  
2. live (this one)- performs all predictions live (to be used e.g. for patient predictions)  

**Usage example:**  
`Rscript performPredictions.R -f /home/NFS/users/l.fanchi/neolution-live/rte_kitchensink.txt -m A0201 -l 9`

The call should be run from the script directory from the Terminal and will start neo-antigen predictions for **rte_kitchensink.txt**, __HLA-A*02:01__ and **9-mer** peptides.

**IMPORTANT:** For additional information regarding the required commandline arguments, read segment below.

**Required commandline arguments:**  

1. full input file path
2. hla/mhc type (e.g. A0201)
3. peptide length (e.g. 9) 

**Optional commandline arguments:**  

1. netMHCpan affinity cutoff
2. netChop processing cutoff
3. rna expression cutoff
4. simple self-similarity check (9-, 10-, 11-mers)
5. extended self-similarity check (9-mers only)
6. single sequence input (not paired tumor-normal)
7. use self-epitope list

**NOTE: self-similarity checking requires predicted self-epitope lists of matching HLA & peptide length**

---

`Rscript performPredictions.R --help`  

**Usage: performPredictions.R [options]**

**Options:**  
`-f FILE, --file=FILE`  
*Full path to file containing variant calls (required)*

`-m MHC, --mhc=MHC`  
*MHC/HLA allele, formatted as follows: A0201 (required)*

`-l LENGTH, --length=LENGTH`  
*Peptide length (required)*

`-a AFFINITY, --affinity=AFFINITY`  
*netMHCpan affinity cutoff (optional, default: <= 500 nM)*

`-p PROCESSING, --processing=PROCESSING`  
*netChop processing score cutoff (optional, default: >= 0.5)*

`-e EXPRESSION, --expression=EXPRESSION`  
*RNA expression cutoff (optional, default: > 0)*

`--single`  
*Single sequence predictions (not paired normal-tumor) (optional, default: FALSE)*

`--selfsim`  
*Perform simple self-similarity check; compatible with 9-, 10-, 11-mers (optional, default: FALSE)*

`--extselfsim`  
*Perform extended self-similarity check; only compatible with 9-mers (optional, default: FALSE)*

`--selflist`  
*Add predicted self-epitopes to self-similarity check, requires length- & HLA-matched selflist (optional, default: FALSE)*

`-h, --help`  
*Show this help message and exit*