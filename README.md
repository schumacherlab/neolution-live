### **Neolution-live pipeline**  
*Pipeline for the prediction of neo-antigens*

---

This pipeline performs live predictions using netMHCpan and netChop to predict peptide affinity and proteasomal processing. RNA expression data and 'similarity-to-self' filters can used to further increase the accuracy of the predictions.

**Usage example:**  
`Rscript performPredictions.R -f /home/NFS/users/username/patient/variants.tsv -m A0201 -l 9`

The call should be run from the script directory from the Terminal and will start neo-antigen predictions for **variants.tsv**, __HLA-A*02:01__ and **9-mer** peptides. 

**By default, netMHCpan v3.0 will be used and 1/4th of the available cores are used per run for parallel computations; four samples can be processed in parallel on one cluster.  
Make sure you *nice* your runs and don't exceed max. HPC load (max. load = # cores)!!**

**IMPORTANT:** For additional information regarding the required commandline arguments, read segment below.

**Required commandline arguments:**  

1. full input file path
2. hla/mhc type (e.g. A0201)
3. peptide length (e.g. 9) 

**Optional commandline arguments:**  

1. netMHCpan affinity cutoff
2. netMHCpan rank cutoff
2. netChop processing cutoff
3. rna expression cutoff
4. single sequence input (fasta input: not paired tumor-normal, no rna expression)
5. simple self-similarity check (9-, 10-, 11-mers)
6. extended self-similarity check (9-mers only)
7. use self-epitope list
8. use database for peptide affinity lookups (9-mers, netMHCpan-2.4 only)
9. netMHCpan version

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

`-r RANK, --rank=RANK`  
*netMHCpan rank cutoff (optional, default: FALSE)*

`-p PROCESSING, --processing=PROCESSING`  
*netChop processing score cutoff (optional, default: >= 0.5)*

`-e EXPRESSION, --expression=EXPRESSION`  
*RNA expression cutoff (optional, default: > 0; use -1 for no filtering)*

`--single`  
*Single sequence predictions (not paired normal-tumor) (optional, default: FALSE)*

`--selfsim`  
*Perform simple self-similarity check; compatible with 9-, 10-, 11-mers (optional, default: FALSE)*

`--extselfsim`  
*Perform extended self-similarity check; only compatible with 9-mers (optional, default: FALSE)*

`--selflist`  
*Add predicted self-epitopes to self-similarity check, requires length- & HLA-matched selflist (optional, default: FALSE)*

`--fasdb`  
*Look up peptide affinity in FASdb, predict if not found; only compatible with 9-mers & netMHCpan-2.4 (optional, default: FALSE)*

`--panversion`  
*Use different version of netMHCpan; must be installed in path specified in runConfig.R (optional, default: 3.0)*

`-h, --help`  
*Show this help message and exit*

---

### Suggested allele-specific cutoffs

|Allele|Affinity <br>cutoff (nM)|Population <br>frequency|
|:------:|:----:|:---:|
|A\*01:01|884|16.2|
|A\*02:01|255|25.2|
|A\*02:03|92|3.3|
|A\*02:06|60|4.9|
|A\*03:01|602|15.4|
|A\*11:01|382|12.9|
|A\*23:01|740|6.4|
|A\*24:02|849|16.8|
|A\*25:01|795|2.5|
|A\*26:01|815|4.7|
|A\*29:02|641|2.9|
|A\*30:01|109|5.1|
|A\*30:02|674|5|
|A\*31:01|329|4.7|
|A\*32:01|131|5.7|
|A\*33:01|606|3.2|
|A\*68:01|197|4.6|
|A\*68:02|259|3.3|
|B\*07:02|687|13.3|
|B\*08:01|663|11.5|
|B\*14:02|700|2.8|
|B\*15:01|528|5.2|
|B\*18:01|732|4.4|
|B\*27:05|584|2|
|B\*35:01|348|6.5|
|B\*35:03|888|1.2|
|B\*38:01|944|2|
|B\*39:01|542|2.9|
|B\*40:01|639|10.3|
|B\*40:02|590|3.5|
|B\*44:02|904|9.2|
|B\*44:03|780|7.6|
|B\*46:01|926|4|
|B\*48:01|887|1.8|
|B\*51:01|939|5.5|
|B\*53:01|538|5.4|
|B\*57:01|716|3.2|
|B\*58:01|446|3.6|


[1 http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions](http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions)