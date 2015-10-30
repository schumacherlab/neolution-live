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

---

### Suggested allele-specific cutoffs

|Allele|Population <br>frequency|Affinity <br>cutoff|
|:------:|:----:|:---:|
|A\*01:01|16.2|884|
|A\*02:01|25.2|255|
|A\*02:03|3.3|92|
|A\*02:06|4.9|60|
|A\*03:01|15.4|602|
|A\*11:01|12.9|382|
|A\*23:01|6.4|740|
|A\*24:02|16.8|849|
|A\*25:01|2.5|795|
|A\*26:01|4.7|815|
|A\*29:02|2.9|641|
|A\*30:01|5.1|109|
|A\*30:02|5|674|
|A\*31:01|4.7|329|
|A\*32:01|5.7|131|
|A\*33:01|3.2|606|
|A\*68:01|4.6|197|
|A\*68:02|3.3|259|
|B\*07:02|13.3|687|
|B\*08:01|11.5|663|
|B\*14:02|2.8|700|
|B\*15:01|5.2|528|
|B\*18:01|4.4|732|
|B\*27:05|2|584|
|B\*35:01|6.5|348|
|B\*35:03|1.2|888|
|B\*38:01|2|944|
|B\*39:01|2.9|542|
|B\*40:01|10.3|639|
|B\*40:02|3.5|590|
|B\*44:02|9.2|904|
|B\*44:03|7.6|780|
|B\*46:01|4|926|
|B\*48:01|1.8|887|
|B\*51:01|5.5|939|
|B\*53:01|5.4|538|
|B\*57:01|3.2|716|
|B\*58:01|3.6|446|

[1 http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions](http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions)