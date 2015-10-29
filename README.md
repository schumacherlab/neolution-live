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

|Allele|Population<br>frequency|Affinity<br>cutoff|
|------|:----:|:---:|
|A*0101|16.2|884|
|A*0201|25.2|255|
|A*0203|3.3|92|
|A*0206|4.9|60|
|A*0301|15.4|602|
|A*1101|12.9|382|
|A*2301|6.4|740|
|A*2402|16.8|849|
|A*2501|2.5|795|
|A*2601|4.7|815|
|A*2902|2.9|641|
|A*3001|5.1|109|
|A*3002|5|674|
|A*3101|4.7|329|
|A*3201|5.7|131|
|A*3301|3.2|606|
|A*6801|4.6|197|
|A*6802|3.3|259|
|B*0702|13.3|687|
|B*0801|11.5|663|
|B*1402|2.8|700|
|B*1501|5.2|528|
|B*1801|4.4|732|
|B*2705|2|584|
|B*3501|6.5|348|
|B*3503|1.2|888|
|B*3801|2|944|
|B*3901|2.9|542|
|B*4001|10.3|639|
|B*4002|3.5|590|
|B*4402|9.2|904|
|B*4403|7.6|780|
|B*4601|4|926|
|B*4801|1.8|887|
|B*5101|5.5|939|
|B*5301|5.4|538|
|B*5701|3.2|716|
|B*5801|3.6|446|

[1 http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions](http://help.iedb.org/entries/23854373-Selecting-thresholds-cut-offs-for-MHC-class-I-and-II-binding-predictions)