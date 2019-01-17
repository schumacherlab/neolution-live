# Neolution-live
## Pipeline for neo-antigen prediction

### Introduction

This pipeline performs neo-antigen predictions using netMHCpan and netChop. RNA expression 
data and 'similarity-to-self' filters can additionally be used to further increase the 
precision of the predictions. Alternatively, we trained a random forest classifier on mass 
spectrometry data, which integrates predicted affinity rank, proteasomal processing scores 
and RNA expression values to a combined model score.  This score lies in the range [0-1] 
and gives the probability of peptide presentation and has been shown to substantially 
increase the precision and sensitivity over the conventional use of binary cutoff values 
for the various predicted parameters.  This is now our preferred prediction method.  
  Suggested model score cutoff values are:

* 0.01 for TIL screens (more inclusive; we picked up low magnitude TIL hits at low 
  probability scores)
* 0.02 for PBMC screens (more stringent; unlikely to pick up low magnitude responses with low prob scores in peripheral blood)

The pipeline can be supplied with different forms of input. Either As input, the pipeline 
expects a tsv file with affected germline and tumor transcripts. Additional 
variant/transcript information can be provided and will be transferred into the output. 
See **Input file format** paragraph for more information regarding input file generation.

### Minimal usage example

The only required argument is the input file, other settings will will take default values if these are not supplied

`neolution --mhc A0201 --model 0.02 --length 9 /path/to/variants.tsv`

The call will start neo-antigen predictions for **variants.tsv**, __HLA-A*02:01__, 
**9-mer** peptides, applying a model prediction cutoff score of 0.02.  

`neolution /path/to/variants.tsv`

Run Neolution using default settings (see below)

##

### Commandline arguments

#### Required

1. full input file path

#### Optional

1. netChop processing cutoff
2. RNA expression cutoff
3. single sequence input (fasta input: not paired tumor-normal, no rna expression)
4. structural variant predictions
5. simple self-similarity check (9-, 10-, 11-mers)
6. extended self-similarity check (9-mers only)
7. use self-epitope list
8. use database for peptide affinity lookups (9-mers, netMHCpan-2.4 only)
9. netMHCpan version

**NOTE: self-similarity checking requires predicted self-epitope lists of matching HLA & peptide length**

### Input file format

Required input format is a wide, tab-separated table with the following columns. 
[Neolution-prep](https://gitlab.nki.nl/l.fanchi/neolution-prep) generates all required 
data and provides required input file.

#### Variant

| variant\_id | chromosome | start\_position | end\_position | variant\_strand | ref\_allele | alt\_allele |
|-------------|------------|-----------------|---------------|-----------------|-------------|-------------|

#### Gene/Transcript

| gene\_id | transcript\_id | transcript\_strand | hugo\_symbol | variant\_classification | transcript\_remark | transcript\_extension | nmd\_status | nmd\_remark |
|----------|----------|---------|----------|-----------|---------|------------|-------------|-------------|

#### DNA/RNA

| dna\_ref\_read\_count | dna\_alt\_read\_count | dna\_total\_read\_count | dna\_vaf | rna\_ref\_read\_count | rna\_alt\_read\_count | rna\_total\_read\_count | rna\_vaf | rna\_alt\_expression | rna\_expression |
|----------|-----------|------------|------------|----------|-----------|-----------|------------|--------|--------|

#### Effect

| codon\_germline | codon\_tumor | aa\_germline | aa\_tumor | aa\_pos\_germline | aa\_pos\_tumor\_start | aa\_pos\_tumor\_stop | peptidecontextnormal | peptidecontexttumor |
|----------|-----------|------------|-------------|-----------|------------|-----------|----------|--------|

---

`Rscript performPredictions.R --help`  

**Usage: performPredictions.R [OPTIONS]**

**OPTIONS**  
`-f FILE, --file=FILE`  
*Full path to file containing variant calls (required)*

`-m MHC, --mhc=MHC`  
*MHC/HLA allele, formatted as follows: A0201 (required)*

`-l LENGTH, --length=LENGTH`  
*Peptide length (required)*

`-d MODEL, --model=MODEL`  
*Random forest model score cutoff (optional, suggested cutoff 0.02 for PBMC, 0.01 for TIL screens)*

`-r RANK, --rank=RANK`  
*netMHCpan rank cutoff (optional, suggested is 2.8 ~500nM for A\*02:01)*

`-a AFFINITY, --affinity=AFFINITY`  
*netMHCpan affinity cutoff (optional)*

`-p PROCESSING, --processing=PROCESSING`  
*netChop processing score cutoff (optional, default: >= 0.5)*

`-e EXPRESSION, --expression=EXPRESSION`  
*RNA expression cutoff (optional, default: > 0; use -1 for no filtering)*

`--single`  
*Single sequence predictions (not paired normal-tumor) (optional, default: FALSE)*

`--structural`  
*Structural variant predictions (optional, default: FALSE)*

`--selfsim`  
*Perform simple self-similarity check; compatible with 9-, 10-, 11-mers (optional, default: FALSE)*

`--extselfsim`  
*Perform extended self-similarity check; only compatible with 9-mers (optional, default: FALSE)*

`--selflist`  
*Add predicted self-epitopes to self-similarity check, requires length- & HLA-matched selflist (optional, default: FALSE)*

`--fasdb`  
*Look up peptide affinity in FASdb, predict if not found; only compatible with 9-mers & netMHCpan-2.4 (optional, default: FALSE)*

`--panversion`  
*Use different version of netMHCpan; must be installed in path specified in runConfig.R (optional, default: 4.0)*

`--verbose`  
*Be chatty (optional, default: FALSE)*

`-h, --help`  
*Show this help message and exit*
