# hydi: calling differential hydroxymethylation

## Description

hydi implements methods for calling differential hydroxymethylation from oxWGBS and oxRRBS experiments.

## Background

Briefly, ox[RR|WG]BS experiments involve sequencing the DNA of each biological sample twice - following two distinct chemical treatments. 
The first treatment with sodium bisulfite (BS) is carried out to measure the sum of 5-methyl-cytosine (5mc) and 5-hydroxymethyl-cytosine (5hmC). 
Thus, we denote the sum of both modifications by 5modC and `5modC=5mC+5hmC`.
The second treatment involves an oxidation reaction prior to the sodium bisulfite conversion (oxBS). Sequencing yields information on 5mC-levels. 
To obtain 5hmC estimates, the signals from BS and oxBS experiments are subtracted. Specifically, in hydi, we model `5hmC=5modC-5mC` employing
binomial distributions directly based on the read counts of each experiment.



## Usage

## Installing

To install hydi, please make sure to have the following dependencies installed on your system:

- zlib data compression library (zlib)
- GNU scientific library (gsl)
- GNU Multiple Precision Arithmetic Library (gmp)

Subsequently, run

```{sh}
make
```

After compilation, hit 

```{sh}
./hydi.x --help
```

to ensure the program compiled correctly.

## Input data

hydi calls differentially hydroxymethylated cytosines (d5hmC) by comparing two groups of samples analysed with
oxRRBS or oxWGBS sequencing experiments. For each group, hydi needs to be given a separate (gzipped) tab-separated file. Each file summarizes 
the count data obtained from the individual sequencing runs in a specific order. 

The first line is a header line with column descriptions and sample identifiers. 
Each of the following lines represents a single CpG in the genome. The first three fields of each line specifiy the coordinates of the CpG.

| chrom      | pos      | strand | 
| :-----------|:--------| :------|
| Chromosome (alphanumeric) | Position of CpG (integer) | strand of CpG ("+" or "-") |   

All following fields contain the count data. **For each sample, there are four fields**, i.e. two for each sequencing run:

|  ID of BS-seq (coverage)                           | ID of BS-Seq (non-conversion)       | ID of oxBS-seq (coverage)  | ID of BS-Seq (non-conversion) |
| :---------------------------------------------------| :------------------------------------| :---------------------------| :----------------------------- |
| number of reads aligned to coordinate in BS-Seq | number of non-converted Cs at coordinate in BS-Seq | number of reads aligned to coordinate in oxBS-Seq | number of non-converted Cs at coordinate in oxBS-Seq |   

where BS and oxBS denote the two distinct treatments briefly described above. For a two samples, 

```
chrom	pos	strand	SRR2074675	SRR2074675	SRR2074676	SRR2074676	SRR2074679	SRR2074679	SRR2074680	SRR2074680
chr1	608564	+	40	37	23	17	31	28	19	12
```

Running examples are provided with the code.

## Output data

Hydi returns the results in a tab-separated text file. Each line represents a single genomic CpGs and summarizes the test data for group 1 (G1) and group 2 (G2). The
frist three fields hold their respective coordinates (see input data). The following fields 18 fields are
described in the following table

|field      | name  | description |
|:----------|-------|-------|
| 4         | cil1  | lower bound of confidence interval of the 5hmC-level (G1) |
| 5         | ml1   | maximum likelihood estimate of 5hmC-level (G1) |
| 6         | ciu1  | upper bound of confidence interval of the 5hmC-level (G1) |
| 7         | pval1 | p-value for test of absence of hydroxymetylation/overshoots (G1) |
| 8         | fdr1  | fdr corrected p-value for test of absence (G1) |
|           |       |                 |
| 9         | cil2  | lower bound of confidence interval of the 5hmC-level (G2) |
| 10        | ml2   | maximum likelihood estimate of 5hmC-level in (G2) |
| 11        | ciu2  | upper bound of confidence interval of the 5hmC-level (G2) |
| 12        | pval2 | p-value for test of absence of hydroxymetylation/overshoots (G2) |
| 13        | fdr2  | fdr corrected p-value for test of absence (G2) |
|           |       |
| 14        | overshoot   | flag (**0**: no overshoot; **1**:overshoot G1; **2**:overshoot G2; **3**: overshoot G1 & G2) |
| 15        | 5hmC        | flag (**0**: no 5hmC; **1**:5hmC in G1; **2**:5hmC in G2; **3**: 5hmC in G1 & G2) |
|           |             |
| 16        | cil_diff | lower bound of confidence interval of 5hmC differences between G1 and G2|     
| 17        | ml_diff |maximum likelihood estimate of 5hmC differences between G1 and G2 |     
| 18        | ciu_diff | upper bound likelihood estimate of difference of 5hmC-levels between G1 and G2 |    
| 19        | pval_diff | p-value for test on equality of hydroxymethylation in G1 and G2 |
| 20        | fdr_diff | fdr corrected p-value for test on equality in G1 and G2|
| 21        | est_mindiff | estimated minimum difference of hydroxymethylation between G1 and G2 | 


An example output line looks like this:

```
chr1	434286	+	-0.069985	0.013117	0.096011	0.752027	0.977618	-0.087720	0.015558	0.117299	0.763707	0.951594	0	0	-0.130453	-0.002441	0.125570	0.9
70594	1.000000	0.000000
```



## Example

```
./hydi.x -a examples/G1.txt.gz -b examples/G2.txt.gz > examples/test.out
```
