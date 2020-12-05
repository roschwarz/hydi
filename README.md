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

All following fields contain the count data. *For each sample, there are four fields*, i.e. two for each sequencing run:

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

Hydi returns the results in a tab-separated text file. Each line represents a single genomic CpGs and the first three fields give their coordinates

| chrom      | pos      | strand | 
| :-----------|:--------| :------|
| Chromosome (alphanumeric) | Position of CpG (integer) | strand of CpG ("+" or "-") |   


The following fields 19 fields are as follows

|field      | description |
|:----------|-------------|
| 4         | lower bound of confidence interval of the 5hmC-level in group 1 |
| 5         | maximum likelihood estimate of 5hmC-level in group 1 |
| 6         | upper bound of confidence interval of the 5hmC-level in group 1 |
| 7         | p-value for test of absence of hydroxymetylation/overshoots in group 1 |
| 8         | fdr corrected p-value for test of absence in group 1 |
|           |                   |
| 9         | lower bound of confidence interval of the 5hmC-level in group 2 |
| 10        | maximum likelihood estimate of 5hmC-level in group 2 |
| 11        | upper bound of confidence interval of the 5hmC-level in group 2 |
| 12        | p-value for test of absence of hydroxymetylation/overshoots in group 2 |
| 13        | fdr corrected p-value for test of absence in group 2 |
|           |              |
| 14        | overshoot flag (0: no overshoot; 1:overshoot in group 1; 2:overshoot in group 2; 3: overshoot in both groups) |
| 15        | 5hmC flag (0: no 5hmC; 1:5hmC in group 1; 2:5hmC in group 2; 3: 5hmC in both groups) |
|           |             |
| 16        | lower bound of confidence interval of 5hmC differences between both groups|     
| 17        | maximum likelihood estimate of 5hmC differences between both groups |     
| 18        | upper bound likelihood estimate of difference of 5hmC-levels |    
| 19        | p-value for test on equality of hydroxymethylation |
| 20        | fdr corrected p-value for test on equality |
| 21        | estimated minimum difference of hydroxymethylation | 



```
chr1	434286	+	-0.069985	0.013117	0.096011	0.752027	0.977618	-0.087720	0.015558	0.117299	0.763707	0.951594	0	0	-0.130453	-0.002441	0.125570	0.9
70594	1.000000	0.000000
```

