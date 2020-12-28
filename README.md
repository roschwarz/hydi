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

## Parameter

```
usage: ./hydi.x [-a <file>] [-b <file>] [-o <file>] [-m] [-p] [-e]
  Testing for the equality of hydroxymethylation from oxbs sequencing

 [INPUT]
 -a, --group1 <file>  path/filename of count file for group 1 (default:none)
 -b, --group2 <file>  path/filename of count file for group 2 (default:none)
 [OUTPUT]
 -o, --output <file>  path/filename of output (default:none)
 [CONTROL]
 -m, --maxiter        maximum iterations in gradient descent (default:100)
 -p, --alpha          significance level for flags and confidence intervals (default:0.050000)
 -e, --epsilon        precision of numerical approximation (default:1.000000e-06)
```

### Input and output files

| parameter  | description      |
|----------- |-------------------------------------|
| group1     | formated input file for the 1st group (G1). The format is described below.                      | 
| group2     | formated input file for the 2nd group (G2). The format is described below.                      | 
| output     | path or filename for the output. If option is omitted, the output will be dumped to standard output  |

Please note, the assignment of the sample files to the groups matter. The calculation of estimates and confidence intervals for hydroxymethylation uses G1 as a reference. Thus, negative values
indicate a loss of hydroxymethylation in G1 as compared to G2. Conversely, positive values indicate a gain of hydroxymethylation G1. In a comparison of a group of cancer samples versus healthy controls, for instance, it is thus
recommendable to assign the cancer data to G1 to make the data easier to interpret.

### Control

| parameter  | description      |
|----------- |-------------------------------------|
| maxiter    | control of the number of iterations for the gradient descent used in the estimation of maximum likelihoods  | 
| alpha      | parameter controlling the level of confidence intervals and used for flags | 
| epsilon    | parameter controlling the numerical precision of all calculations | 

Please note, while increasing maxiter and decreasing epsilon may provide a better precision, the runtime of hydi may be affected substantially. Convenience flags of hydis output are set based on the value of `alpha`. Furthermore, hydi calculates `1-alpha` confidence intervals.


## Input data format

hydi calls differentially hydroxymethylated cytosines (d5hmC) by comparing two groups of samples analysed with
oxRRBS or oxWGBS sequencing experiments. For each group, hydi needs to be given a separate (gzipped) tab-separated file. Each file summarizes 
the count data obtained from the individual sequencing runs in a specific order. 

The first line is a header line with column descriptions and sample identifiers. 
Each of the following lines represents a single C in the genome. The first three fields of each line specify the coordinates of the C.


All following fields contain the count data. **For each sample, there are four fields**, i.e. two for each sequencing run:

|  column  | description |
| :---------------------------------------------------| :------------------------------------| 
| chrom                            | chromosome (alphanumeric) | 
| pos                              | position of C (integer) | 
| strand                           | strand of C ("+" or "-") |   
|                                  |                            |
| ID of BS-seq (coverage)          |  number of reads aligned to coordinate in BS-Seq | 
| ID of BS-Seq (non-conversion)    | number of non-converted Cs at coordinate in BS-Seq | 
| ID of oxBS-seq (coverage)        |  number of reads aligned to coordinate in oxBS-Seq | 
| ID of BS-Seq (non-conversion)    | number of non-converted Cs at coordinate in oxBS-Seq |   

where BS and oxBS denote the two distinct treatments briefly described above. For a two samples, 

```
chrom	pos	strand	SRR2074675	SRR2074675	SRR2074676	SRR2074676	SRR2074679	SRR2074679	SRR2074680	SRR2074680
chr1	608564	+	40	37	23	17	31	28	19	12
```

Running examples are provided with the code.

## Output data

Hydi returns the results in a tab-separated text file. Each line represents a single genomic Cs and summarizes the test data for group 1 (G1) and group 2 (G2). The
first three fields hold their respective coordinates (see input data). The following fields 18 fields are
described in the following table

|field      | name  | description |
|:----------|:-------|:-------|
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


All confidence interval bounds and estimates are given as rates, i.e. on the interval of [0,1]. Confidence intervals are calculated to the `1-alpha`-level. 
The convenience flag `overshoot` is set for one or both groups if the maximum likelihood estimate for 5hmC is negative and significantly different from 0 (based on the chosen `alpha`).
Similarly, the `5hmC`-flag is set for one or both groups if the estimates are positive and significantly different from 0 (base on the chosen `alpha`). 

An example output line looks like this:

```
chr1	434286	+	-0.069985	0.013117	0.096011	0.752027	0.977618	-0.087720	0.015558	0.117299	0.763707	0.951594	0	0	-0.130453	-0.002441	0.125570	0.9
70594	1.000000	0.000000
```

## Recommendations

### Analysis of differential 5hmC

For the analysis of differential hydroxymethylation it is recommended to filter out all Cs with overshoots, i.e. a statistically significant
**negative** methylation. This effect may be caused by coverage or alignment problems. To do this, the convenience flag `overshoot` (field 14) may be
used. This can be done, for instance, by calling

```{sh}
awk '{if(NR == 1 || $14 == 0) print }' examples/output.txt > examples/output.noovershoot.txt
```

For all downstream analyses it is recommended to filter the data for a desired fdr cutoff, e.g. fdr <= 0.1, using `fdr_diff` (field 20) and to use the convenience field `est_mindiff` (field 21) to eliminate biologically irrelevant differences. Both values are calculated based on the difference p-value (field 19) and the confidence interval (field 16-18), respectively. For instance, the shell command

```{sh}
awk '{if(NR == 1 || ($14 == 0 && $20 <= 0.1 && $21 >= 0.1)) print }' examples/output.txt > examples/output.noovershoot.filtered.txt
```
extracts all sites without significant overshoots, significant fdr-corrected differential hydroxymethylation and a minimum difference of at least 0.1.   

### Analysis of 5hmC 

If hydi is used for identifying 5hmCs in either group, the convenience flag `5hmC` (field 15) may come in handy. For instance, to identify
sites that are hydroxymethylated in both groups run

```{sh}
awk '{if(NR == 1 || $15 == 3) print }' examples/output.txt > examples/output.5hmC.txt
```

## Example

```
./hydi.x -a examples/G1.txt.gz -b examples/G2.txt.gz > examples/test.out
```

### Complaint department

steve hoffmann leibniz minus fli de
robert schwarz leibniz minus fli de

## Tutorial

The following tutorial shows the complete workflow for dhmC analysis from (ox)BS-Seq fastq data to calling differenital hydroxymethylation with hydi.
