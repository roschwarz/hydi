# hydi: calling differential hydroxymethylation

## Description

hydi implements methods for calling differential hydroxymethylation from oxWGBS and oxRRBS experiments.

## Background

Briefly, ox[RR|WG]BS experiments involve sequencing the DNA of each biological sample twice - following two distinct chemical treatments. 
The first treatment with sodium bisulfite (BS) is carried out to measure the sum of 5-methyl-cytosine (5mc) and 5-hydroxymethyl-cytosine (5hmC). 
Thus, we denote the sum of both modifications by 5modC and $5modC=5mC+5hmC$.
The second treatment involves an oxidation reaction prior to the sodium bisulfite conversion (oxBS). Sequencing yields information on 5mC-levels. 
To obtain 5hmC estimates, the signals from BS and oxBS experiments are subtracted. Specifically, in hydi, we model $5hmC=5modC-5mC$ employing
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

The first line is a header line with column descriptions and sample identifiers, e.g.

```
chrom	pos	strand	SRR2074675_N	SRR2074675_C	SRR2074676_N	SRR2074676_C
```

Each of the following lines represents a single CpG in the genome. The first three fields of each line specifiy the coordinates of the CpG.

| chrom      | pos      | strand | 
| -----------|:--------:| ------:|
| Chromosome (alphanumeric) | Position of CpG (integer) | strand of CpG (character; [+-]) |   

All following fields contain the count data. For each sample, there are four fields, i.e. two for each sequencing run:

| coverage (BS) | non-converted Cs (BS) | coverage (oxBS)  | converted Cs (oxBS) |
| --------------| ------------------------ | -----------------| -------------------------- |
| number of all reads aligned in BS-Seq | number of non-converted Cs in BS-Seq | number of all reads aligned in oxBS-Seq | number of non-converted Cs in oxBS-Seq |   

where BS and oxBS denote the two distinct treatments briefly described above. For a single sample

```
chrom	pos	strand	SRR2074675_N	SRR2074675_C	SRR2074676_N	SRR2074676_C
chr1	608564	+	40	37	23	17	31	28	19	12
```

Example inputs are provided with the code.

