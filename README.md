# hydi: calling differential DNA hydroxymethylation

## Description

hydi implements methods for calling differential DNA hydroxymethylation from oxWGBS-Seq and oxRRBS-Seq experiments.

## Background

Briefly, ox[RR|WG]BS experiments involve sequencing the DNA of each biological sample twice - following two distinct chemical treatments. 
The first treatment with sodium bisulfite (BS) is carried out to measure the sum of 5-methylcytosine (5mC) and 5-hydroxymethyl-cytosine (5hmC). 
Thus, we denote the sum of both modifications by 5modC and `5modC=5mC+5hmC`.
The second treatment involves an oxidation reaction prior to the sodium bisulfite conversion (oxBS). Subsequent sequencing yields information on 5mC-levels. 
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
git clone https://github.com/Hoffmann-Lab/hydi.git
cd hydi
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
| group1     | formatted input file for the 1st group (G1). The format is described below.                      | 
| group2     | formatted input file for the 2nd group (G2). The format is described below.                      | 
| output     | path or filename for the output. If option is omitted, the output will be dumped to standard output  |

Please note, the assignment of the sample files to the groups matters. The calculation of estimates and confidence intervals for differential hydroxymethylation uses G1 as a reference. Thus, negative values
indicate a loss of hydroxymethylation in G1 as compared to G2. Conversely, positive values indicate a gain of hydroxymethylation in G1. For instance, in a comparison of a group of cancer samples versus healthy controls it is recommended to assign the cancer data to G1 to make the interpretation more intuitive.

### Control

| parameter  | description      |
|----------- |-------------------------------------|
| maxiter    | control of the number of iterations for the gradient descent used in the estimation of maximum likelihoods  | 
| alpha      | parameter controlling the level of confidence intervals and used for flags | 
| epsilon    | parameter controlling the numerical precision of all calculations | 

Please note, while increasing maxiter and decreasing epsilon may provide a better precision, the runtime of hydi may be affected substantially. Convenience flags of hydis output are set based on the value of `alpha`. Furthermore, hydi calculates `1-alpha` confidence intervals.


## Input data format

hydi calls differentially hydroxymethylated cytosines (dhmC) by comparing two groups of samples analysed with
oxRRBS or oxWGBS sequencing experiments. For each group, hydi needs to be given a separate (gzipped) tab-separated file. Each file summarizes 
the count data obtained from the individual sequencing libraries in a specific order. 

The first line is a header with column descriptions and sample identifiers. 
Each of the following lines represents a single C in the genome. The first three fields of each line specify the coordinates of the C.


All following fields contain the count data. **For each sample, there are four fields**, i.e. two for each sequencing run:

|  column  | description |
| :---------------------------------------------------| :------------------------------------| 
| chrom                            | chromosome (alphanumeric) | 
| pos                              | position of C (integer) | 
| strand                           | strand of C ("+" or "-") |   
|                                  |                            |
| ID of BS-Seq (coverage)          |  number of reads aligned to coordinate in BS-Seq | 
| ID of BS-Seq (non-conversion)    | number of non-converted Cs at coordinate in BS-Seq | 
| ID of oxBS-Seq (coverage)        |  number of reads aligned to coordinate in oxBS-Seq | 
| ID of oxBS-Seq (non-conversion)    | number of non-converted Cs at coordinate in oxBS-Seq |   

where BS and oxBS denote the two distinct treatments briefly described above. For a two samples, 

```
chrom	pos	strand	SRR2074675	SRR2074675	SRR2074676	SRR2074676	SRR2074679	SRR2074679	SRR2074680	SRR2074680
chr1	608564	+	40	37	23	17	31	28	19	12
```
Should counts be unavailable for a sample, e.g. when a specific C is not covered by one of the sequencing runs, the value `NA` may be used. Running examples are provided with the code.

## Output data

Hydi returns the results in a tab-separated text file. Each line represents a single genomic Cs and summarizes the test data for group 1 (G1) and group 2 (G2). The
first three fields hold their respective coordinates (see input data). The remaining 18 fields are
described in the following table:

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
Similarly, the `5hmC`-flag is set for one or both groups if the estimates are positive and significantly different from 0 (based on the chosen `alpha`). 

An example output line looks like this:

```
chr1	434286	+	-0.069985	0.013117	0.096011	0.752027	0.977618	-0.087720	0.015558	0.117299	0.763707	0.951594	0	0	-0.130453	-0.002441	0.125570	0.9
70594	1.000000	0.000000
```

## Recommendations

### Analysis of differential 5hmC

For the analysis of differential hydroxymethylation it is recommended to filter out all Cs with overshoots, i.e. a statistically significant
**negative** hydroxymethylation. Overshoots may be caused by coverage or alignment problems. To do this, the convenience flag `overshoot` (field 14) may be
used. This can be done, for instance, by calling

```{sh}
awk '{if(NR == 1 || $14 == 0) print }' examples/output.txt > examples/output.noovershoot.txt
```

For all downstream analyses, it is recommended to filter the data for a desired fdr cutoff, e.g. fdr <= 0.1, using `fdr_diff` (field 20) and to use the convenience field `est_mindiff` (field 21) to eliminate biologically irrelevant differences. Both values are calculated based on the difference p-value (field 19) and the confidence interval (field 16-18), respectively. For instance, the shell command

```{sh}
awk '{if(NR == 1 || ($14 == 0 && $20 <= 0.1 && $21 >= 0.1)) print }' examples/output.txt > examples/output.noovershoot.filtered.txt
```
excludes sites with significant overshoots and extracts all sites with significant fdr-corrected differential hydroxymethylation and a minimum difference of at least 0.1.   

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

## Complaint department

- steve hoffmann leibniz minus fli de
- robert schwarz leibniz minus fli de

## Tutorial

The following tutorial shows the complete workflow for dhmC analysis from (ox)BS-Seq fastq data to calling differenital hydroxymethylation with hydi.

**Tools that are needed**

The following list contains all tools/scripts that are needed for this tutorial

- SRA Toolkit
- Cutadapt
- segemehl.x
- samtools
- BamUtil
- haarz
- picard
- vcfs2tab.py
- hydi.x

### Preparation and download data 

**Download reference genome**

```bash
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#rename the reference genome
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.p12.fa
```

**Generate segemehl indices**

```bash
segemehl.x -d GRCh38.p12.fa -x GRCh38.p12.fa.segemehl.ctidx -y GRCh38.p12.fa.segemehl.gaidx -F 1
```

**Download data**

Download of the whole genome oxBS-Seq and BS-Seq data of normal and malignant human liver (GEO database, accession number GSE70090).

The liver samples from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA287622&o=acc_s%3Aa are selected and the accession list is downloaded (SraAccList.txt).
The sra-files are downloaded with prefetch and extracted with fastq-dump. For each selected sample a sra file is downloaded. In the following, sample_x is used for simplicity, which means that the commands have to be done for each sample.

```bash

prefetch --option-file SraAccList.txt
fastq-dump --split-files --gzip sample_x.sra

```

### Analysing Data

**Quality trimming with cutadapt**

Clipping of adapters and bases with bad quality. The downloaded fastq-files serve as input.

```bash
cutadapt -a NNAGATCGGAAGAGC -A NNAGATCGGAAGAGC -q 20 -O 10 -m 25 -j 10 -o sample_x_1_clipped.fastq.gz -p sample_x_2_clipped.fastq.gz sample_x_1.fastq.gz sample_x_2.fastq.gz
```

**Mapping with segemehl**

The clipped reads are aligned with `segemehl` and the results are stored in alignment files (.bam-format). The alignment files are sorted with `samtools sort` and subsequently indexed with `samtools index`. 

```bash
segemehl.x -d GRCh38.p12.fa -i GRCh38.p12.fa.segemehl.ctidx -j GRCh38.p12.fa.segemehl.gaidx -q sample_x_1_clipped.fastq.gz -p sample_x_fastq_2_clipped.fastq.gz -F 1 -t 100 -b -o sample_x.bam

samtools sort -T tmp/tmp --threads 20 -o sample_x.sort.bam sample_x.bam

samtools index sample_x.sort.bam
```

**Clipping of Overlaps**

The overlaps within the alignment files are removed with `bam clipOverlap` and the clipped .bam files are sorted and indexed with `samtools`.

```bash
bam clipOverlap --in sample_x.sort.bam --out sample_x.clipped.bam

samtools sort -T tmp/tmp --threads 20 -o sample_x.cl.sort.bam sample_x.clipped.bam

samtools index sample_x.cl.sort.bam
```

**Generating vcf-files with haarz**

For each sample the methylation is called for each cytosine with `haarz` and the results are sorted with `picard`. At the end the sample specific vcfs are merged with `bcftools merge`

```bash
haarz.x callmethyl -d GRCh38.p12.fa -u -b  sample_x.cl.sort.bam > sample_x.unique.vcf 2>>vcfs.stderr

java -jar picard.jar SortVcf I=sample_x.unique.vcf O=sample_x.unique.sorted.vcf.gz

bcftools merge `ls -lm *unique.sorted.vcf.gz | tr -d '\n' | tr ',' ' '`  -Oz -o samples.merged.vcf.gz
```

**Generating tables needed for hydi**

`vcfs2tab.py` needs a text file that contains the assignments of samples to group/stage. In the following is the content of the file for this specific tutorial example.

```
sample_BS	sample_oxBS	stage
SRR2074675	SRR2074676	normal
SRR2074679	SRR2074680	normal
SRR2074683	SRR2074684	normal
SRR2074687	SRR2074688	normal
SRR2074677	SRR2074678	tumor
SRR2074681	SRR2074682	tumor
SRR2074685	SRR2074686	tumor
SRR2074689	SRR2074690	tumor
```

`vcfs2tab.py` runs with python 2.7 and needs as input the sample assignment file, the merged vcf file and a minimal coverage (-c). The minimal coverage determines how much reads must at least cover the respective C. The script does extract all Cs in CpG context and generates two tables in a zipped format that contains in their name the respective group/stage.
You will find the sampleAssignment.txt in the examples directory and `vcfs2tab.py` in the scripts directory. 

```bash
python vcfs2tab.py -s sampleAssignment.txt -v samples.merged.vcf.gz -c 10
```

**Run Hydi**

Finally, run hydi with the two generated tables by `vcfs2tab.py`.

```bash
./hydi.x -a tumor.dat.gz -b normal.dat.gz > hydi.result.out
```
