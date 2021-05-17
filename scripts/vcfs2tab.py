# Script to extract count data from individual sequencing runs and storing
# them in group separated tables that are required for hydi.
# The script needs a text file where the sample names are assigned to the
# respective group and the table where the counts of the individual sequencing
# results are stored in a merged manner (vcf format). 
# With -c you can set a minimal Coverage to filter out samples that do not
# reach that coverage. The Flag -m determines the number of samples that
# hast to be at least in both groups, which means that if in a group are
# to less samples, because of low coverage for example, the CG is not 
# considered in both groups.


import gzip
import re
import argparse

def parse_options():

    parser=argparse.ArgumentParser(prog="vcfs2tab")
    parser.add_argument('-s', type = argparse.FileType('r'), dest = 'samples')
    parser.add_argument('-v', type = argparse.FileType('r'), dest = 'vcfFile')
    parser.add_argument('-c', type=int, default = 10, dest='minCov')
    parser.add_argument('-m', type = int, default = 3, dest='minSam')

    args = parser.parse_args() 

    return parser, args

def loadSample(assignment):

    samples = {}
    types = []
    with open(assignment, 'r') as input_handle:
        s = 0
        for line in input_handle:
        
            if s != 0:
                samples[s] = line.strip() +'\tNA\tNA'
                types.append(line.split('\t')[2].strip())
            s +=  1

    return samples, list(set(types))


def vcfGenerator(vcfFile):
    
    with gzip.open(vcfFile, 'r') as f:
        for line in f:
            yield line


def detectColumn(header, samples):
    i = -1
    for col in header.split('\t'):
        i += 1
        for entry in samples:
            
            entrySplitted = samples[entry].split('\t')
           
            if re.search(entrySplitted[0] , col):
                entrySplitted[3] = str(i)
            if re.search(entrySplitted[1] , col):
                entrySplitted[4] = str(i)
            
            samples[entry] = '\t'.join(entrySplitted)

def generateNewLineSingleSamp(line, samples, types, minCov, minSam):


    line = line.split('\t')
    strand = line[7].split(';')[1][-1]
    newLine = [line[0], line[1], strand]
    counter = 0

    for entry in samples:
        
        entry = samples[entry].split('\t')
        bs_reads = line[int(entry[3])].split(':')[3]
        bs_non_conv = line[int(entry[3])].split(':')[4]
        oxbs_reads = line[int(entry[4])].split(':')[3]
        oxbs_non_conv = line[int(entry[4])].split(':')[4]

        # Count samples with a true count
        if bs_reads != '.' and oxbs_reads != '.':
            if int(bs_reads) > minCov and int(oxbs_reads) > minCov:
                    counter += 1

        # substitute dots or counts that not reach the min Coverage 
        # with NA for better readability and to let hydi know
        # to ignore that samples
        if bs_reads == '.' or int(bs_reads) <= minCov:
            bs_reads = 'NA'
            bs_non_conv = 'NA'
        if oxbs_reads =='.' or int(oxbs_reads) <= minCov:
            oxbs_reads = 'NA'
            oxbs_non_conv = 'NA'

        newLine.append(bs_reads)
        newLine.append(bs_non_conv)
        newLine.append(oxbs_reads)
        newLine.append(oxbs_non_conv)
    
    if counter >= minSam:
        return '\t'.join(newLine)
    else:
        return None

def generateNewLine(line, samples, types, minCov, minSam):

    line = line.split('\t')
    strand = line[7].split(';')[1][-1]
    newLineType1 = [line[0], line[1], strand]
    newLineType2 = [line[0], line[1], strand]
    type1Counter = 0
    type2Counter = 0
    
    for entry in samples:
        
        entry = samples[entry].split('\t')
        bs_reads = line[int(entry[3])].split(':')[3]
        bs_non_conv = line[int(entry[3])].split(':')[4]
        oxbs_reads = line[int(entry[4])].split(':')[3]
        oxbs_non_conv = line[int(entry[4])].split(':')[4]

        # Count samples with a true count
        if bs_reads != '.' and oxbs_reads != '.':
            if int(bs_reads) > minCov and int(oxbs_reads) > minCov:
                
                if entry[2] == types[0]:

                    type1Counter += 1
                elif entry[2] == types[1]:
                    type2Counter += 1

        # substitute dots or counts that not reach the min Coverage 
        # with NA for better readability and to let hydi know
        # to ignore that samples
        if bs_reads == '.' or int(bs_reads) <= minCov:
            bs_reads = 'NA'
            bs_non_conv = 'NA'
        if oxbs_reads =='.' or int(oxbs_reads) <= minCov:
            oxbs_reads = 'NA'
            oxbs_non_conv = 'NA'


        if entry[2] == types[0]:
            newLineType1.append(bs_reads)
            newLineType1.append(bs_non_conv)
            newLineType1.append(oxbs_reads)
            newLineType1.append(oxbs_non_conv)
        elif entry[2] == types[1]:
            newLineType2.append(bs_reads)
            newLineType2.append(bs_non_conv)
            newLineType2.append(oxbs_reads)
            newLineType2.append(oxbs_non_conv)
    
    if type1Counter >= minSam and type2Counter >= minSam:
        return ['\t'.join(newLineType1), '\t'.join(newLineType2)]
    else:
        return None

def generateHeader(samples, types):

    header1 = ['chrom', 'pos', 'strand']
    header2 = ['chrom', 'pos', 'strand']
    
    for entry in samples:

        entry = samples[entry].split('\t')
        
        if len(types) == 1:
                
                header1.append('{}_N'.format(entry[0]))
                header1.append('{}_C'.format(entry[0]))
                header1.append('{}_N'.format(entry[1]))
                header1.append('{}_C'.format(entry[1]))
        else:
            if entry[2] == types[0]:
                header1.append('{}_N'.format(entry[0]))
                header1.append('{}_C'.format(entry[0]))
                header1.append('{}_N'.format(entry[1]))
                header1.append('{}_C'.format(entry[1]))
            else:

                header2.append('{}_N'.format(entry[0]))
                header2.append('{}_C'.format(entry[0]))
                header2.append('{}_N'.format(entry[1]))
                header2.append('{}_C'.format(entry[1]))
        
    
    if len(types) == 1:
        return '\t'.join(header1)
    else:
        return ['\t'.join(header1), '\t'.join(header2)]


def writeSingleSample(samples, vcfFile, minCov, minSam, types):

    lineCounter = 0 

    with gzip.open('{}.dat.gz'.format(types[0]), 'w') as output:

            header = generateHeader(samples, types)

            output.write(header+'\n')

            for line in vcfGenerator(vcfFile):

                lineCounter += 1

                if lineCounter % 1000000 == 0:
                    print('{} lines processed ... '.format(lineCounter))
                
                if re.search('^#CHROM', line):
        
                    detectColumn(line, samples)

                if re.search('CC=CG;', line):
            
                    newLine = generateNewLineSingleSamp(line, samples, types, minCov, minSam)
                     
                    if newLine:
                        output.write(newLine + '\n')

def writeGroupFiles(samples, vcfFile, minCov, minSam):
    
    samples, types = loadSample(samples)
    lineCounter = 0
    print(types)
    if len(types) == 2:
        with gzip.open('{}.dat.gz'.format(types[0]), 'w') as out_group1, \
                gzip.open('{}.dat.gz'.format(types[1]), 'w') as out_group2:

            header = generateHeader(samples, types)

            out_group1.write(header[0]+'\n')
            out_group2.write(header[1]+'\n')
    
            for line in vcfGenerator(vcfFile):
                lineCounter += 1

                if lineCounter % 1000000 == 0:
                    print('{} lines processed ... '.format(lineCounter))
                if re.search('^#CHROM', line):
        
                    detectColumn(line, samples)

                if re.search('CC=CG;', line):
            
                    newLine = generateNewLine(line, samples, types, minCov, minSam)
                
                    if newLine:

                        out_group1.write(newLine[0] + '\n')
                        out_group2.write(newLine[1] + '\n')
    else:
        
        writeSingleSample(samples, vcfFile, minCov, minSam, types)

parser, args = parse_options() 

print('Command: vcfs2tab.py -s {} -v {} -c {} -m {}'.format(args.samples.name,
    args.vcfFile.name, args.minCov, args.minSam))
writeGroupFiles(args.samples.name, args.vcfFile.name, args.minCov, args.minSam)


