
import gzip
import re
import argparse

def parse_options():

    parser=argparse.ArgumentParser(prog="vcfs2tab")
    parser.add_argument('-s', type = argparse.FileType('r'), dest = 'samples')
    parser.add_argument('-v', type = argparse.FileType('r'), dest = 'vcfFile')
    parser.add_argument('-c', type=int, default = 10, dest='minCov')

    args = parser.parse_args() 

    return parser, args

def loadSample(assignment):

    samples = {}
    types = []
    with open('sampleAssignment.txt', 'r') as input_handle:
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

   
def generateNewLine(line, samples, types, minCov):

    line = line.split('\t')
    strand = line[7].split(';')[1][-1]
    newLineType1 = [line[0], line[1], strand]
    newLineType2 = [line[0], line[1], strand]
    
    for entry in samples:
        
        entry = samples[entry].split('\t')
        bs_reads = line[int(entry[3])].split(':')[3]
        bs_non_conv = line[int(entry[3])].split(':')[4]
        oxbs_reads = line[int(entry[4])].split(':')[3]
        oxbs_non_conv = line[int(entry[4])].split(':')[4]
        
        if bs_reads == '.' or oxbs_reads == '.':
            return None
        elif int(bs_reads) < minCov or int(oxbs_reads) < minCov:
            return None

        else:
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
    
    return ['\t'.join(newLineType1), '\t'.join(newLineType2)]

def generateHeader(samples, types):

    header1 = ['chrom', 'pos', 'strand']
    header2 = ['chrom', 'pos', 'strand']
    
    for entry in samples:

        entry = samples[entry].split('\t')
        
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

    return ['\t'.join(header1), '\t'.join(header2)]

def writeGroupFiles(samples, vcfFile, minCov):
    
    samples, types = loadSample(samples)
    lineCounter = 0

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
            
                newLine = generateNewLine(line, samples, types, minCov)
                
                if newLine:

                    out_group1.write(newLine[0] + '\n')
                    out_group2.write(newLine[1] + '\n')

parser, args = parse_options() 

print('Command: vcfs2tab.py -s {} -v {} -c {}'.format(args.samples.name,
    args.vcfFile.name, args.minCov))
writeGroupFiles(args.samples.name, args.vcfFile.name, args.minCov)


