from Bio import SeqIO
import subprocess
import sys
import re
from pathlib import Path

# Investigate unmapped reads
# From Devon Jensen @ BD

# Assumes that Samtools and bowtie2 are in path

# Gencode transcriptome bowtie2 index is pre-built:
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz
# bowtie2-build -f --threads 4 gencode.v35.transcripts.fa gencode35transcriptome


# Input, either a BAM file (from which unmapped cellular reads will be extracted) (assumes index .bai is already created)
    # or, a fasta file of unmapped reads
inputFile=sys.argv[1]

# Set this to None or empty string to skip the transcriptome alignment
#transcriptomeBowtie2Index = '/resources/Gencode/human/gencode35transcriptome'
transcriptomeBowtie2Index = None


print(f'InputFile: {inputFile}')
print('-----')

baseInputName = Path(inputFile).stem

fastaFile = f'{baseInputName}_cellular_unmapped.fasta'

if inputFile.lower().endswith('bam'):
    print(f'Extracting unmapped reads with valid cell label')
    samtoolsRun = subprocess.run([f'samtools view -f 4 {inputFile} | grep -v "CB:Z:0" | cut -f1,10 | sed "s/^/>/" | tr "\\t" "\\n" > {fastaFile}'], shell=True)
else:
    fastaFile = inputFile



print(f'Analying reads for known features')

sumPercent = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
countGR40 = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
poly15BaseCount = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
CTCTcount = 0
AGAGcount = 0
ILMNadaptorFwd = 0
ILMNadaptorRev = 0

beadv2CL = 0
beadv2CLrc = 0
fiveprimeCL = 0
fiveprimeCLrc = 0
re_beadv2CL = re.compile('GTGA.........GACA')
re_beadv2CLrc = re.compile('TGTC.........TCAC')
re_fiveprimeCL = re.compile('AATG.........CCAC')
re_fiveprimeCLrc = re.compile('GTGG.........CATT')

tsoForward = 0
tsoReverse = 0
anyBeadSeq = 0

countReads = 0
for record in SeqIO.parse(fastaFile, 'fasta'):
    countReads += 1
    isBeadSeq = False
    readLength = len(record.seq)

    countA = record.seq.count('A')
    countC = record.seq.count('C')
    countG = record.seq.count('G')
    countT = record.seq.count('T')

    sumPercent['A'] += countA/readLength
    sumPercent['C'] += countC/readLength
    sumPercent['G'] += countG/readLength
    sumPercent['T'] += countT/readLength

    if countA/readLength > 0.4:
        countGR40['A'] += 1
    if countC/readLength > 0.4:
        countGR40['C'] += 1
    if countG/readLength > 0.4:
        countGR40['G'] += 1
    if countT/readLength > 0.4:
        countGR40['T'] += 1

    if 'AAAAAAAAAAAAAAA' in record.seq:
        poly15BaseCount['A'] += 1 
    if 'CCCCCCCCCCCCCCC' in record.seq:
        poly15BaseCount['C'] += 1 
    if 'GGGGGGGGGGGGGGG' in record.seq:
        poly15BaseCount['G'] += 1 
    if 'TTTTTTTTTTTTTTT' in record.seq:
        poly15BaseCount['T'] += 1 

    if 'CTCTCTCT' in record.seq:
        CTCTcount += 1
        isBeadSeq = True
    if 'AGAGAGAG' in record.seq:
        AGAGcount += 1
        isBeadSeq = True

    if 'CTCTTCCGATCT' in record.seq:
        ILMNadaptorFwd += 1
        isBeadSeq = True
    if 'AGAGCGTCGTGT' in record.seq:
        ILMNadaptorRev += 1
        isBeadSeq = True


    if re_beadv2CL.search(str(record.seq)):
        beadv2CL += 1
        isBeadSeq = True
    if re_beadv2CLrc.search(str(record.seq)):
        beadv2CLrc += 1
        isBeadSeq = True
    if re_fiveprimeCL.search(str(record.seq)):
        fiveprimeCL += 1
        isBeadSeq = True
    if re_fiveprimeCLrc.search(str(record.seq)):
        fiveprimeCLrc += 1
        isBeadSeq = True

    if 'TATGCGTAGTAGGTATGGG' in record.seq:
        tsoForward += 1
        isBeadSeq = True
    if 'CCCATACCTACTACGCATA' in record.seq:
        tsoReverse += 1
        isBeadSeq = True

    if isBeadSeq:
        anyBeadSeq += 1

if transcriptomeBowtie2Index:
    print('Checking for alignments to transcriptome')
    samFilename = f'{baseInputName}_align.sam'
    bowtierun = subprocess.run([f'bowtie2 --local -p 4 --no-hd --score-min G,40,8 -x {transcriptomeBowtie2Index} -f -U {fastaFile} -S {samFilename}'], shell=True)


    alignGeneCount = {}
    alignReadCount = 0
    with open(samFilename, 'r') as samfile:

        for line in samfile:
            _, flag, alignTarget, _, _, _, _, _, _, _, _, *tags = line.split('\t')
            if flag != '4':
                alignReadCount += 1
                _, _, _, _, _, gene, _, _, _ = alignTarget.split('|')
                alignGeneCount[gene] = alignGeneCount.get(gene, 0) + 1

    topGenes = sorted(alignGeneCount, key=alignGeneCount.get, reverse=True)


with open(f'{baseInputName}_stats.csv', 'w') as outfile:

    outfile.write(f'Total Unmapped Reads,{countReads}\n')

    for base in ['A', 'C', 'G', 'T']:
        outfile.write(f'Average percent composition of base {base} per read,{round(sumPercent[base]/countReads*100, 2)}\n')
    outfile.write('\n')

    for base in ['A', 'C', 'G', 'T']:
        outfile.write(f'Percent reads where base {base} greater than 40% of read,{round(countGR40[base]/countReads*100, 2)}\n')
    outfile.write('\n')

    for base in ['A', 'C', 'G', 'T']:
        outfile.write(f'Percent reads which contain 15 consecutive {base},{round(poly15BaseCount[base]/countReads*100, 2)}\n')
    outfile.write('\n')

    outfile.write(f'Percent reads containing CTCTCTCT,{round(CTCTcount/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing AGAGAGAG,{round(AGAGcount/countReads*100, 2)}\n')
    outfile.write('\n')

    outfile.write(f'Percent reads containing ILMNadaptorFwd CTCTTCCGATCT,{round(ILMNadaptorFwd/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing ILMNadaptorRev AGAGCGTCGTGT,{round(ILMNadaptorRev/countReads*100, 2)}\n')
    outfile.write('\n')

    outfile.write(f'Percent reads containing BeadV2 CL         (GTGA.........GACA),{round(beadv2CL/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing BeadV2 CL revComp (TGTC.........TCAC),{round(beadv2CLrc/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing 5prime CL         (AATG.........CCAC),{round(fiveprimeCL/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing 5prime CL revComp (GTGG.........CATT),{round(fiveprimeCLrc/countReads*100, 2)}\n')
    outfile.write(f'Sum percent reads containing cell label,{round((beadv2CL+beadv2CLrc+fiveprimeCL+fiveprimeCLrc)/countReads*100, 2)}\n')

    outfile.write('\n')

    outfile.write(f'Percent reads containing TSO forward TATGCGTAGTAGGTATGGG,{round(tsoForward/countReads*100, 2)}\n')
    outfile.write(f'Percent reads containing TSO revcomp CCCATACCTACTACGCATA,{round(tsoReverse/countReads*100, 2)}\n')
    outfile.write('\n')

    outfile.write(f'Contains any bead based seqeunce,{round(anyBeadSeq/countReads*100, 2)}\n')
    outfile.write('\n')

    if transcriptomeBowtie2Index:
        outfile.write(f'Percent unmapped reads aligning to transcriptome,{round(alignReadCount/countReads*100, 2)}\n')
        outfile.write(f'Number of different genes,{len(alignGeneCount)}\n')
        outfile.write('Top genes:\n')


        for count, gene in enumerate(topGenes, 0):
            if count == 20:
                break
            outfile.write(f',{gene} : {alignGeneCount[gene]}\n')

















