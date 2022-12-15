import pandas as pd
import io
import os
from Bio import SeqIO
import re
import warnings
import argparse
import sys

#warning surpressor for debugging
warnings.simplefilter(action='ignore', category=FutureWarning)

#function used to open vcf
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

#Function used to translated codons into proteins
def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

#Function checks if degenerative bases have possible translations
def degenerative(codon):

    #must use or to make phrase 4 characters to not make it codon anymore
    if "R" in codon:
        codon = codonReplace(codon,"R","G","A")
    elif "Y" in codon:
        codon = codonReplace(codon,"Y","C","T")
    elif "K" in codon:
        codon = codonReplace(codon,"K","G","T")
    elif "M" in codon:
        codon = codonReplace(codon,"M","A","C")
    elif "W" in codon:
        codon = codonReplace(codon,"W","A","T")
    elif "S" in codon:
        codon = codonReplace(codon,"S","G","C")
    elif "B" in codon:
        codon = codonReplace2(codon,"B","G","T","C")
    elif "D" in codon:
        codon = codonReplace2(codon,"D","G","A","T")
    elif "H" in codon:
        codon = codonReplace2(codon,"H","A","C","T")
    elif "V" in codon:
        codon = codonReplace2(codon,"V","G","C","A")
    elif "N" in codon:
        codon = codonReplace3(codon,"N","A","G","C","T")
    return codon

#if degen base has two options
def codonReplace(codon, degen, aminoA, aminoB):

    codonOne = codon.replace(degen , aminoA)
    codonTwo = codon.replace(degen , aminoB)

    codonOne = translate(codonOne)
    codonTwo = translate(codonTwo)

    if codonOne == codonTwo:
        codon = codonOne
    else:
        codon = aminoA + "or" + aminoB

    return codon

#if degen base has three options
def codonReplace2(codon, degen, aminoA, aminoB, aminoC):

    codonOne = codon.replace(degen , aminoA)
    codonTwo = codon.replace(degen , aminoB)
    codonThree = codon.replace(degen , aminoC)

    codonOne = translate(codonOne)
    codonTwo = translate(codonTwo)
    codonThree = translate(codonThree)

    if codonOne == codonTwo and codonOne == codonThree:
        codon = codonOne
    else:
        codon = aminoA + "or" + aminoB + "or" + aminoC
    return codon

#if degen base has four options
def codonReplace3(codon, degen, aminoA, aminoB, aminoC, aminoD):

    codonOne = codon.replace(degen , aminoA)
    codonTwo = codon.replace(degen , aminoB)
    codonThree = codon.replace(degen , aminoC)
    codonFour = codon.replace(degen , aminoD)

    codonOne = translate(codonOne)
    codonTwo = translate(codonTwo)
    codonThree = translate(codonThree)
    codonFour = translate(codonFour)

    if codonOne == codonTwo and codonOne == codonThree and codonOne == codonFour:
        codon = codonOne
    else:
        codon =  aminoA + "or" + aminoB + "or" + aminoC + "or" + aminoD

    return codon

parser = argparse.ArgumentParser()
parser.add_argument('avinput')
parser.add_argument('GFF')
parser.add_argument('FASTA')

args = parser.parse_args()

gff = pd.read_csv(args.GFF, sep='\t',header = None, skiprows=[i for i in range(0,2)])

variantFunction = pd.read_csv(args.avinput, sep="\t", header=None)

variantFunctionName = args.avinput
variantFunctionName = os.path.basename(variantFunctionName)

variantFunctionName = re.sub('.fastq.avinput', '', str(variantFunctionName))

for fasta in SeqIO.parse(args.FASTA, "fasta"):
 print("")
variantFunction.insert(0, "", "")
variantFunction.insert(0, " ", "")

numRegions = 0

for i in range(len(gff.index)):
    if(pd.isnull(gff.iloc[i,3]) != True):
        numRegions = numRegions +1
    else:
        break

test = 0
df = None

#annotates the names of the protein on the section that mutated
for j in range(numRegions):
    if(gff.iloc[j,2] == "CDS"):
        for k in range(len(variantFunction)):

            proteinName = gff.iloc[j, 8]
            proteinName = proteinName.replace("ID=CDS:", "gene:")
            sep = ';'
            proteinName = proteinName.split(sep, 1)[0]

            if(variantFunction.iloc[k,8] >= int(gff.iloc[j,3]) and variantFunction.iloc[k,8] <= int(gff.iloc[j,4])):
                if(variantFunction.iloc[k, 0] == ""):
                    variantFunction.iloc[k, 0] = "exonic"

                    variantFunction.iloc[k, 1] = proteinName

                elif (variantFunction.iloc[k, 0] != ""):
                    if (test == 0):
                        df = pd.DataFrame(variantFunction.loc[[k]])
                        df.iloc[0, 1] = proteinName

                    df2 = pd.DataFrame(variantFunction.loc[[k]])
                    df2.iloc[0, 1] = proteinName

                    if(test != 0):
                        df = df.append(df2)

                    test = test + 1

if df is not None:

    df = df.iloc[1: , :]

    variantFunction = pd.concat([variantFunction, df])

variantFunction.to_csv('varriantfunction.csv',index = False,header = False)
variantFunction.to_csv(variantFunctionName + '.varriant_function', index = False, sep='\t', encoding='utf-8',header=False)

variantFunction.insert(0, "   ", "")
variantFunction.insert(0, "    ", "")


#This section annotates the type of mutation
for l in range(len(variantFunction)):
    if((variantFunction.iloc[l,2]) == "exonic"):
        variantFunction.iloc[l, 0] = "line" + str(l+1)
        for m in range(len(gff.index)):
            if(gff.iloc[m, 2] == "gene"):
                if(variantFunction.iloc[l,5] >= gff.iloc[m,3] and variantFunction.iloc[l,5] <= gff.iloc[m,4]):

                    proteinName = variantFunction.iloc[l, 3]
                    proteinName2 = proteinName.replace("gene:", "transcript:")

                    if((variantFunction.iloc[l,3] in gff.iloc[m,8])):

                        variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c."

                        # adds deletions and insertions, detects if frameshift or not
                        if (len(variantFunction.iloc[l, 12]) > len(variantFunction.iloc[l, 13])):
                            if (len(variantFunction.iloc[l, 7]) % 3 == 0):
                                variantFunction.iloc[l, 1] = "nonframeshift deletion"
                            else:
                                variantFunction.iloc[l, 1] = "frameshift deletion"

                        if (len(variantFunction.iloc[l, 12]) < len(variantFunction.iloc[l, 13])):
                            if (len(variantFunction.iloc[l, 8]) % 3 == 0):
                                variantFunction.iloc[l, 1] = "nonframeshift insertion"
                            else:
                                variantFunction.iloc[l, 1] = "frameshift insertion"

                        if (len(variantFunction.iloc[l, 12]) == len(variantFunction.iloc[l, 13])):

                            aminoNum = int(variantFunction.iloc[l, 10]) + 1 - int(gff.iloc[m, 3])

                            for o in range(len(gff.index)):

                                if (gff.iloc[o, 2] == "gene"):
                                    if((gff.iloc[o, 3] < gff.iloc[m, 3]) and (gff.iloc[o, 8] == gff.iloc[m, 8])):

                                        slipageNum = (int(gff.iloc[m, 3]) - int(gff.iloc[o, 4]))

                                        if(slipageNum > 1):

                                            aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3]))

                                        if(slipageNum == 0):

                                            aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3])) + 1

                                        if(slipageNum < 0):

                                            aminoNum = (int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3])) + 1 + slipageNum)

                            #checks amino acid and protein position
                            proteinNum = str(int(aminoNum/3) + ((aminoNum % 3 ) > 0))
                            aminoNum = str(aminoNum)

                            #If mutation is on first amino acid of codon
                            if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 0):

                                if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-" and fasta[int(variantFunction.iloc[l,10]+1)] !=  "-"):

                                    if("0" in variantFunction.iloc[l,8]):
                                        variantFunction.iloc[l,8] = (variantFunction.iloc[l,8]).replace("0" , variantFunction.iloc[l,13])
                                    if("0" in variantFunction.iloc[l,7]):
                                        (variantFunction.iloc[l,7]) = (variantFunction.iloc[l,7]).replace("0" , variantFunction.iloc[l,12])

                                    before = fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]
                                    after = variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)]

                                    #deals with degenerative bases
                                    before = degenerative(before)
                                    after = degenerative(after)

                                    if (len(before) == 3):
                                        before = translate(fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)])
                                    if (len(after) == 3):
                                        after = translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)])

                                    #formats and renames output excell values of vcf
                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    #checks change to annotate type of mutation
                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*" and before != "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

                            #If mutation is on second amino acid of codon
                            if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 1):

                                if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10]-2)] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-"):

                                    if("0" in variantFunction.iloc[l,8]):
                                        variantFunction.iloc[l,8] = (variantFunction.iloc[l,8]).replace("0" , variantFunction.iloc[l,13])
                                    if("0" in variantFunction.iloc[l,7]):
                                        (variantFunction.iloc[l,7]) = (variantFunction.iloc[l,7]).replace("0" , variantFunction.iloc[l,12])

                                    before = fasta[int(variantFunction.iloc[l, 10]) - 2] + fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])]
                                    after = fasta[int(variantFunction.iloc[l, 10]) - 2] + variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])]

                                    before = degenerative(before)
                                    after = degenerative(after)

                                    if (len(before) == 3):
                                        before = translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])])
                                    if (len(after) == 3):
                                        after = translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])])

                                    #formats and renames output excell values of vcf
                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    #checks change to annotate type of mutation
                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*" and before != "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

                            #If mutation is on third amino acid of codon
                            if ((int(variantFunction.iloc[l, 10]) - int(gff.iloc[m, 3])) % 3 == 2):

                                if (variantFunction.iloc[l, 8] != "-" and fasta[int(variantFunction.iloc[l, 10] - 2)] != "-" and fasta[int(variantFunction.iloc[l, 10] - 3)] != "-"):

                                    if("0" in variantFunction.iloc[l,8]):
                                        variantFunction.iloc[l,8] = (variantFunction.iloc[l,8]).replace("0" , variantFunction.iloc[l,13])
                                    if("0" in variantFunction.iloc[l,7]):
                                        (variantFunction.iloc[l,7]) = (variantFunction.iloc[l,7]).replace("0" , variantFunction.iloc[l,12])

                                    before = fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1]
                                    after = fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8]

                                    before = degenerative(before)
                                    after = degenerative(after)

                                    if (len(before) == 3):
                                        before = translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1])
                                    if (len(after) == 3):
                                        after = translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8])

                                    #formats and renames output excell values of vcf
                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    #checks change to annotate type of mutation
                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*" and before != "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

#remove any proteins that are not exonic
variantFunction = variantFunction[variantFunction.iloc[:,2] == "exonic"]

variantFunction.drop(variantFunction.columns[2], axis=1,inplace=True)

#renaming
variantFunctionName = re.sub('.avinput', '', str(variantFunctionName))
variantFunction.to_csv('varriantfunction_exonic.csv', index = False,header=False)

#produces file name the same if fastq gzipped or not
if '.fastq' in variantFunctionName:

    variantFunction.to_csv(variantFunctionName + '.exonic_variant_function', index = False, sep='\t', encoding='utf-8',header=False)

else:

    variantFunction.to_csv(variantFunctionName + '.fastq.exonic_variant_function', index = False, sep='\t', encoding='utf-8',header=False)
