import pandas as pd
import io
import os
from Bio import SeqIO
import re
import warnings
import argparse
import sys

warnings.simplefilter(action='ignore', category=FutureWarning)

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

parser = argparse.ArgumentParser()

parser.add_argument('avinput')

parser.add_argument('GFF')

parser.add_argument('FASTA')

args = parser.parse_args()

#meta = pd.read_csv(args.metadata)

#df = pd.read_csv(args.all_summary_stats)

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

for l in range(len(variantFunction)):
    if((variantFunction.iloc[l,2]) == "exonic"):
        variantFunction.iloc[l, 0] = "line" + str(l+1)
        for m in range(len(gff.index)):
            if(gff.iloc[m, 2] == "gene"):
                if(variantFunction.iloc[l,5] >= gff.iloc[m,3] and variantFunction.iloc[l,5] <= gff.iloc[m,4]):

                    proteinName = variantFunction.iloc[l, 3]
                    proteinName2 = proteinName.replace("gene:", "transcript:")

                    if((variantFunction.iloc[l,3] in gff.iloc[m,8])):

                        # proteinName = variantFunction.iloc[l, 3]
                        # proteinName2 = proteinName.replace("gene:", "transcript:")

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

                    #check if stopgain, not checking if stopgain in frameshift
                        if (len(variantFunction.iloc[l, 12]) == len(variantFunction.iloc[l, 13])):

                            # proteinName = variantFunction.iloc[l, 3]
                            # proteinName2 = proteinName.replace("gene:", "transcript:")

                            aminoNum = int(variantFunction.iloc[l, 10]) + 1 - int(gff.iloc[m, 3])

                            for o in range(len(gff.index)):

                                if (gff.iloc[o, 2] == "gene"):
                                    #if((gff.iloc[o, 3] < gff.iloc[m, 3]) and (gff.iloc[o, 8] == gff.iloc[m, 8])):
                                    if((gff.iloc[o, 3] < gff.iloc[m, 3]) and (gff.iloc[o, 8] == gff.iloc[m, 8])):

                                        #slipageNum = (int(gff.iloc[o, 4]) - int(gff.iloc[m, 3])) + 1
                                        slipageNum = (int(gff.iloc[m, 3]) - int(gff.iloc[o, 4]))

                                        if(slipageNum > 1):

                                            aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3]))

                                        if(slipageNum == 0):

                                            aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3])) + 1

                                        if(slipageNum < 0):

                                            aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3])) + 1 + slipageNum

                                        # aminoNum = int(aminoNum) + (int(gff.iloc[o, 4]) - int(gff.iloc[o, 3])) + slipageNum

                            #proteinNum = str(int(((variantFunction.iloc[l, 10]) + 1 - int(gff.iloc[m, 3])) / 3) + (int(((variantFunction.iloc[l, 10]) + 1 - int(gff.iloc[m, 3])) % 3 > 0)))
                            proteinNum = str(int(aminoNum/3) + ((aminoNum % 3 ) > 0))
                            aminoNum = str(aminoNum)

                            if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 0):
                                if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-" and fasta[int(variantFunction.iloc[l,10]+1)] !=  "-"):

                                    if 'N' in fasta[int(variantFunction.iloc[l,10])-1] or 'N' in fasta[int(variantFunction.iloc[l,10])] or 'N' in fasta[int(variantFunction.iloc[l,10])+1] or 'N' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    elif 'R' in fasta[int(variantFunction.iloc[l,10])-1] or 'R' in fasta[int(variantFunction.iloc[l,10])] or 'R' in fasta[int(variantFunction.iloc[l,10])+1] or 'R' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    else:

                                        before = translate(fasta[int(variantFunction.iloc[l,10])-1] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)])
                                        after = translate(variantFunction.iloc[l,8] + fasta[int(variantFunction.iloc[l,10])] + fasta[int(variantFunction.iloc[l,10]+1)])

                                    if(before == "STOP"):
                                        before = "*"
                                    if (after == "STOP"):
                                        after = "*"

                                    #variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","


                            if((int(variantFunction.iloc[l,10]) - int(gff.iloc[m,3]))%3 == 1):
                                if(variantFunction.iloc[l,8] !=  "-" and fasta[int(variantFunction.iloc[l,10]-2)] !=  "-" and fasta[int(variantFunction.iloc[l,10])] !=  "-"):

                                    if 'N' in fasta[int(variantFunction.iloc[l,10])-2] or 'N' in fasta[int(variantFunction.iloc[l,10])-1] or 'N' in fasta[int(variantFunction.iloc[l,10])] or 'N' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    elif 'R' in fasta[int(variantFunction.iloc[l,10])-2] or 'R' in fasta[int(variantFunction.iloc[l,10])-1] or 'R' in fasta[int(variantFunction.iloc[l,10])] or 'R' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    else:

                                        before = translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + fasta[int(variantFunction.iloc[l, 10]) - 1] + fasta[int(variantFunction.iloc[l, 10])])
                                        after = translate(fasta[int(variantFunction.iloc[l, 10]) - 2] + variantFunction.iloc[l, 8] + fasta[int(variantFunction.iloc[l, 10])])

                                    if(before == "STOP"):
                                        before = "*"
                                    if (after == "STOP"):
                                        after = "*"

                                    #variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","


                            if ((int(variantFunction.iloc[l, 10]) - int(gff.iloc[m, 3])) % 3 == 2):
                                if (variantFunction.iloc[l, 8] != "-" and fasta[int(variantFunction.iloc[l, 10] - 2)] != "-" and fasta[int(variantFunction.iloc[l, 10] - 3)] != "-"):

                                    if 'N' in fasta[int(variantFunction.iloc[l,10])-3] or 'N' in fasta[int(variantFunction.iloc[l,10])-2] or 'N' in fasta[int(variantFunction.iloc[l,10])-1] or 'N' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    elif 'R' in fasta[int(variantFunction.iloc[l,10])-3] or 'R' in fasta[int(variantFunction.iloc[l,10])-2] or 'R' in fasta[int(variantFunction.iloc[l,10])-1] or 'R' in (variantFunction.iloc[l,8]):                                    #
                                        print('no')

                                    else:

                                        before = translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + fasta[int(variantFunction.iloc[l, 10])- 1])
                                        after = translate(fasta[int(variantFunction.iloc[l, 10]) - 3] +  fasta[int(variantFunction.iloc[l, 10])- 2] + variantFunction.iloc[l, 8])

                                    if(before == "STOP"):
                                        before = "*"
                                    if (after == "STOP"):
                                        after = "*"

                                    #variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","

                                    if(before == after):
                                        variantFunction.iloc[l, 1] = "synonymous SNV"

                                    if(before != after):
                                        variantFunction.iloc[l, 1] = "nonsynonymous SNV"

                                    if(after == "*"):
                                        variantFunction.iloc[l, 1] = "stopgain"

                                    if(before == "*" and after != "*"):
                                        variantFunction.iloc[l, 1] = "stoploss"

                                    variantFunction.iloc[l, 3] = proteinName + ":" + proteinName2 + ":exon1:c." + variantFunction.iloc[l, 7] + aminoNum + variantFunction.iloc[l, 8] + ":p." + before + proteinNum + after + ","


variantFunction = variantFunction[variantFunction.iloc[:,2] == "exonic"]

variantFunction.drop(variantFunction.columns[2], axis=1,inplace=True)

variantFunctionName = re.sub('.avinput', '', str(variantFunctionName))


# variantFunction.to_csv('varriantfunction_exonic.csv',index = False,header = False)
variantFunction.to_csv('varriantfunction_exonic.csv', index = False,header=False)

if '.fastq' in variantFunctionName:

    variantFunction.to_csv(variantFunctionName + '.exonic_variant_function', index = False, sep='\t', encoding='utf-8',header=False)

else:

    variantFunction.to_csv(variantFunctionName + '.fastq.exonic_variant_function', index = False, sep='\t', encoding='utf-8',header=False)
