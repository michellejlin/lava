import argparse
import sys


def read_file(file_path):
    cds_list = []
    start_list = []
    stop_list = []
    for line in open(file_path):
        if line[0] != '#':
            if line.split()[2].upper() == 'CDS':
                cds_list.append(line.split()[8].split('=')[1].split('CDS')[0])
                start_list.append(line.split()[3])
                stop_list.append(line.split()[4])
            name=line.split()[0]
    return cds_list, start_list, stop_list, name


if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Use this tool to autocorrect genious GFF output into the shit that '
                                                 'lava only accepts')
    parser.add_argument('file_path', help='Just list the gff file that genious outputs, as long as the cds annotations'
                                          ' were correct in genious this will output a new gff that you can use with '
                                          'lava')
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
    file_path = args.file_path

    cds_list, start_list, stop_list, name = read_file(file_path)

    g = open('new_gff.gff', 'w')
    g.write('##gff-version 3\n##source-version geneious 9.1.7\n')
    for x in range(0, len(cds_list)):
        g.write(name + '\tGeneious\tgene\t' + start_list[x] + '\t' + stop_list[x] + '\t.\t+\t.\tID=gene:' + cds_list[x] + ';biotype=protein_coding\n')
        g.write(name + '\tGeneious\tCDS\t' +  start_list[x] + '\t' + stop_list[x] + '\t.\t+\t0\tID=CDS:' + cds_list[x] + ';Parent=transcript:' + cds_list[x] + ';biotype=protein_coding\n')
        g.write(name + '\tGeneious\ttranscript\t' + start_list[x] + '\t' + stop_list[x] + '\t.\t+\t.\tID=transcript:' + cds_list[x] + ';Parent=gene:' + cds_list[x] + ';biotype=protein_coding\n')