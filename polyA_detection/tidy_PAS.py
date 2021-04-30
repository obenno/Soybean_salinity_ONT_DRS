#! /usr/bin/env python

import sys
import argparse
import subprocess as sp
import pandas as pd
from datetime import datetime

def get_timeStamp():
    return str(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))

def group_PAS(distance, PAS_positionList):
    ## This function return grouping formation
    ## index of PAS will be grouped into list

    ## PAS_positionList is the list of PAS coordinates and must be sorted

    ## returned groups will be:
    ## e.g. [[0],[1,2],[3]]
    groups = []
    for i, pas_A in enumerate(PAS_positionList):
        found_group = False
        for group in groups:
            for j in group:
                pas_B = PAS_positionList[j]
                if abs(int(pas_A) - int(pas_B)) <= int(distance):
                    group.append(i)
                    found_group = True
                    break
            ## Only assign to one group
            if found_group:
                break
        if not found_group:
            groups.append([i])
    return groups

def merge_PAS(PASList, groups):
    mergedList = []
    for group in groups:
        tmplist = []
        for idx in group:
            tmplist.append(PASList[idx])
        tmpdf = pd.DataFrame(tmplist,
                             columns = ['trans', 'read', 'chr', 'pos'])
        PAS_Start = tmpdf['pos'].min()
        PAS_End = tmpdf['pos'].max()
        PAS_Chr = tmpdf['chr'][0]
        tmpdf_callapsed = tmpdf.groupby('trans', as_index = False).agg(','.join)
        Trans = ",".join(list(map(str, tmpdf_callapsed['trans'].tolist())))
        tmpdf_callapsed = tmpdf.groupby('trans', as_index = False).agg('|'.join)
        reads = ",".join(tmpdf_callapsed.iloc[:,1].tolist())
        mergedList.append([PAS_Chr, PAS_Start, PAS_End,
                           Trans,
                           reads])
    return mergedList

def main():
    parser = argparse.ArgumentParser(description='Extract Alternative Splicing Events.')
    parser.add_argument('-b', '--bam', required=True,
                        help='Input bam file (Nanoporre reads vs transcriptome)')
    parser.add_argument('-p', '--pas', required=True,
                        help='Input PAS tsv file, containing ReadID, PolyA_len, Chr, PAS_positon (extract_PAS.sh ouput)')
    parser.add_argument('-g', '--gene', required=True,
                        help='Input gene transcript correspondence file, including gene coordinates (Tab separated)')
    parser.add_argument('-d', '--distance', default = 10,
                        help='Difine a distance to merge neighbour PASs [default: 10]')
    parser.add_argument('-o', '--output', required=True,
                        help='Output PAS informaiotn')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    print(get_timeStamp(), ": Analysis started.")
    bf= open(args.bam, 'r')
    ## use triple quote to escape
    pre_command = '''samtools view -F 0x4 -F 0x10 -F 0x100 -F 0x800 - | awk '{print $3"\\t"$1}' | awk 'NR==FNR{a[$1]=$3"\\t"$4}NR>FNR{if($2 in a){print $1"\\t"$2"\\t"a[$2]}}' ''' + args.pas + ''' - | awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]"\\t"$0}' ''' + args.gene + ''' - '''
    ##print(pre_command)
    preprocessing = sp.Popen(pre_command, shell = True, stdin = bf, stdout = sp.PIPE)
    stdout_value = preprocessing.communicate()[0]
    processed = stdout_value.decode()
    processedLines = processed.strip().split('\n')
    ## processedLines headers:
    ## geneID transcriptID readsID chr PAS_position
    bf.close()
    print(get_timeStamp(), ": Finished preprocessing.")

    ## Read gene coordinates
    gene_coordinates = {}
    with open(args.gene, 'r') as gf:
        for line in gf.readlines():
            line = line.strip()
            lineList = line.split('\t')
            gene = lineList[1]
            chromosome = lineList[2]
            start = lineList[3]
            end = lineList[4]
            gene_coordinates[gene] = [chromosome, start, end]

    ## Read data into dic
    gene_PAS = {}
    for line in processedLines:
        fields = line.split('\t')
        gene = fields[0]
        gene_chr = gene_coordinates[gene][0]
        gene_start = gene_coordinates[gene][1]
        gene_end = gene_coordinates[gene][2]
        read_id = fields[2]
        read_chr = fields[3]
        read_PAS = fields[4]
        ## Filter reads, discard reads with inconstant
        ## alignment position between vs_genome result and
        ## vs_transcriptome result
        if(read_chr == gene_chr and
           read_PAS >= gene_start and
           read_PAS <= gene_end):
            if gene not in gene_PAS:
                gene_PAS[gene] = []
                gene_PAS[gene].append(fields[1:])
            else:
                gene_PAS[gene].append(fields[1:])
        else:
            pass

    merged_PASList = []
    for gene in gene_PAS:
        positionList = []
        ## Sort PAS[gene] by PAS coordinates
        gene_PAS[gene].sort(key=lambda x: x[3])
        for d in gene_PAS[gene]:
            positionList.append(d[3])
        group_info = group_PAS(int(args.distance), positionList) # PAS will be grouped by 10 nt
        gene_PAS_merged = merge_PAS(gene_PAS[gene], group_info)
        for record in gene_PAS_merged:
            merged_PASList.append([gene] + record)
    print(get_timeStamp(), ": Finished merging PAS.")

    with open(args.output, 'w') as of:
        for line in merged_PASList:
            of.write("%s\n" % "\t".join(line))
    print(get_timeStamp(), ": All Finished.")


if __name__ == "__main__":
    main()
