#! /usr/bin/env python

## This script is to extract alternative splicing events from gtf file,
## including alternative TSS, retated interon (RI), skipped exon (SE),
## alternative 5'SS and alternative 3'SS.

## Please note complicate events containing
## multiple categories (e.g. both A5'SS and A3'SS)
## was not reported in the result, and excluded
## from the analysis

## Alternative PAS will be infered separately by
## nanopore reads polyA signal

import os, sys, argparse
import datetime
import gffutils
import tempfile
import subprocess as sp
import re
import pandas as pd

def make_GFFdb(File_in):
    if ".gtf" in File_in:
        db = gffutils.create_db(data=File_in,
                                dbfn=":memory:",
                                force=True,
                                id_spec={"gene": "gene_id",
                                         "transcript": "transcript_id"},
                                keep_order=True,
                                merge_strategy="merge",
                                disable_infer_transcripts=True)
        ## Infer introns and added to existing db
        db.update(db.create_introns())
        return db
    else:
        sys.exit("Please provide GTF file as input")


def extract_feature(db, transID, featuretype='intron'):
    featureList = []
    for f in db.children(transID, featuretype = featuretype,
                         order_by='start'):
        featureList.append([f.start, f.end])

    ##if db[transID].strand == "-":
    ##    featureList.reverse()
    ## reverse will make it more complicate
    return featureList

def calc_overlap(A_start, A_end,
                 B_start, B_end):
    A_region = set(range(A_start, A_end+1))
    B_region = set(range(B_start, B_end+1))
    lenOverlap = len(A_region & B_region)
    pA = lenOverlap/len(A_region)
    pB = lenOverlap/len(B_region)

    return lenOverlap, pA, pB

def define_ASS_events(A_exonList, A_intronList,
                      B_exonList, B_intronList,
                      A_transID, B_transID,
                      strand):
    ## Define intron events: A5SS A3SS
    ## Input is start and end of intron A, B
    outList = []
    for i, A_intron in enumerate(A_intronList):
        A_start = A_intron[0]
        A_end = A_intron[1]
        ## If A_intron contains any B_exon
        ## not treat as ASS, but should be SE
        B_exon_in = 0
        for B_exon in B_exonList:
            B_exonStart = B_exon[0]
            B_exonEnd = B_exon[1]
            if(B_exonStart >= A_start and
               B_exonEnd <= A_end): ## B_exon was fully covered by A_intron
                B_exon_in = 1
                break
        if B_exon_in ==1:
            break
        for j, B_intron in enumerate(B_intronList):
            B_start = B_intron[0]
            B_end = B_intron[1]
            A_exon_in = 0
            for A_exon in A_exonList:
                A_exonStart = A_exon[0]
                A_exonEnd = A_exon[1]
                if(A_exonStart >= B_start and
                   A_exonEnd <= B_end):
                    A_exon_in = 1
                    break
            if A_exon_in ==1:
                break
            ## Check A, B overlap
            lenOverlap, pA, pB = calc_overlap(A_start, A_end,
                                              B_start, B_end)

            if lenOverlap > 0:
                ## The tolerance of intron boundary is 10 nt
                if abs(A_start - B_start) < 10 and A_end - B_end > 10:
                    if strand == "+":
                        ## Assign A as downstream usage 3' SS, exclusive feature
                        out_record= ["A3SS",
                                     str(B_start), str(B_end),
                                     str(A_start), str(A_end),
                                     B_transID, A_transID]
                        outList.append(out_record)
                    else:
                        ## Assign A as upsteam usage 5' SS, inclusive feature
                        out_record= ["A5SS",
                                     str(A_start), str(A_end),
                                     str(B_start), str(B_end),
                                     A_transID, B_transID]
                        outList.append(out_record)
                elif abs(A_start - B_start) < 10 and B_end - A_end > 10:
                    if strand == "+":
                        ## Assign B as downstream usage 3' SS, exclusive feature
                        out_record= ["A3SS",
                                     str(A_start), str(A_end),
                                     str(B_start), str(B_end),
                                     A_transID, B_transID]
                        outList.append(out_record)
                    else:
                        ## Assign B as upsteam usage 5' SS, inclusive feature
                        out_record= ["A5SS",
                                     str(B_start), str(B_end),
                                     str(A_start), str(A_end),
                                     B_transID, A_transID]
                        outList.append(out_record)
                elif A_start - B_start > 10 and abs(A_end - B_end) < 10:
                    if strand == "+":
                        ## Assign A as downstream usage 5' SS, exclusive feature
                        out_record= ["A5SS",
                                     str(B_start), str(B_end),
                                     str(A_start), str(A_end),
                                     B_transID, A_transID]
                        outList.append(out_record)
                    else:
                        ## Assign A as upstream usage 3' SS, inclusive feature
                        out_record= ["A3SS",
                                     str(A_start), str(A_end),
                                     str(B_start), str(B_end),
                                     A_transID, B_transID]
                        outList.append(out_record)
                elif B_start - A_start > 10 and abs(A_end - B_end) < 10:
                    if strand == "+":
                        ## Assign B as downstream usage 5' SS, exclusive feature
                        out_record= ["A5SS",
                                     str(A_start), str(A_end),
                                     str(B_start), str(B_end),
                                     A_transID, B_transID]
                        outList.append(out_record)
                    else:
                        ## Assign B as upstream usage 3' SS, inclusive feature
                        out_record= ["A3SS",
                                     str(B_start), str(B_end),
                                     str(A_start), str(A_end),
                                     B_transID, A_transID]
                        outList.append(out_record)
            else:
                pass
    return outList

def define_RI_events(A_exonList, A_intronList,
                     B_exonList, B_intronList,
                     A_transID, B_transID):
    ## This define RI events
    ## exonList is a list [start, end] of exons
    ## e.g.:
    ## A: [[1,10],[30,50],[70,90]]
    ## always ordered by their genomic coordinates

    ## Events were defined as any intron fully covered
    ## by one of the exons from the other isoform's exonList
    outList = []
    ## Search A's intron from B's exon list
    for A_intron in A_intronList:
        intron_start = A_intron[0]
        intron_end = A_intron[1]
        for exon in B_exonList:
            exon_start = exon[0]
            exon_end = exon[1]
            ## Generate regions
            intron_region = set(range(intron_start, intron_end+1))
            exon_region = set(range(exon_start, exon_end+1))
            lenOverlap = len(intron_region & exon_region)

            if lenOverlap/len(intron_region) == 1:
                ## According to MISO paper, proportion of isoform with
                ## retained intron is PSI (percentage of inclusion isoform)
                out_record= ["RI",
                             str(exon_start), str(exon_end),
                             str(intron_start), str(intron_end),
                             B_transID, A_transID]
                outList.append(out_record)

    ## Search B's intron from A's exon list
    for B_intron in B_intronList:
        intron_start = B_intron[0]
        intron_end = B_intron[1]
        for exon in A_exonList:
            exon_start = exon[0]
            exon_end = exon[1]
            ## Generate regions
            intron_region = set(range(intron_start, intron_end+1))
            exon_region = set(range(exon_start, exon_end+1))
            lenOverlap = len(intron_region & exon_region)

            if lenOverlap/len(intron_region) == 1:
                out_record= ["RI",
                             str(exon_start), str(exon_end),
                             str(intron_start), str(intron_end),
                             A_transID, B_transID]
                outList.append(out_record)

    return outList


def define_SE_events(A_exonList, A_intronList,
                     B_exonList, B_intronList,
                     A_transID, B_transID):
    ## Events were defined as any exon fully covered
    ## by one of the introns from the other isoform's intronList
    outList = []
    for A_exon in A_exonList:
        exon_start = A_exon[0]
        exon_end = A_exon[1]
        for intron in B_intronList:
            intron_start = intron[0]
            intron_end = intron[1]
            ## Generate regions
            exon_region = set(range(exon_start, exon_end+1))
            intron_region = set(range(intron_start, intron_end+1))
            lenOverlap = len(exon_region & intron_region)

            if lenOverlap/len(exon_region) == 1:
                out_record= ["SE",
                             str(exon_start), str(exon_end),
                             str(intron_start), str(intron_end),
                             A_transID, B_transID]
                outList.append(out_record)

    for B_exon in B_exonList:
        exon_start = B_exon[0]
        exon_end = B_exon[1]
        for intron in A_intronList:
            intron_start = intron[0]
            intron_end = intron[1]
            exon_region = set(range(exon_start, exon_end+1))
            intron_region = set(range(intron_start, intron_end+1))
            lenOverlap = len(exon_region & intron_region)
            if lenOverlap/len(exon_region) == 1:
                out_record= ["SE",
                             str(exon_start), str(exon_end),
                             str(intron_start), str(intron_end),
                             B_transID, A_transID]
                outList.append(out_record)

    return outList

def define_ATSS_events(A_exonList, B_exonList,
                       A_transID, B_transID,
                       strand):
    ## Start of the first exon will be compared
    ## 25 nt difference cut-off will be used
    outList = []
    if strand == "+":
        A_start = A_exonList[0][0]
        A_end = A_exonList[0][1]
        B_start = B_exonList[0][0]
        B_end = B_exonList[0][1]
        A_region = set(range(A_start, A_end+1))
        B_region = set(range(B_start, B_end+1))
        lenOverlap = len(A_region & B_region)
        if (abs(A_start - B_start) > 25 and
            (lenOverlap/len(A_region) < 0.9 or
             lenOverlap/len(B_region) < 0.9)):
            if A_start < B_start:
                out_record = ["ATSS",
                              str(A_start), str(A_start+1),
                              str(B_start), str(B_start+1),
                              A_transID, B_transID]
                outList.append(out_record)
            else:
                out_record = ["ATSS",
                              str(B_start), str(B_start+1),
                              str(A_start), str(A_start+1),
                              B_transID, A_transID]
                outList.append(out_record)
    else:
        A_start = A_exonList[-1][0]
        A_end = A_exonList[-1][1]
        B_start = B_exonList[-1][0]
        B_end = B_exonList[-1][1]
        A_region = set(range(A_start, A_end+1))
        B_region = set(range(B_start, B_end+1))
        lenOverlap = len(A_region & B_region)
        if (abs(A_end - B_end) > 25 and
            (lenOverlap/len(A_region) < 0.9 or
             lenOverlap/len(B_region) < 0.9)):
            if A_end > B_end:
                out_record = ["ATSS",
                              str(A_end-1), str(A_end),
                              str(B_end-1), str(B_end),
                              A_transID, B_transID]
                outList.append(out_record)
            else:
                out_record = ["ATSS",
                              str(B_end-1), str(B_end),
                              str(A_end-1), str(A_end),
                              B_transID, A_transID]
                outList.append(out_record)

    return outList

def detect_events(db):
    eventList = []
    for gene in db.all_features(featuretype="gene", order_by="seqid"):
        geneID = gene['gene_id'][0]
        geneStrand = gene.strand
        transcriptList = []
        for transcript in db.children(gene, featuretype="transcript", order_by="start"):
            transID = transcript['transcript_id'][0]
            transStrand = transcript.strand
            transcriptList.append([transID, transStrand])

        ## pairwise comparison of transcripts
        for i, transcript_1 in enumerate(transcriptList):
            transID_t1 = transcript_1[0]
            strand_t1 = transcript_1[1]
            for j, transcript_2 in enumerate(transcriptList):
                transID_t2 = transcript_2[0]
                strand_t2 = transcript_2[1]
                if i > j and strand_t1 == strand_t2:
                    exonList_t1 = extract_feature(db, transID_t1, 'exon')
                    intronList_t1 = extract_feature(db, transID_t1, 'intron')
                    exonList_t2 = extract_feature(db, transID_t2, 'exon')
                    intronList_t2 = extract_feature(db, transID_t2, 'intron')
                    ## Extract ATSS
                    ATSS_out = define_ATSS_events(exonList_t1, exonList_t2,
                                                  transID_t1, transID_t2, strand_t1)
                    ## Extract ASS
                    ASS_out = define_ASS_events(exonList_t1, intronList_t1,
                                                exonList_t2, intronList_t2,
                                                transID_t1, transID_t2, strand_t1)
                    ## Extract RI
                    RI_out = define_RI_events(exonList_t1, intronList_t1,
                                              exonList_t2, intronList_t2,
                                              transID_t1, transID_t2)
                    ## Extract SE
                    SE_out = define_SE_events(exonList_t1, intronList_t1,
                                              exonList_t2, intronList_t2,
                                              transID_t1, transID_t2)
                    ## Add events to final list
                    locus_eventList = ATSS_out + ASS_out + RI_out + SE_out
                    locus_eventList = list(map(lambda x: [geneID] + x, locus_eventList))
                    eventList = eventList + locus_eventList

    return eventList # redundant event list

def compare_events(eventType,
                   event_A, event_B):
    ## event_A and event_B is a list of
    ## event informaiont.
    ## e.g.:
    ## A: [InclusionStart, InclusionEnd,
    ##     ExclusionStart, ExclusionEnd,
    ##     InclusionTrans, ExclusionTrnas]
    ## It is the value of events_dic
    A_inclusionStart = int(event_A[0])
    A_inclusionEnd = int(event_A[1])
    A_exclusionStart = int(event_A[2])
    A_exclusionEnd = int(event_A[3])
    B_inclusionStart = int(event_B[0])
    B_inclusionEnd = int(event_B[1])
    B_exclusionStart = int(event_B[2])
    B_exclusionEnd = int(event_B[3])

    if eventType == "ATSS":
        if(abs(B_inclusionStart - A_inclusionStart)<=10 and
           abs(B_exclusionStart - A_exclusionStart)<=10):
            ## merge
            return True
        else:
            return False
    elif eventType == "SE":
        if(B_inclusionStart == A_inclusionStart and
           B_inclusionEnd == A_inclusionEnd and
           abs(B_exclusionStart - A_exclusionStart)<=10 and
           abs(B_exclusionEnd - A_exclusionEnd)<=10):
            # merge
            return True
        else:
            return False
    elif eventType == "RI":
        if(abs(B_exclusionStart - A_exclusionStart)<=10 and
           abs(B_exclusionEnd - A_exclusionEnd)<=10):
            # merge
            return True
        else:
            return False
    else:
        return False


def group_events(eventType, eventList):
    ## This function return grouping formation
    ## index of events will be grouped into list
    ## returned groups will be:
    ## e.g. [[0],[1,2],[3]]
    groups = []
    for i, event_B in enumerate(eventList):
        found_group = False
        for group in groups:
            for j in group:
                event_A = eventList[j]
                compare_result = compare_events(eventType,
                                                event_A, event_B)
                if compare_result:
                    group.append(i)
                    found_group = True
                    break
            ## Only assign to one group
            if found_group:
                break
        if not found_group:
            groups.append([i])
    return groups


def merge_events(eventList, groups):
    mergedList = []
    for group in groups:
        tmplist = []
        for idx in group:
            tmplist.append(eventList[idx])
        tmpdf = pd.DataFrame(tmplist)
        inclusionStart = tmpdf[0].min()
        inclusionEnd = tmpdf[1].max()
        exclusionStart = tmpdf[2].min()
        exclusionEnd = tmpdf[3].max()
        inclusionTrans = ",".join(list(map(str, tmpdf[4].tolist())))
        exclusionTrans = ",".join(list(map(str, tmpdf[5].tolist())))
        mergedList.append([inclusionStart, inclusionEnd,
                           exclusionStart, exclusionEnd,
                           inclusionTrans, exclusionTrans])
    return mergedList


def main():
    parser = argparse.ArgumentParser(description='Extract Alternative Splicing Events.')
    parser.add_argument('-i', '--input',
                        help='Input GTF file')
    parser.add_argument('-t', '--tmpfile',
                        help='Output a tempfile for transcript pairwise comparison')
    parser.add_argument('-o', '--output',
                        help='Output events list')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    db = make_GFFdb(args.input)
    ## Generate redundant event list
    eventList = detect_events(db)
    ## Generate collapsed events list with bedtools
    ## Identical events will be collapsed
    if args.tmpfile:
        with open(args.tmpfile, 'w') as ttf:
            for record in eventList:
                ttf.write("%s\n" % "\t".join(record))

    with tempfile.NamedTemporaryFile('wt') as tf:
        for record in eventList:
            tf.write("%s\n" % "\t".join(record))
        tf.flush()
        tf.seek(0) # Reset file point to the beginning of the file
        proc = sp.Popen('cat - | sort -k 1 -k 2 -k 3n -k 4n -k 5n -k 6n | bedtools groupby -g 1,2,3,4,5,6 -c 7,8 -o distinct',
                shell = True, stdin = tf, stdout = sp.PIPE)
        stdout_value = proc.communicate()[0]
        collapsed_events = stdout_value.decode()
        collapsed_events_list = []
        collapsed_events_lineSplit = collapsed_events.strip().split('\n')
        for line in collapsed_events_lineSplit:
            collapsed_events_list.append(line.split('\t'))

    ## Covert data into dic
    events_dic = {}
    for line in collapsed_events_list:
        gene_locus = line[0]
        eventType = line[1]
        InclusionStart = line[2]
        InclusionEnd = line[3]
        ExclusionStart = line[4]
        ExclusionEnd = line[5]
        InclusionTrans = line[6]
        ExclusionTrnas = line[7]
        if (gene_locus, eventType) in events_dic:
            events_dic[(gene_locus, eventType)].append([InclusionStart, InclusionEnd,
                                                        ExclusionStart, ExclusionEnd,
                                                        InclusionTrans, ExclusionTrnas])
        else:
            events_dic[(gene_locus, eventType)] = []
            events_dic[(gene_locus, eventType)].append([InclusionStart, InclusionEnd,
                                                        ExclusionStart, ExclusionEnd,
                                                        InclusionTrans, ExclusionTrnas])

    merged_events_list = []
    ## group and merge close events
    for key, value in events_dic.items():
        eventType = key[1]
        eventList = value
        group_info = group_events(eventType, eventList)
        mergedList = merge_events(eventList, group_info)
        for record in mergedList:
            ## Add gene and event type information
            merged_events_list.append([key[0]] + [eventType] + record)

    ## Write to output file
    with open(args.output, 'w') as of:
        ## Add header
        of.write("%s\n" % "\t".join(['geneID', 'EventType',
                                     'InclusionFeature_start', 'InclusionFeature_end',
                                     'ExclusionFeature_start', 'ExclusionFeature_end',
                                     'InclusionTransID', 'ExclusionTransID']))
        for record in merged_events_list:
            record = list(map(str, record))
            ## Get unique TransID in field 7,8:
            ## inclusion trans, exclusion trans
            record[6]=",".join(list(set(re.split(',|\|', record[6]))))
            record[7]=",".join(list(set(re.split(',|\|', record[7]))))
            of.write("%s\n" % "\t".join(record))


if __name__ == "__main__":
    main()
