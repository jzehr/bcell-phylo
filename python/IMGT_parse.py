import os
import itertools as it
import re
import collections
import json
import itertools as it

import numpy as np
import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq


def extract_imgt_records(path_to_imgt_db, gene,outfile):
    ## make sure that this string is actually getting a viable file ##
    print(gene)
    #patient = str(patient)
    data = open(full_file)
    imgt = list(Bio.SeqIO.InsdcIO.ImgtIterator(data))
    
    # with open(path_to_imgt_db) as full_file:
    #     text = ''
    #     should_save = False
    #     for line in full_file:
    #         try:
    #             if line.split()[0] == 'ID':
    #                 imgt_id = line.split()[1][:-1]
    #                 #print(imgt_id)
    #         except Exception as e:
    #             continue
    #         hit_delimiter = line[:2] == '//'
    #         if not hit_delimiter:
    #             text += line
    #             is_v = str(gene) in line
    #             #is_human = "Homo sapiens" in line or "Human" in line or "H.sapiens" in line
    #             is_human = "Human" in line 
    #             is_immunoglobulin = "immunoglobulin" in line
    #             is_heavy_chain = "heavy chain" in line or "heavy-chain" in line
    #             is_variable = "variable region" in line
    #             is_gene = "gene" in line
    #             is_exon = "exon" in line
    #             is_desired = is_v and is_human and is_immunoglobulin and is_heavy_chain and is_variable and is_gene and is_exon
    #             if is_desired:
    #                 should_save = True
    #                 print('this is the good line: ', line)
    #         elif should_save:
    #             should_save = False
    #             #os.makedirs('data/imgt/%s' % v_gene, exist_ok=True)
    #             #filename = 'data/imgt/%s/V%s.txt' % (patient, gene)
    #             with open(outfile, 'w') as record_file:
    #                 record_file.write(text)
    #         if hit_delimiter:
    #             text = ''
        

## lets set this up so this large function calls smaller ones so we dont have to keep adding rules 
## to the snakemake file
## try to do this earlier as well to reduce riles

# def sequence(line, ):
#     pass

# def exon():
#     pass

# def CDR3():
#     pass

# def FR1():
#     pass

# def inside_FR2():
#     pass

# def inside_FR3():
#     pass

# ## not going to work yet because the indicies are nutty ## 
# def info_getter(line, info):
#     info = str(info)
#     if info in line:

#     pass

def parse_imgt_record(input_record, output_nucleotide_fasta, output_protein_fasta, output_json, v_gene):
    inside_sequence = False
    inside_exon = False
    inside_cdr3 = False
    inside_fr3 = False

    sequence_string = ''
    exon_lines = []
    cdr3_lines = []
    fr3_lines = []
    with open(input_record) as raw_file:
        for line in raw_file:
            if 'V-EXON' in line:
                inside_exon = True
            elif 'L-PART2' in line:
                inside_exon = False
            if 'CDR3-IMGT' in line:
                inside_cdr3 = True
            elif "3'UTR" in line:
                inside_cdr3 = False
            if 'FR3-IMGT' in line:
                inside_fr3 = True
            elif '2nd-CYS' in line:
                inside_fr3 = False

            if inside_sequence:
                sequence_string += ''.join(line.split()[:-1])
            if inside_exon:
                exon_lines.append(line)
            if inside_cdr3:
                cdr3_lines.append(line)
            if inside_fr3:
                fr3_lines.append(line)

            if line[:2] == 'SQ':
                inside_sequence = True
    sequence = Seq(sequence_string)

    exon_start = int(exon_lines[0].split()[2].split('..')[0])
    exon_start += 1
    exon_end = int(exon_lines[0].split()[2].split('..')[1])
    exon_end -= 2
    exon_codon_start = int(exon_lines[1].split('=')[-1])
    exon_translation = exon_lines[2].split('"')[-1]
    exon_translation += exon_lines[3].split()[-1]
    exon_translation += exon_lines[4].split()[1][:-1]
    exon_translation = "".join(exon_translation.split())
    assert str(sequence[exon_start: exon_end].translate()) == exon_translation

    cdr3_nucleotide_start = int(cdr3_lines[0].split()[2].split('..')[0])
    cdr3_nucleotide_start -= 1
    cdr3_nucleotide_end = int(cdr3_lines[0].split()[2].split('..')[1])
    cdr3_nucleotide_end -= 2
    cdr3_translation = cdr3_lines[1].split('"')[1]
    cdr3_nucleotides = sequence[cdr3_nucleotide_start: cdr3_nucleotide_end]
    cdr3_check = str(cdr3_nucleotides.translate())
    assert cdr3_check == cdr3_translation
    cdr3_protein_start = int((cdr3_nucleotide_start - exon_start) / 3)
    cdr3_protein_width = int((cdr3_nucleotide_end - cdr3_nucleotide_start) / 3)
    cdr3_protein_end = cdr3_protein_start + cdr3_protein_width
    assert exon_translation[cdr3_protein_start: cdr3_protein_end] == cdr3_translation

    fr3_nucleotide_start = int(fr3_lines[0].split()[2].split('..')[0])
    fr3_nucleotide_start -= 1
    fr3_nucleotide_end = int(fr3_lines[0].split()[2].split('..')[1])
    fr3_translation = fr3_lines[1].split('"')[1]
    fr3_check = str(sequence[fr3_nucleotide_start: fr3_nucleotide_end].translate())
    assert fr3_check == fr3_translation
    fr3_protein_start = int((fr3_nucleotide_start - exon_start) / 3)
    fr3_protein_width = int((fr3_nucleotide_end - fr3_nucleotide_start) / 3)
    fr3_protein_end = fr3_protein_start + fr3_protein_width
    assert exon_translation[fr3_protein_start: fr3_protein_end] == fr3_translation

    record_information = {
        "exon_start": exon_start,
        "exon_finish": exon_end,
        "exon_codon_start": exon_codon_start,
        "exon_translation": exon_translation,
        "cdr3_nucleotide_start": cdr3_nucleotide_start,
        "cdr3_nucleotide_end": cdr3_nucleotide_end,
        "cdr3_protein_start": cdr3_protein_start,
        "cdr3_protein_end": cdr3_protein_end,
        "cdr3_translation": cdr3_translation,
        "fr3_nucleotide_start": fr3_nucleotide_start,
        "fr3_nucleotide_end": fr3_nucleotide_end,
        "fr3_protein_start": fr3_protein_start,
        "fr3_protein_end": fr3_protein_end,
        "fr3_translation": fr3_translation
    }

    header = "Germline_V%s" % v_gene
    with open(output_nucleotide_fasta, 'w') as nucleotide_fasta_file:
        nucleotide_fasta_file.write('>%s\n%s\n' % (header, sequence))
    with open(output_protein_fasta, 'w') as protein_fasta_file:
        protein_fasta_file.write('>%s\n%s\n' % (header, exon_translation))
    with open(output_json, 'w') as json_file:
        json.dump(record_information, json_file, indent=4)




