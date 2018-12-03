import os
import itertools as it
import re
import collections
import json
import itertools as it

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


## had to manually get rid of "V", "V1", "V2", etc ##
def get_unique_vs(input, patients, clones, all_vs):
    #giant_vs = []
    ## have and and clause that says if the re.split is not in this list ##
    bad_eggs = ['V', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7']
    for patient_id in patients:
        vs = []
        for clone in clones:
            json_filename = 'data/input/%s_%s_clone.json' % (patient_id, clone) 
            with open(json_filename) as json_file:
                data = json.load(json_file)
                all_entries = it.chain.from_iterable(data)
                vs += [re.split(',|\*|\|', entry['tag'])[0] for entry in all_entries if entry['size'] > 30]
                #giant_vs += [re.split(',|\*|\|', entry['tag'])[0] for entry in all_entries]

        unique_vs = list(set(vs))
        unique_vs.sort()

        os.makedirs('data/%s/' %  patient_id, exist_ok=True)
        with open('data/%s/unique_vs.json' %  patient_id, 'wt') as output_file:
            json.dump(unique_vs, output_file)

    # unique_giant = list(set(giant_vs))
    # unique_giant.sort()
    # with open('data/unique_vs.json', 'wt') as output_file:
    #     json.dump(unique_vs, output_file)


def clone_json_to_unaligned_fasta(input, output, clone):
    with open(input, 'r') as input_file:
        data = json.load(input_file)
    bad_sequences = 0
    size = 30
    with open(output, 'w') as out:
        for i, item in enumerate(it.chain.from_iterable(data)):
            if int(item["size"]) > size:
                t = Seq(item["tag"]).split('|')[1] #this is nuc seq
                j = str(item["tag"]).split('|')[0] #this is v(d)jc
                c_list = str(item["centroid"]).split(':')
                c = c_list[-1]
            try:
                translated = str((t).translate()[1:-2])
                out.write(''.join(['>seq' + str(i) + '_time-' + str(clone) + '_size-' + str(item["size"]) + '_' + translated + '_' + str(j) + '_' +str(c)]))
            except:
                bad_sequences += 1
    if bad_sequences == 0:
        print('No bad sequencees.')
    else:
        print('%d sequences that would not translate.' % bad_sequences)


def separate_into_regions(input, output, v_gene):
    v_gene_regex = re.compile('._V(' + v_gene + ').')
    bad_searches = []
    bad_search_text = []
    v_sequences = []
    for f in input:
      current_clone = SeqIO.parse(f, 'fasta')
      for sequence in current_clone:
          regex_search = v_gene_regex.search(sequence.id)
          if not regex_search is None:
              v_gene = regex_search.group(1)
              v_sequences.append(sequence)
          else:
              bad_searches.append(sequence.id)
      SeqIO.write(v_sequences, output, 'fasta')
      print('%d searches were bad: %s' % ( len(bad_searches), bad_search_text))


def collapse_identical_sequences(input_fasta, output_fasta):
    ### return true of false, then exit the snakemake ### 
    records = list(SeqIO.parse(input_fasta, 'fasta'))
    collapsed_records = []
    enumerated_records = list(enumerate(records))
    deleted_records = []
    for i, record_i in enumerated_records:
        if not i in deleted_records:
            for j, record_j in enumerated_records[i+1:]:
                time_i = int(record_i.name.split('_')[1].split('-')[1])
                time_j = int(record_j.name.split('_')[1].split('-')[1])
                if record_i.seq == record_j.seq and time_i == time_j:
                    id_i = record_i.name.split('_')[0][3:]
                    id_j = record_j.name.split('_')[0][3:]
                    size_i = int(record_i.name.split('_')[2].split('-')[1])
                    size_j = int(record_j.name.split('_')[2].split('-')[1])
                    new_id = id_i+id_j
                    new_size = size_i+size_j
                    header_portion = '_'.join(record_i.name.split('_')[3:])
                    header_parameters = (new_id, time_i, new_size, header_portion)
                    new_header = 'seq%s_time-%d_size-%d_%s' % header_parameters 
                    record_i.name = new_header
                    deleted_records.append(j)
            collapsed_records.append(record_i)
    with open(output_fasta, 'w') as output_file:
        for record in collapsed_records:
            output_file.write('>%s\n%s\n' % (record.name, str(record.seq)))


def protein_and_corrected_dna(input, output_aa, output_nuc):
    records = list(SeqIO.parse(input, 'fasta'))
    k = []
    j = []
    q = []
    r = []
    s = []
    seqs_with_stops = []
    for seq_record in records:
        x = len(seq_record) % 3.0
        y = int(len(seq_record.seq.translate()))
        z = int(len(seq_record.seq.translate(to_stop=True)))
        # something might be going on here
        s = int(len(seq_record.id.split('_')[4].split('*')[0].replace('/', '-')))
        if x == 0 and y == z and s > 2:
            k.append(seq_record.description)
            q.append(seq_record.description)
            r.append(Seq(str(seq_record.seq)))
            j.append(Seq(str(seq_record.seq)).translate())
        else:
            seqs_with_stops.append('1')

    lines = zip(k,j)
    with open(output_aa, 'w') as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))

    lines = zip(q,r)
    with open(output_nuc, 'w' ) as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))


def protein_alignment_to_codon_alignment(protein_alignment, nucleotide_fasta, output):
    aligned_AA = list(SeqIO.parse(protein_alignment, 'fasta'))
    un_nuc = list(SeqIO.parse(nucleotide_fasta, 'fasta'))
    lines = zip(aligned_AA, un_nuc)

    for line in lines:
      u = []
      me = [str(line[1].seq[i:i+3]) for i in range(0, len(line[1]), 3)]
      this = iter(me)
      for i, j in enumerate(line[0]):
          if '-' in j:
              u.append('-'*3)
          else:
              u.append(next(this))
                
      with open(output, 'a') as out:
          out.write('>' + str(line[0].description) +'\n')
          out.write(''.join(u))
          out.write('\n')


def gap_trimmer(input, output):
    codon_file = list(SeqIO.parse(input, 'fasta'))
    id_codon = []
    seq_codon = []
    for i in codon_file:
        id_codon.append(str(i.id))
        seq_codon.append(list(i.seq))

    temp = np.array(seq_codon)
    gap_fraction = np.sum(temp=='-', axis=0)/temp.shape[0]
    good_sites = gap_fraction < .4
    j = [''.join(row) for row in list(temp[:, good_sites])]

    lines = zip(id_codon,j)
    with open(output, 'w' ) as file:
        for line in lines:
            file.write("{}{}\n{}\n".format( '>', line[0], line[1]))


def indicial_mapper(in_fasta, in_json, out_json):
    profile = list(SeqIO.parse(in_fasta, 'fasta'))
    germline = None
    for seq_record in profile:
        print(seq_record)
        if 'Germline_V' in seq_record.description:
            print('found it')
          #if 'Human' in seq_record.description:
            germline = seq_record
            #print(seq_record)

    germline_np = np.array(list(str(germline.seq)), dtype='<U1')
    is_gap = germline_np == '-'
    profile_indices = np.arange(len(is_gap))
    index_map = profile_indices[~is_gap]

    with open(in_json) as v_gene_json_file:
        imgt_vgene_data = json.load(v_gene_json_file)

    CDR3_profile_coords = (
        int(index_map[imgt_vgene_data['cdr3_protein_start']])+1,
        int(index_map[imgt_vgene_data['cdr3_protein_end']-1])+1
    )

    FR3_profile_coords = (
        int(index_map[imgt_vgene_data['fr3_protein_start']])+1,
        int(index_map[imgt_vgene_data['fr3_protein_end']-1])+1
    )

    output_dict = {
        'CDR3': CDR3_profile_coords,
        'FR3': FR3_profile_coords
    }

    with open(out_json, 'w') as output_json_file:
        json.dump(output_dict, output_json_file)

