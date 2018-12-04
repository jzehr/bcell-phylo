import os
from itertools import product
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Memory


JOBLIB_CACHE = os.path.join('data', 'joblib')
MEMORY = Memory(JOBLIB_CACHE, verbose=0)


def is_relevant_imgt_record(record):
    """Determines whether or not a given IMGT record is of interest."""
    is_human = record.annotations['organism'] == 'Homo sapiens (human)'
    is_heavy_chain = 'IG-Heavy' in record.annotations['keywords']
    is_vgene = 'V-gene' in record.annotations['keywords']
    is_relevant = is_human and is_heavy_chain and is_vgene
    return is_relevant


def extract_imgt_records(full_imgt_data):
    with open(full_imgt_data) as full_file:
        text = ''
        for line in full_file:
            if line[:2] != '//':
                text += line
                if line[:2] == 'ID':
                    imgt_id = line.split()[1].split(';')[0]
            else:
                filename = 'data/imgt/%s.txt' % imgt_id
                text += '//\n'
                with open(filename, 'w') as record_file:
                    record_file.write(text)
                text = ''


@MEMORY.cache
def get_relevant_imgt_subset(full_imgt_data):
    with open(full_imgt_data) as f:
        imgt_records = SeqIO.InsdcIO.ImgtIterator(f)
        relevant_subset = list(filter(is_relevant_imgt_record, imgt_records))
    return relevant_subset


def harvest_information(record, information, types, key):
    has_feature = key in types
    if has_feature:
        index = types.index(key)
        feature = record.features[index]
        found_match = False
        has_translation = 'translation' in feature.qualifiers
        if has_translation:
            translation = feature.qualifiers['translation'][0]
            for start_adjustment, end_adjustment in product(range(-5, 5), range(-5, 5)):
                start = feature.location.start + start_adjustment
                end = feature.location.end + end_adjustment
                if record.seq[start:end].translate() == translation:
                    found_match = True
                    break
            if found_match:
                information[key + '-START'] = start
                information[key + '-END'] = end
                information[key + '-TRANSLATION'] = translation
    if not has_feature or not has_translation or not found_match:
        information[key + '-START'] = None
        information[key + '-END'] = None 
        information[key + '-TRANSLATION'] = None


def get_imgt_information(full_imgt_data, human_imgt_fasta):
    records = get_relevant_imgt_subset(full_imgt_data)
    SeqIO.write(records, human_imgt_fasta, 'fasta')
    for record in records:
        information = {}
        types = [feature.type for feature in record.features]
        harvest_information(record, information, types, 'V-EXON')
        harvest_information(record, information, types, 'CDR3-IMGT')
        harvest_information(record, information, types, 'FR3-IMGT')
        record_path = 'data/imgt/%s.json' % record.id
        with open(record_path, 'w') as record_file:
            json.dump(information, record_file, indent=4)


def process_blast_result(blast_result, nuc_fasta, aa_fasta, json_information):
    with open(blast_result) as json_file:
        result = json.load(json_file)
    hits = result["BlastOutput2"]["report"]["results"]["search"]["hits"]
    write_empty_file = False
    if len(hits) > 0:
        top_id = hits[0]["description"][0]["id"]
        imgt_json_path = 'data/imgt/%s.json' % top_id
        with open(imgt_json_path) as imgt_json_file:
            imgt_json = json.load(imgt_json_file)
            with open(json_information, 'w') as imgt_information_file:
                json.dump(imgt_json, imgt_information_file, indent=4)

        record = SeqIO.to_dict(SeqIO.parse('data/imgt/sequences.fasta', 'fasta'))[top_id]
        SeqIO.write(record, nuc_fasta, 'fasta')

        top_id = top_id.split('.')[0]
        imgt_path = 'data/imgt/%s.txt' % top_id
        with open(imgt_path) as imgt_file:
            record = next(SeqIO.InsdcIO.ImgtIterator(imgt_file))
        types = [feature.type for feature in record.features]
        if 'V-EXON' in types:
            index = types.index('V-EXON')
            feature = record.features[index]
            if 'translation' in feature.qualifiers:
                exon = Seq(feature.qualifiers['translation'][0])
                record = SeqRecord(
                    seq=exon,
                    id='Germline_' + top_id,
                    description=''
                )
                SeqIO.write(record, aa_fasta, 'fasta')
            else:
                write_empty_file = True
        else:
            write_empty_file = True
    else:
        write_empty_file = True
    if write_empty_file:
        for empty_file_path in [nuc_fasta, aa_fasta, json_information]:
            with open(empty_file_path, 'w') as empty_file:
                pass

