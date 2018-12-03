from python import *


PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
#PATIENT_IDS = "77612"
CLONES = ["1", "2", "3", "4", "5", "6"]

## full set of unique genes pulled out by the script and put into the json ## 
# GENES = ["V1", "V1-18", "V1-2", "V1-24", "V1-3", "V1-38-4", "V1-45", "V1-46", "V1-58", "V1-69", "V1-69-2", "V1-69D", "V1-8", "V1/OR15", "V1/OR15-1", "V1/OR15-5", "V1/OR15-9",
#   "V2", "V2-26", "V2-5", "V2-70", "V2-70D", "V2/OR16-5", 
#   "V3", "V3-11", "V3-13", "V3-15", "V3-16", "V3-20", "V3-21", "V3-23", "V3-25", "V3-30", "V3-30-3", "V3-33", "V3-35", "V3-38", "V3-38-3", "V3-43", "V3-43D", "V3-48", "V3-49", "V3-53", "V3-64", "V3-64D", "V3-66", "V3-7", "V3-72", "V3-73", "V3-74", "V3-9", "V3-NL1", "V3/OR15-7", "V3/OR16", "V3/OR16-10", "V3/OR16-12", "V3/OR16-13", "V3/OR16-6", "V3/OR16-8", "V3/OR16-9", 
#   "V4", "V4-28", "V4-30", "V4-30-2", "V4-30-4", "V4-31", "V4-34", "V4-38-2", "V4-39", "V4-4", "V4-59", "V4-61", "V4/OR15-8", 
#   "V5", "V5-10-1", "V5-51", 
#   "V6", "V6-1", 
#   "V7", "V7-4-1", "V7-81"]


# "3-NL1","3-OR16", "3-OR16-10", "3-OR16-12", "3-OR16-13", "3-OR16-6", "3-OR16-9", "3-38-3", "3-43D", 3, "3-73", "3-66", "collapse to 1: 3-38, 3-16"
GENES = ("77612", ["3-11", "3-13", "3-15", "3-20", 
  "3-21", "3-23", "3-30", "3-30-3", "3-33", 
  "3-43", "3-48", "3-49", "3-53", "3-64", "3-7",
  "3-72", "3-74", "3-9",]
)


rule all:
  input:
    expand(
      "data/{patient_id}/V{v_gene}_dashboard.json",
      patient_id=PATIENT_IDS,
      v_gene=GENES,
      # patient_id=GENES[0],
      # v_gene=GENES[1],
    )

rule unpacked:
  input:
    "data/bcell-phylo_Ver4.tar.gz"
  output:
    #expand("data/input/{patient_id}_{clone}_clone.json", patient_id=GENES[0], clone=CLONES),
    expand("data/input/{patient_id}_{clone}_clone.json", patient_id=PATIENT_IDS, clone=CLONES),
    "data/input/imgt.dat",
    "data/input/imgt.fasta"
  shell:
    "tar xvzf {input} -C data/input"

## good up to here ## 
rule clone_json_to_unaligned_fasta:
  input:
    json="data/input/{patient_id}_{clone}_clone.json"
  output:
    fasta="data/{patient_id}/clone_{clone}_unaligned.fasta"
  run:
    clone_json_to_unaligned_fasta(input.json, output.fasta, wildcards.clone)
 
rule unique_vs:
  input:
    #expand("data/input/{patient_id}_{clone}_clone.json", patient_id=GENES[0], clone=CLONES),
    expand("data/input/{patient_id}_{clone}_clone.json", patient_id=PATIENT_IDS, clone=CLONES),
    #"data/input/{patient_id}_{clone}_clone.json"
  output:
    all_vs="data/{patient_id}/unique_vs.json",
  run:
    get_unique_vs(input, PATIENT_IDS, CLONES, output.all_vs)

rule imgt_records:
  input:
    dat="data/input/imgt.dat"
  output:
    txt="data/imgt/V{v_gene}/V{v_gene}.txt"
  run:
    extract_imgt_records(input.dat, wildcards.v_gene, output.txt)


rule imgt_information:
  input:
    raw = rules.imgt_records.output.txt
  output:
    nucleotide_fasta="data/imgt/V{v_gene}/nucleotide.fasta",
    protein_fasta="data/imgt/V{v_gene}/protein.fasta",
    json="data/imgt/V{v_gene}/data.json"
  run:
    parse_imgt_record(input.raw, output.nucleotide_fasta, output.protein_fasta, output.json, wildcards.v_gene)

######### 

rule separate_into_regions:
  input:
    fasta=expand("data/{{patient_id}}/clone_{clone}_unaligned.fasta", clone=CLONES)
  output:
    fasta="data/{patient_id}/V{v_gene}_unaligned.fasta"
  run:
    separate_into_regions(input.fasta, output.fasta, wildcards.v_gene)

rule collapse_identical_sequences:
  input:
    fasta=rules.separate_into_regions.output.fasta,
  output:
    fasta="data/{patient_id}/V{v_gene}_unaligned_collapsed.fasta"
  run:
    collapse_identical_sequences(input.fasta, output.fasta)

rule protein_and_corrected_dna:
  input:
    fasta=rules.collapse_identical_sequences.output.fasta
  output:
    aa="data/{patient_id}/V{v_gene}_unaligned_corrected_AA.fasta",
    nuc="data/{patient_id}/V{v_gene}_unaligned_corrected_nuc.fasta"
  run:
    protein_and_corrected_dna(input.fasta, output.aa, output.nuc)

## need to add a try and catch here in bash to just exit when the ##
## files are empty and it is trying to align ##
rule alignments:
  input:
    rules.protein_and_corrected_dna.output.aa
  output:
    "data/{patient_id}/V{v_gene}_aligned_AA.fasta",
  shell:
    "mafft --amino {input} > {output}"

rule codon_maker:
  input:
    aligned_aa=rules.alignments.output[0],
    unaligned_nucleotide=rules.protein_and_corrected_dna.output.nuc
  output:
    codon="data/{patient_id}/V{v_gene}_codon.fasta",
  run:
    protein_alignment_to_codon_alignment(input.aligned_aa, input.unaligned_nucleotide, output.codon)


rule profile_alignment:
  input:
    codon=rules.alignments.output[0],
    germline=rules.imgt_information.output.protein_fasta
  output:
    fasta="data/{patient_id}/V{v_gene}_profile.fasta"
  shell:
    "mafft --add {input.germline} --reorder {input.codon} > {output.fasta}"  

rule gap_trimmer:
  input:
    codon=rules.codon_maker.output.codon
  output:
    trimmed="data/{patient_id}/V{v_gene}_ungapped.fasta",
  run:
    gap_trimmer(input.codon, output.trimmed)

rule indicial_mapper:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.imgt_information.output.json
  output:
    json="data/{patient_id}/V{v_gene}_indices.json"
  run:
    indicial_mapper(input.fasta, input.json, output.json)

rule trees:
  input:
    fasta=rules.gap_trimmer.output.trimmed
  output:
    tree="data/{patient_id}/V{v_gene}.new"
  shell:
    "FastTree -nt {input.fasta} > {output.tree}"

rule v_gene_json:
  input:
    fasta=rules.profile_alignment.output.fasta,
    json=rules.indicial_mapper.output.json,
    tree=rules.trees.output.tree,
  output:
    json="data/{patient_id}/V{v_gene}_dashboard.json"
  run:
    json_for_dashboard(input.fasta, input.json, input.tree, output.json, wildcards)

