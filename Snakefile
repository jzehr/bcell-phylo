from python import *


PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]
CLONES = ["1", "2", "3", "4", "5", "6"]
PATIENT_V_PAIRS = get_patient_vgene_pairs(PATIENT_IDS, CLONES)

dashboard_filename = lambda pair: "data/%s/%s/dashboard.json" % (pair['patient_id'], pair['v_gene'])
cluster_filename = lambda pair: "data/%s/%s/cluster.json" % (pair['patient_id'], pair['v_gene'])
rule all:
  input:
    [dashboard_filename(pair) for pair in PATIENT_V_PAIRS],
    [cluster_filename(pair) for pair in PATIENT_V_PAIRS]
  run:
    cleanup(PATIENT_IDS)

rule imgt_information:
  input:
    data="data/input/imgt.dat"
  output:
    fasta="data/imgt/sequences.fasta"
  run:
    extract_imgt_records(input.data)
    get_imgt_information(input.data, output.fasta)

rule blast_database:
  input:
    rules.imgt_information.output.fasta
  output:
    "data/blast"
  shell:
    "makeblastdb -in {input} -dbtype nucl -parse_seqids -out data/blast/db"

rule clone_json_to_unaligned_fasta:
  input:
    json="data/input/{patient_id}_{clone}_clone.json"
  output:
    fasta="data/{patient_id}/clone_{clone}_unaligned.fasta"
  run:
    clone_json_to_unaligned_fasta(input.json, output.fasta, wildcards.clone)

rule separate_into_regions:
  input:
    fasta=expand("data/{{patient_id}}/clone_{clone}_unaligned.fasta", clone=CLONES)
  output:
    fasta="data/{patient_id}/V{v_gene}/unaligned.fasta"
  run:
    separate_into_regions(input.fasta, output.fasta, wildcards.v_gene)

rule collapse_identical_sequences:
  input:
    fasta=rules.separate_into_regions.output.fasta,
  output:
    fasta="data/{patient_id}/V{v_gene}/unaligned_collapsed.fasta"
  run:
    collapse_identical_sequences(input.fasta, output.fasta)

rule protein_and_corrected_dna:
  input:
    fasta=rules.collapse_identical_sequences.output.fasta
  output:
    aa="data/{patient_id}/V{v_gene}/unaligned_corrected_AA.fasta",
    nuc="data/{patient_id}/V{v_gene}/unaligned_corrected_nuc.fasta"
  run:
    protein_and_corrected_dna(input.fasta, output.aa, output.nuc)

rule alignments:
  input:
    nuc=rules.protein_and_corrected_dna.output.nuc,
    aa=rules.protein_and_corrected_dna.output.aa
  output:
    nuc="data/{patient_id}/V{v_gene}/aligned.fasta",
    aa="data/{patient_id}/V{v_gene}/aligned_AA.fasta",
  shell:
    """
      mafft {input.nuc} > {output.nuc}
      mafft --amino {input.aa} > {output.aa}
    """

rule codon_maker:
  input:
    aligned_aa=rules.alignments.output.aa,
    unaligned_nucleotide=rules.protein_and_corrected_dna.output.nuc
  output:
    codon="data/{patient_id}/V{v_gene}/codon.fasta",
  run:
    protein_alignment_to_codon_alignment(input.aligned_aa, input.unaligned_nucleotide, output.codon)

rule consensus:
  input:
    fasta=rules.codon_maker.output.codon
  output:
    fasta="data/{patient_id}/V{v_gene}/consensus.fasta"
  run:
    get_consensus_sequence(input.fasta, output.fasta)

rule blast_imgt:
  input:
    rules.consensus.output.fasta
  output:
    summary="data/{patient_id}/V{v_gene}/blast/result.json",
    consensus="data/{patient_id}/V{v_gene}/blast/result_1.json"
  shell:
    "blastn -db data/blast/db -outfmt 13 -query {input} -out {output.summary}"

rule blast_result:
  input:
    blast_result=rules.blast_imgt.output.consensus
  output:
    nuc_fasta="data/{patient_id}/V{v_gene}/imgt_reference.fasta",
    aa_fasta="data/{patient_id}/V{v_gene}/imgt_reference_AA.fasta",
    json="data/{patient_id}/V{v_gene}/imgt.json"
  run:
    process_blast_result(
      input.blast_result,
      output.nuc_fasta,
      output.aa_fasta,
      output.json
    )

rule profile_alignment:
  input:
    aa=rules.alignments.output.aa,
    germline_aa=rules.blast_result.output.aa_fasta,
    nuc=rules.alignments.output.nuc,
    germline_nuc=rules.blast_result.output.nuc_fasta
  output:
    fasta_nuc="data/{patient_id}/V{v_gene}/profile.fasta",
    fasta_aa="data/{patient_id}/V{v_gene}/profile_AA.fasta"
  shell:
    """
      mafft --add {input.germline_nuc} --reorder {input.nuc} > {output.fasta_nuc} || 
        cp {input.nuc} {output.fasta_nuc}
      mafft --add {input.germline_aa} --reorder {input.aa} > {output.fasta_aa} ||
        cp {input.aa} {output.fasta_aa}
    """

rule gap_trimmer:
  input:
    codon=rules.codon_maker.output.codon
  output:
    trimmed="data/{patient_id}/V{v_gene}/codon_ungapped.fasta",
  run:
    gap_trimmer(input.codon, output.trimmed)

rule indicial_mapper:
  input:
    fasta=rules.profile_alignment.output.fasta_aa,
    json=rules.blast_result.output.json
  output:
    json="data/{patient_id}/V{v_gene}/indices.json"
  run:
    indicial_mapper(input.fasta, input.json, output.json, wildcards.v_gene)

rule trees:
  input:
    fasta=rules.gap_trimmer.output.trimmed
  output:
    tree="data/{patient_id}/V{v_gene}/tree.new"
  shell:
    "FastTree -nt {input.fasta} > {output.tree}"

rule v_gene_json:
  input:
    fasta=rules.profile_alignment.output.fasta_aa,
    json=rules.indicial_mapper.output.json,
    tree=rules.trees.output.tree,
  output:
    json="data/{patient_id}/V{v_gene}/dashboard.json"
  run:
    json_for_dashboard(input.fasta, input.json, input.tree, output.json, wildcards)

rule cluster_data:
  input:
    newick=rules.trees.output.tree
  output:
    json="data/{patient_id}/V{v_gene}/cluster.json"
  run:
    cluster_data(input.newick, output.json)

