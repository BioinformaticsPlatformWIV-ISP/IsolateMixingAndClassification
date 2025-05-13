from pathlib import Path


rule get_isolate_sequences_length:
    """
    Extract length of each sequence.
    Output consist of two columns with seq id and length.
    """
    input:
        lambda wildcards: config['location_isolates'][wildcards.isolate_sample]
    output:
        isolates_data_length = Path('{ROOT}') / 'isolates' / '{isolate_sample}.length'
    shell:
        """
        ml seqkit/2.3.1;
        seqkit fx2tab -nil \
        {input} > {output:q}
        """



rule subsample_isolates_random:
    input:
        rules.get_isolate_sequences_length.output.isolates_data_length
    output:
        subsampeled_isolate = Path('{ROOT}') / '{mix_n}' / 'subsample' / 'random' / '{isolate_sample}.txt'
    params:
        total_yield = config["total_yield"],
        relative_abundance = lambda wildcards: config['mixes'][wildcards.mix_n][wildcards.isolate_sample],
        seed = lambda wildcards: sum(ord(char) for char in wildcards.mix_n)
    script:
        "../scripts/random_subsampling.py"


rule sample_isolates:
    input:
        isolate_fasta = lambda wildcards: config["location_isolates"][wildcards.isolate_sample],
        selected_seqids = rules.subsample_isolates_random.output.subsampeled_isolate
    output:
        subsampeled_isolate_fasta=Path('{ROOT}') / '{mix_n}' / 'subsample' / 'random' / '{isolate_sample}.fq'
    shell:
        """
        #!TODO
        ml seqtk/1.4;
        ml seqkit/2.3.1;
        seqtk subseq {input.isolate_fasta:q} {input.selected_seqids:q} | seqkit rename > {output.subsampeled_isolate_fasta:q}
        # ml seqkit/2.3.1;
        # seqkit grep -f <(cut -f 1 {input.selected_seqids}) {input.isolate_fasta} -o {output.subsampeled_isolate_fasta};
        """



