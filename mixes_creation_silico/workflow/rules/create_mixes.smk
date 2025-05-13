from pathlib import Path

def input_for_merge_sampled_isolates(wildcards):
    # requires a config containing switches for the whole workflow
    if config.get("length_model"):
        return expand(Path('{ROOT}') / '{mix_n}' / 'subsample' / 'model' / '{isolate_sample}.fq',
                                ROOT=config['output_dir'],
                                mix_n=wildcards.mix_n,
                                isolate_sample=[isolate for isolate, abundance in config['mixes'][wildcards.mix_n].items() if abundance != 0])
    else:
        return expand(Path('{ROOT}') / '{mix_n}' / 'subsample' / 'random' / '{isolate_sample}.fq',
                                ROOT=config['output_dir'],
                                mix_n=wildcards.mix_n,
                                isolate_sample=[isolate for isolate, abundance in config['mixes'][wildcards.mix_n].items() if abundance != 0])


rule merge_sampled_isolates:
    input:
       input_for_merge_sampled_isolates
    output:
        Path(config['output_dir']) / '{mix_n}' / '{mix_n}.fq'
    shell:
        """
        cat {input:q} > {output}
        """

rule gzip_fastq:
    input:
        rules.merge_sampled_isolates.output
    output:
        Path('{ROOT}') / '{mix_n}' / '{mix_n}.fq.gz'
    shell:
        """
        gzip {input}
        """


rule ground_truth:
    output:
        Path('{ROOT}') / '{mix_n}' / 'ground_truth.names'
    params:
        truth = lambda wildcards: {isolate: abundance for isolate, abundance in config['mixes'][wildcards.mix_n].items() if abundance != 0}
    run:
        with open(output[0], 'a') as f:
            for key, value in params.truth.items():
                f.write(f"{key}: {value}\n")

