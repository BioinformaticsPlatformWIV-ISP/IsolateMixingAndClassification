from pathlib import Path

rule get_training_sequences_length:
    """
    Extract length of each sequence.
    Output consist of two columns with seq id and length.
    """
    input:
        config.get("length_model")
    output:
        training_data_length = Path('{ROOT}') / 'model' / (Path(config.get("length_model")).stem + '.length')
    shell:
        """
        ml seqkit/2.3.1;
        seqkit fx2tab -nil \
        {input} > {output}
        """


rule train_model:
    """
    Create probabilities
    """
    input:
        training_data_length  = rules.get_training_sequences_length.output.training_data_length
    output:
        trained_model = Path('{ROOT}') / 'model' / 'trained_model.tsv',
        trained_model_barplot = Path('{ROOT}') / 'model' / 'trained_model.html'
    params:
        BIN_WIDTH = 100
    run:
        import pandas as pd
        import numpy as np
        import plotly.express as px

        df = pd.read_csv(input.training_data_length, sep='\t', names=['id', 'length'])
        # Create bins based on sequence length with a width of 100
        # bins are n = (n*100, n*100+100]
        df['bin'] = pd.cut(df['length'], bins=np.arange(0, df['length'].max() + params.BIN_WIDTH, params.BIN_WIDTH), labels=False)

        # Calculate the probability of each bin
        bin_probabilities = df['bin'].value_counts(normalize=True)
        df['bin_probability'] = df['bin'].map(bin_probabilities)
        bin_probabilities.to_csv(output.trained_model, sep='\t')

        fig = px.bar(bin_probabilities, y='proportion')
        fig.write_html(output.trained_model_barplot)


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

rule subsample_isolates:
    input:
        training_data_length = rules.get_isolate_sequences_length.output.isolates_data_length,
        trained_model = rules.train_model.output.trained_model
    output:
        subsampeled_isolate = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}.bed',
        subsampeled_isolate_w_bins = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}.bed.full',
    params:
        mix_n = lambda wildcards: wildcards.get("mix_n"),
        isolate_sample= lambda wildcards: wildcards.get("isolate_sample"),
        BIN_WIDTH = 100,
        TOTAL_YIELD = config["total_yield"]
    run:
        import numpy as np
        import pandas as pd
        import bisect
        from collections import defaultdict

        random_number = int(sum(i * ord(char) for i, char in enumerate(wildcards.isolate_sample))) + (params.TOTAL_YIELD - 1_000_000)
        random_generator = np.random.default_rng(seed=random_number)

        abundance = config['mixes'][params.mix_n][params.isolate_sample]
        total_n_bases = round(params.TOTAL_YIELD * (abundance / 100))

        # Index is bin n and value is probability
        model_hist_prob = pd.read_table(input.trained_model, index_col='bin').squeeze()

        df_to_select = pd.read_csv(input[0], sep='\t', names=['id', 'length'])
        # Create bins based on sequence length with a width of 100
        df_to_select['bin'] = pd.cut(df_to_select['length'], bins=np.arange(0, df_to_select['length'].max() + params.BIN_WIDTH,
                                                                            params.BIN_WIDTH),
                                                                            labels=False)

        # Create dict with bin as key and lengths as values
        def generate_index_group(group):
            return random_generator.permutation(group.index).tolist()

        # Dict with bin number as key and indices of reads as values (in list)
        bin_index_pair = df_to_select.groupby('bin').apply(generate_index_group).to_dict()
        # Dict with read index as key and information on the read as value (in dict)
        read_idx_length =  dict(zip(df_to_select.index, [{'start': 0, 'end': length, 'length': length}
                                                         for length in df_to_select['length']
                                                        ]
                                   )
                               )

        def infinite_random_bins(model_hist_prob):
            while True:
                yield random_generator.choice(model_hist_prob.index, p=model_hist_prob)

        def find_nearest(search_in, value):
            # Filter out keys that have empty list values
            valid_keys = [k for k, v in search_in.items() if v]

            # Ensure the keys list is sorted
            sorted_keys = sorted(valid_keys)

            # Find the insertion point for the given number using binary search
            index = bisect.bisect_right(sorted_keys, value)

            if index < len(sorted_keys):
                # There is a key greater than the given number
                higher_keys = sorted_keys[index:]
                return random_generator.choice(higher_keys)
            else:
                # No key greater than the given number, return the closest lower key
                return sorted_keys[-1] if sorted_keys else None


        sampled = {}
        selected_total_n_bases, needed_total_n_bases = 0, 0
        random_bin_gen = infinite_random_bins(model_hist_prob)
        second_round_bins = []
        idx = 0
        track_debug = defaultdict(list)
        while selected_total_n_bases < total_n_bases:
            random_bin = next(random_bin_gen)
            # random_bin is not available in the reads, save it for second round
            if random_bin not in bin_index_pair or not bin_index_pair[random_bin]:
                random_length_from_random_bin = random_generator.integers(low=(random_bin * 100) + 1,
                    high=(random_bin * 100) + 100,
                    endpoint=True)
                second_round_bins.append((random_bin, random_length_from_random_bin, idx))
                selected_total_n_bases += random_length_from_random_bin
            # random_bin is available in the reads. Take the first read in that bin fully
            else:
                selected_random_read = bin_index_pair[random_bin].pop()
                start, end, _ = read_idx_length[selected_random_read].values()
                selected_random_read_name = df_to_select.loc[selected_random_read, 'id']
                sampled[idx] = dict(
                    start=start,
                    end=end,
                    length=end,
                    random_bin=random_bin,
                    updated_random_bin=-1,
                    read_idx=selected_random_read,
                    read_name=selected_random_read_name
                )
                selected_total_n_bases += end
                needed_total_n_bases -= end

            idx += 1
            # print(selected_total_n_bases)

        needed_total_n_bases += selected_total_n_bases
        for random_bin, random_length_from_random_bin, idx in reversed(sorted(second_round_bins)):
            updated_random_bin = find_nearest(bin_index_pair, random_bin)
            selected_random_read = bin_index_pair[updated_random_bin].pop()
            start, end, length = read_idx_length[selected_random_read].values()
            selected_random_read_name = df_to_select.loc[selected_random_read, 'id']
            if updated_random_bin < random_bin:
                sampled[idx] = dict(
                    start=start,
                    end=end,
                    length=end,
                    random_bin=random_bin,
                    updated_random_bin=updated_random_bin,
                    read_idx=selected_random_read,
                    read_name=selected_random_read_name
                )
            else:
                sampled[idx] = dict(
                    start=start,
                    end=start + random_length_from_random_bin,
                    length=random_length_from_random_bin,
                    random_bin=random_bin,
                    updated_random_bin=updated_random_bin,
                    read_idx=selected_random_read,
                    read_name=selected_random_read_name
                )
                updated_start = start+random_length_from_random_bin
                update_length = length - random_length_from_random_bin
                read_idx_length[selected_random_read] = {'start': updated_start, 'end': end, 'length': update_length}
                move_to_bin = int((update_length - 1) / 100)
                if move_to_bin not in bin_index_pair:
                    bin_index_pair[move_to_bin] = []
                bin_index_pair[move_to_bin].append(selected_random_read)
                track_debug[selected_random_read].append({'start': updated_start, 'end': end, 'length': update_length})

            needed_total_n_bases -= random_length_from_random_bin
            # print(needed_total_n_bases)

        results = pd.DataFrame.from_dict(sampled, orient='index')

        results.loc[:, ['read_name', 'start', 'end']].to_csv(output.subsampeled_isolate, sep='\t', index=False, header=False)
        results.loc[:, ['read_name', 'start', 'end', 'length', 'random_bin', 'updated_random_bin']].to_csv(output.subsampeled_isolate_w_bins, sep='\t', index=False)


rule sample_isolates:
    input:
        isolate_fasta = lambda wildcards: config["location_isolates"][wildcards.isolate_sample],
        selected_seqids = rules.subsample_isolates.output.subsampeled_isolate
    output:
        subsampeled_isolate_fasta=Path('{ROOT}') / '{mix_n}' / 'subsample' / 'model' / '{isolate_sample}.fq'
    shell:
        """
        #!TODO
        ml seqtk/1.4;
        ml seqkit/2.3.1;
        seqtk subseq {input.isolate_fasta:q} {input.selected_seqids:q} | seqkit rename > {output.subsampeled_isolate_fasta:q}
        # ml seqkit/2.3.1;
        # seqkit grep -f <(cut -f 1 {input.selected_seqids}) {input.isolate_fasta} -o {output.subsampeled_isolate_fasta};
        """



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

