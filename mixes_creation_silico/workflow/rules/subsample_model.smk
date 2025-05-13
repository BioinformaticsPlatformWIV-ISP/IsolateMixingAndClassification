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
        training_data_length = Path('{ROOT}') / 'model' / '{isolate_sample}.length'
    shell:
        """
        ml seqkit/2.3.1;
        seqkit fx2tab -nil \
        {input} > {output:q}
        """

rule sample_length_from_isolate_based_on_training_set:
    input:
        training_data_length = rules.get_isolate_sequences_length.output.training_data_length,
        trained_model = rules.train_model.output.trained_model
    output:
        subsampeled_isolate = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}.bed',
        subsampeled_isolate_w_bins = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}.bed.full',
        subsampeled_isolate_barplot = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}_seq_per_bin.html',
        subsampeled_isolate_barplot_copies = Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}_copies_per_bin.html',
        subsampeled_isolate_lineplot= Path('{ROOT}') / '{mix_n}'/ 'subsample' / 'model'  / '{isolate_sample}_lines.html',

    params:
        mix_n = lambda wildcards: wildcards.get("mix_n"),
        isolate_sample= lambda wildcards: wildcards.get("isolate_sample"),
        BIN_WIDTH = 100,
        TOTAL_YIELD = config["total_yield"]
    run:
        import numpy as np
        import pandas as pd
        from collections import Counter
        import plotly.express as px
        import plotly.graph_objects as go
        import sys

        random_number = int(sum(i * ord(char) for i, char in enumerate(wildcards.isolate_sample))) + (params.TOTAL_YIELD - 1_000_000)


        abundance = config['mixes'][params.mix_n][params.isolate_sample]
        total_n_reads = round(params.TOTAL_YIELD * (abundance / 100))

        # Index is bin n and value is probability
        model_hist_prob = pd.read_table(input.trained_model, index_col='bin').squeeze()

        df_to_select = pd.read_csv(input[0], sep='\t', names=['id', 'length'])
        # Create bins based on sequence length with a width of 100
        df_to_select['bin'] = pd.cut(df_to_select['length'], bins=np.arange(0, df_to_select['length'].max() + params.BIN_WIDTH,
                                                                            params.BIN_WIDTH),
                                                                            labels=False)

        # Create dict with bin as key and lengths as values
        def generate_index_group(group):
            random_generator = np.random.default_rng(seed=random_number)
            return random_generator.permutation(group.index).tolist()

        bin_index_pair = df_to_select.groupby('bin').apply(generate_index_group, include_groups=False).to_dict()


        def find_nearest(search_in, value):
            # Search for closest higher bin. If not existant, take the closest lower bin
            array = np.fromiter([key for key, val in search_in.items() if val], dtype=int)
            # the only scenario where a higher bin is not existant, is if value is higher than all elements in array.
            # just return the max value is that is
            max_val = np.max(array)
            if value > max_val:
                # print(f'{value} swapped for {max_val}')
                return max_val
            # else return the closest bin lower than value
            closest_lower_bin = np.min(array[value < array])
            # print(f'{value} swapped for {closest_lower_bin}')
            return closest_lower_bin


        sampled = []
        original_bin = []
        start_coordinates = []
        random_generator = np.random.default_rng(seed=random_number)
        random_bins = random_generator.choice(model_hist_prob.index, size=total_n_reads, p=model_hist_prob)
        random_bins_counts = Counter(random_bins)
        # First fetch from bins that exist
        # For the first 10 bins, just take randomly, doesn't matter which sequence, they will be cut later at length 1000
        for random_bin, n in random_bins_counts.items():
            if random_bin < 10:
                isolate_lengths_selected_bin = bin_index_pair[random_bin]
                random_generator = np.random.default_rng(seed=random_number)
                sampled.extend(random_generator.choice(isolate_lengths_selected_bin, size=n))
                original_bin.extend([random_bin] * n)
                random_bins_counts[random_bin] -= n
            else:
                for _ in range(n):
                    if random_bin not in bin_index_pair or not bin_index_pair[random_bin]:
                        print(f'bin {random_bin} {"not present" if random_bin not in bin_index_pair else "is empty"}')
                        break

                    isolate_lengths_selected_bin = bin_index_pair[random_bin].pop()
                    sampled.append(isolate_lengths_selected_bin)
                    original_bin.append(random_bin)
                    random_bins_counts[random_bin] -= 1
        # second round, these random_bins don't have a corresponding isolate bin
        # First search for the closest higher bin, if not there, take the closest lower bin
        for random_bin in sorted(random_bins_counts,reverse=True):
            updated_random_bin = None  # Initialize outside the loop
            while random_bins_counts[random_bin] > 0:
                if updated_random_bin is None:
                    updated_random_bin = find_nearest(bin_index_pair,random_bin)
                    if updated_random_bin < 10:
                        sys.exit("DEBUG")
                selected_bin_index_pair = bin_index_pair[updated_random_bin]
                if selected_bin_index_pair:
                    isolate_lengths_selected_bin = bin_index_pair[updated_random_bin].pop()
                    sampled.append(isolate_lengths_selected_bin)
                    original_bin.append(random_bin)
                    random_bins_counts[random_bin] -= 1
                else:
                    updated_random_bin = None  # Reset if bin_pair is empty

        df_selected = df_to_select.loc[sampled].reset_index().rename(columns={'count': 'old_index'})
        df_selected['original_bin'] = original_bin

        def get_start_coordinate(row):
            if row['original_bin'] >= row['bin']:
                return np.nan
            else:
                random_generator = np.random.default_rng(seed=random_number)
                return random_generator.integers(0, (row['length'] - (row['original_bin'] * 100 + 100)) + 1)

        df_selected['start'] = df_selected.apply(get_start_coordinate, axis=1)
        # below line does not work, because randint is executed for all, resulting in high <=0 for some rows
        # df_selected['start'] = np.where(df_selected['original_bin'] >= df_selected['bin'], 0, np.random.randint(0, high=(df_selected['length'] - (df_selected['original_bin'] * 100 + 100)) + 1))
        random_generator = np.random.default_rng(seed=random_number)
        df_selected['end'] = np.where(df_selected['original_bin'] >= df_selected['bin'],
                                      np.nan,
                                      df_selected['start'] + random_generator.integers(
                                          df_selected['original_bin'] * 100, high=df_selected['original_bin'] * 100 + 100 + 1))
        df_selected['start'] = df_selected['start'].astype('Int64')
        df_selected['end'] = df_selected['end'].astype('Int64')
        df_selected['final_length'] = df_selected['end'] - df_selected['start']

        df_selected.loc[:, ['id', 'start', 'end']].to_csv(output.subsampeled_isolate, sep='\t', header=False, index=False)
        df_selected.to_csv(output.subsampeled_isolate_w_bins,sep='\t',header=False,index=False)

        fig = px.histogram(df_selected, x='bin',
            category_orders={'bin': sorted(df_selected.loc[:, 'bin'].unique())})

        fig.update_xaxes(type='category')

        fig.write_html(output.subsampeled_isolate_barplot)

        duplicates_bin = df_selected.loc[:, ['id', 'bin']].value_counts().reset_index().rename(columns={'count': 'copies'})

        fig = px.histogram(duplicates_bin, x='bin',
            color="copies",
            category_orders={'copies': sorted(duplicates_bin['copies'].unique()),
                             'bin': sorted(duplicates_bin.loc[:, 'bin'].unique())})

        fig.update_xaxes(type='category')

        fig.write_html(output.subsampeled_isolate_barplot_copies)


        df_selected_difference = df_selected[df_selected['original_bin'] != df_selected['bin']]
        # TODO: write this cleaner
        try:
            max_bin = max(df_selected_difference['original_bin'])
        except:
            max_bin = 0

        fig = go.Figure(data=
        go.Parcoords(dimensions=list([
            dict(range=[0, max_bin], label='original_bin',values=df_selected_difference['original_bin']),
            dict(range=[0, max_bin],
                label='bin',values=df_selected_difference['bin'])
        ])
        )
        )

        fig.write_html(output.subsampeled_isolate_lineplot)

rule extract_selected_seqids_from_isolates:
    input:
        isolate_fasta = lambda wildcards: config["location_isolates"][wildcards.isolate_sample],
        selected_seqids = rules.sample_length_from_isolate_based_on_training_set.output.subsampeled_isolate
    output:
        subsampeled_isolate_fasta = Path('{ROOT}') / '{mix_n}' / 'subsample' / 'model' / '{isolate_sample}.fq'
    shell:
        """
        #!TODO
        ml seqtk/1.4;
        ml seqkit/2.3.1;
        seqtk subseq {input.isolate_fasta:q} {input.selected_seqids:q} | seqkit rename > {output.subsampeled_isolate_fasta:q}
        # ml seqkit/2.3.1;
        # seqkit grep -f <(cut -f 1 {input.selected_seqids}) {input.isolate_fasta} -o {output.subsampeled_isolate_fasta};
        """
