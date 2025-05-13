from pathlib import Path



rule filter_fastq:
    input:
        lambda wildcards: config['samples'][wildcards.sample_name]['input']
    output:
        filtered_fastq = Path('{ROOT}') / '{sample_name}' / 'filtered_fastq' / 'filtered.fq'
    params:
        length=lambda wildcards: config['filtering']['length'],
        quality=lambda wildcards: config['filtering']['quality']
    threads: 1
    shell:
        """
        ml seqkit/2.3.1; \
        seqkit seq \
        -m {params.length} \
        -Q {params.quality} \
        {input} > {output}
        """

rule kma_mapping:
    input:
        # lambda wildcards: config['samples'][wildcards.sample_name]['input']
        rules.filter_fastq.output.filtered_fastq
    output:
        mapstat = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma.mapstat',
        tsv = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma.tsv'
    params:
        DB = config['kma_db'],
        prefix = lambda wildcards: Path(wildcards.ROOT) / wildcards.sample_name / 'classification' / 'kma',
        MP = 20,
        MRS = 0.0,
        BC = 0.7,
        TMP = lambda wildcards: f"{Path(wildcards.ROOT) / wildcards.sample_name / 'classification'}/"
    threads: 4
    shell:
        """
        mkdir -p {params.TMP};
        ml kma/1.4.12a;
        -i {input} \
        -o {params.prefix} \
        -t_db {params.DB} \
        -t {threads} \
        -mrs {params.MRS} \
        -bcNano \
        -bc {params.BC} \
        -ef \
        -a \
        -mem_mode \
        -tmp {params.TMP} \
        -1t1 \
        -matrix \
        -tsv \
        -shm
        """

rule merge_kma_tsv_and_mapstat:
    input:
        tsv = rules.kma_mapping.output.tsv,
        mapstat = rules.kma_mapping.output.mapstat
    output:
        tsv_mapstat = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma_tsv_mapstat.tsv'
    run:
        import pandas as pd

        tsv = pd.read_table(input.tsv, usecols=['Template_Name', 'Template_Length', 'Template_Identity',
       'Template_Coverage', 'Template_Depth', 'Query_Identity',
       'Query_Coverage', 'Query_Depth', 'Read_Count_Map', 'Read_Count_Aln'])

        mapstat = pd.read_table(input.mapstat, skiprows=6)
        mapstat.rename(columns={'# refSequence': 'Template_Name'}, inplace=True)

        tsv_mapstat = tsv.merge(mapstat.loc[:, ['Template_Name', "refConsensusSum", "bpTotal"]],
                                on = 'Template_Name')

        tsv_mapstat.to_csv(output.tsv_mapstat,sep='\t',header=True,index=False)

rule remove_plasmids_from_kma:
    input:
        tsv_mapstat = rules.merge_kma_tsv_and_mapstat.output.tsv_mapstat,
    output:
        tsv_mapstat_wo_plasmids = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma_tsv_mapstat_wo_pl.tsv',
    params: lookup=config['lookup']
    run:
        import pandas as pd

        tsv_mapstat = pd.read_table(input.tsv_mapstat)
        df_lookup = pd.read_table(params.lookup,
                                    names=['full_name', 'length', 'accession', 'seq_type', 'taxid',
                                           'species_name', 'species_taxid', "best_accession"])

        plasmid_mask = tsv_mapstat['Template_Name'].map(df_lookup.set_index('full_name')['seq_type']).eq('plasmid')
        tsv_mapstat = tsv_mapstat.loc[~plasmid_mask, :]
        tsv_mapstat.to_csv(output.tsv_mapstat_wo_plasmids, sep='\t', header=True, index=False)

rule collapse_kma_tsv_mapstat:
    input:
        tsv_mapstat = rules.remove_plasmids_from_kma.output.tsv_mapstat_wo_plasmids,
    output:
        tsv_mapstat = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma_tsv_mapstat_wo_pl_collapsed.tsv',
    params: lookup=config['lookup']
    run:
        import pandas as pd

        lookup = {}
        with open(params.lookup) as f:
            for line in f:
                full_name, length, assembly_accession, seq_type, taxid, species_name, species_taxid, _ = line.strip().split('\t')
                if assembly_accession not in lookup:
                    lookup[assembly_accession] = {"total_length": 0,
                                         "total_length_wo_pl_length": 0,
                                         "pl_length": 0,
                                         "species_name": species_name}

                length = int(length)
                lookup[assembly_accession]["total_length"] += length
                if seq_type == 'plasmid':
                    lookup[assembly_accession]["pl_length"] += length
                else:
                    lookup[assembly_accession]["total_length_wo_pl_length"] += length


        tsv_mapstat = pd.read_table(input.tsv_mapstat)
        tsv_mapstat.insert(1, "genome", tsv_mapstat["Template_Name"].str.split(" ").str[1].str.split("__").str[0])

        tsv_mapstat_collapsed = tsv_mapstat.groupby("genome", sort=True).apply(
            lambda x: pd.Series({"species_name": lookup.get(x.name).get('species_name'),
                                 "Genome_length": lookup.get(x.name).get("total_length"),
                                 "Read_Count_Aln": sum(x["Read_Count_Aln"]),
                                 "bpTotal": sum(x["bpTotal"]),
                                 "Template_Identity": 100 * sum(x["refConsensusSum"]) / lookup.get(x.name).get("total_length_wo_pl_length"),
                                 "Depth": sum(x["bpTotal"]) / lookup.get(x.name).get("total_length_wo_pl_length")
                                 }
                                )).reset_index()

        tsv_mapstat_collapsed.to_csv(output.tsv_mapstat, sep='\t', header=True, index=False)

rule highest_id_per_species:
    input:
        rules.collapse_kma_tsv_mapstat.output.tsv_mapstat
    output:
        tsv_mapstat_highest_id = Path('{ROOT}') / '{sample_name}' / 'classification' / 'kma_tsv_mapstat_wo_pl_collapsed_highest_id.tsv',
    run:
        import pandas as pd

        tsv_mapstat = pd.read_table(input[0])

        tsv_mapstat_highest_id = tsv_mapstat.iloc[tsv_mapstat.groupby('species_name', sort=False)['Template_Identity'].idxmax(), :]

        tsv_mapstat_highest_id.to_csv(output.tsv_mapstat_highest_id, sep='\t', header=True, index=False)


