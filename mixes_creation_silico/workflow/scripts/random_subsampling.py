import random
import argparse
import numpy as np
from collections import defaultdict


def read_tsv(file_path):
    input_data = {}
    with open(file_path) as f:
        for idx, line in enumerate(f):
            id_read, length = line.rstrip().split('\t')
            input_data[idx] = (id_read, int(length))
    return input_data


def subsample_until_sum_without_replacement(input_data, target_sum, seed):
    total_length = 0
    subsampled_read_id = defaultdict(int)
    idxs = list(input_data.keys())
    rng = np.random.default_rng(seed)
    idx_shuffeled = rng.permutation(idxs)
    for idx in idx_shuffeled:
        read_id, read_length = input_data[idx]
        subsampled_read_id[read_id] += 1
        total_length += read_length
        if total_length >= target_sum:
            break
    return subsampled_read_id


def subsample_until_sum_with_replacement(input_data, target_sum, seed):
    total_length = 0
    subsampled_read_id = defaultdict(int)
    idx_max = len(input_data) - 1
    random.seed(seed)
    while total_length < target_sum:
        idx_random = random.randint(0, idx_max)
        read_id, read_length = input_data[idx_random]
        subsampled_read_id[read_id] += 1
        total_length += read_length
    return subsampled_read_id


def subsample_until_sum(input_data, target_sum, seed):
    total_bases_input = sum([length for read_id, length in input_data.values()])
    if target_sum <= total_bases_input:
        print(f"Requested bases ({target_sum}) < Input bases ({total_bases_input}) -- Sampling WITHOUT replacement")
        subsampled_reads = subsample_until_sum_without_replacement(input_data, target_sum, seed)
    else:
        print(f"Requested bases ({target_sum}) > Input bases ({total_bases_input}) -- Sampling WITH replacement")
        subsampled_reads = subsample_until_sum_with_replacement(input_data, target_sum, seed)
    return subsampled_reads


def save_txt(subsampled_reads, output_path):
    with open(output_path, 'w') as f:
        for read_id, count in subsampled_reads.items():
            for _ in range(count):
                f.write(f"{read_id}\n")


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Subsample TSV file until a sum of lengths is reached")
    # parser.add_argument('input_file', type=str, help="Input TSV file path")
    # parser.add_argument('output_file', type=str, help="Output TSV file path")
    # parser.add_argument('target_sum', type=int, help="Target sum of lengths")
    #
    # args = parser.parse_args()

    SEED = snakemake.params[2]
    target_sum = round(snakemake.params[0] * (snakemake.params[1] / 100))
    input_file=snakemake.input[0]
    output_file=snakemake.output[0]


    input_data = read_tsv(input_file)
    subsample = subsample_until_sum(input_data, target_sum, SEED)
    save_txt(subsample, output_file)
