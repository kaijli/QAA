#!/usr/bin/env python
'''
python script to plot trimmed read length distributions for two fastq files on the same plot. 
'''
import argparse
import numpy as np
import matplotlib.pyplot as plt
import gzip

# def get_args():
# 	parser = argparse.ArgumentParser(description="")
# 	parser.add_argument("-f", "--filename", help="filename", type=str, required=True)
# 	return parser.parse_args()

# args = get_args()
# filename = args.filename

filename = "/projects/bgmp/kli8/bioinformatics/Bi623/QAA/test.fastq"
read1_file = "/projects/bgmp/kli8/bioinformatics/Bi623/QAA/trimmomatic_results/paired_fwd_4_2C_mbnl_S4_L008_R1.fastq.gz"
read2_file = "/projects/bgmp/kli8/bioinformatics/Bi623/QAA/trimmomatic_results/paired_rev_4_2C_mbnl_S4_L008_R2.fastq.gz"

def count_len(filename: str, size: int) -> np.ndarray:
    '''
    takes a fq file and the number of records. 
    returns the length of each read in an array. 
    '''
    arr = np.zeros(int(size), dtype=np.float64)
    with gzip.open(filename, "r") as fh:
        file_line = 0
        read_num = 0
        for line in fh:
            line = line.strip()
            # print(f"{line=}")
            if file_line % 4 == 1:
                # print(f"{read_num=}")
                arr[read_num] = len(line)
                read_num += 1
            file_line += 1
    return arr

def len_freq(arr: np.ndarray) -> dict:
    counts = dict()
    sorted_counts = dict()
    for item in arr:
        if item not in counts.keys():
            counts[item] = 1
        else:
            counts[item] += 1
    for key in sorted(counts.keys()):
        sorted_counts[key] = counts[key]        
    return sorted_counts

def plot_dists(arr1: np.ndarray, arr2: np.ndarray, suffix: str):
    '''
    takes an array and plots a distribution plot
    '''
    fig, ax = plt.subplots()
    ax.hist(arr1, facecolor='g')
    ax.hist(arr2, facecolor='b')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.yscale("log")  
    plt.title(f"Proportion of Read Lengths\n{filename}")
    plt.savefig(f"hist_{suffix}.png")

def plot_dicts(dict1: dict, dict2: dict, suffix: str):
    '''
    takes an array and plots a distribution plot
    '''
    fig, ax = plt.subplots()
    x1 = list(dict1.keys())
    y1 = list(dict1.values())
    x2 = list(dict2.keys())
    y2 = list(dict2.values())
    ax.plot(x1, y1, color='g', label='Read 1')
    ax.plot(x2, y2, color='b', label='Read 2')
    plt.legend()
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.yscale("log")  
    plt.title(f"Proportion of Read Lengths\nPost Quality and Adapter Trimming")
    plt.savefig(f"hist_{suffix}.png")

# tester = count_len(filename, 100)
# print(f"{tester=}")
# plot_dists(tester, 'test')

r1 = count_len(read1_file, 8980380)
r2 = count_len(read2_file, 8980380)
r1_dict = len_freq(r1)
r2_dict = len_freq(r2)
plot_dicts(r1_dict, r2_dict, 'r1r2_dict')
# plot_dists(r1, r2, 'r1r2')