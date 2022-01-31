import subprocess
import numpy
from Bio import SeqIO
import sys #to remove after debugging


def get_largest_cluster(cdhit_clstr_file):
    #clusters = {}
    largest_cluster_len = 0
    largest_cluster = ""
    curr_sequences = []
    curr_cluster = []
    with open(cdhit_clstr_file, "r") as f:
        for line in f:
            if line[0] == ">":
                if len(curr_sequences) > largest_cluster_len:
                    largest_cluster = curr_cluster
                    largest_cluster_len = len(curr_sequences)
                    largest_cluster_seqs = curr_sequences
                curr_cluster = line.strip().replace(">","")
                curr_sequences = []
                #clusters[cluster_id] = 0
            else:
                curr_sequences.append(line.strip())
        # catch the last cluster        
        if len(curr_sequences) > largest_cluster_len:
            largest_cluster = curr_cluster
            largest_cluster_len = len(curr_sequences)
            largest_cluster_seqs = curr_sequences
    
    for sequence in largest_cluster_seqs:
        if sequence[-1] == "*":
            representative_seq = sequence

    return (largest_cluster, representative_seq)

def get_repr_contig(contigs_fasta, threads="1"):
    c_threshold = "0.8"
    wordsize = "4"
    cdhit_out = "cdhit.out"
    cdhit_out_clstr = cdhit_out + ".clstr"
    cdhit_cmd = ["cd-hit-est", "-i", contigs_fasta, "-c", c_threshold, "-n", wordsize, "-o", cdhit_out, "-T", str(threads)]
    subprocess.run(cdhit_cmd, shell=False)

    repr_contig_cluster, repr_contig_info = get_largest_cluster(cdhit_out_clstr)
    repr_contig_id = repr_contig_info.split()[2].replace(">","").split("_")[0]
    return (repr_contig_id, repr_contig_cluster)

if __name__ == "__main__":
    repr_contig_id, repr_contig_cluster = get_repr_contig(sys.argv[1], sys.argv[2])
    print("Representative contig is {} that belongs to {}".format(repr_contig_id, repr_contig_cluster))
