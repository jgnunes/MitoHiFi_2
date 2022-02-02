import subprocess
import numpy
from Bio import SeqIO
import sys #to remove after debugging


def get_largest_cluster(cdhit_clstr_file):
    """Reads CDHIT cluster file and returns the cluster with highest number of seqs.

    Args:
        cdhit_clst_file (str): file containing information about CDHIT clusters (*.clstr)

    Returns:
        tuple: Largest cluster ID, Representative sequence from the largest cluster
    """

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
    """Gets representative contig from a multifasta file

    Args:
        contigs_fasta (str): file containing all sequences 
        threads (str): number of threads to be used when running CDHIT

    Returns:
        tuple: Representative contig ID, CDHIT cluster where the representative contig came from
    """

    c_threshold = "0.8"
    wordsize = "4"
    cdhit_out = "cdhit.out"
    cdhit_out_clstr = cdhit_out + ".clstr"
    cdhit_cmd = ["cd-hit-est", "-i", contigs_fasta, "-d", "0", "-c", c_threshold, "-n", wordsize, "-o", cdhit_out, "-T", str(threads), "-M", "0"]
    subprocess.run(cdhit_cmd, shell=False, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    repr_contig_cluster, repr_contig_info = get_largest_cluster(cdhit_out_clstr)
    if "_rc_" in repr_contig_info:
        repr_contig_id = repr_contig_info.split()[2].replace(">","").split("_rc_rotated")[0]
    else:
        repr_contig_id = repr_contig_info.split()[2].replace(">","").split("_rotated")[0]
    return (repr_contig_id, repr_contig_cluster)

if __name__ == "__main__":
    repr_contig_id, repr_contig_cluster = get_repr_contig(sys.argv[1], sys.argv[2])
    print("Representative contig is {} that belongs to {}".format(repr_contig_id, repr_contig_cluster))
