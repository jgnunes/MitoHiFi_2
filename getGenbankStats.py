from Bio import SeqIO
import re
import sys

def get_contig_stats(contig_gb):
    frameshift_genes = set()
    all_genes = [] # defines all_genes as list because we want to detect duplicated genes
    for record in SeqIO.parse(contig_gb, "gb"):
        contig_len = len(record.seq) # calculates the contig length
        contig_num_genes = len(record.features)
        for feature in record.features:
            if feature.type == 'gene':
                gene_name = feature.qualifiers['gene'][0]
                all_genes.append(gene_name)
            if feature.type == 'CDS':
                gene_name = feature.qualifiers['gene'][0]
                gene_amino = feature.qualifiers['translation'][0]
                if re.search("[A-Z]\*[A-Z]", gene_amino):
                    print("Gene {} contains frameshift".format(gene_name))
                    frameshift_genes.add(gene_name)
    return (contig_len, contig_num_genes, all_genes, frameshift_genes) 

def main():
    stats = get_contig_stats(sys.argv[1])
    print("stats: ", stats)

if __name__ == "__main__":
    main()
