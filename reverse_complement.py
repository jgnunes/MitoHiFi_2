import sys
from Bio import SeqIO

def reverse_complement(in_gb, out_gb):
    
    # open original record and convert it to 
    # its reverse complement (rc_record)
    record = SeqIO.read(in_gb, "genbank")
    rc_record = record.reverse_complement(id=record.id + "_rc")
    # needs to set molecule type as DNA because Biopython deletes
    # this info when parsing
    rc_record.annotations["molecule_type"] = "DNA"
    
    with open(out_gb, "w") as f:
        SeqIO.write(rc_record, f, "genbank")
    
    out_fasta = out_gb.replace(".gb", ".fasta")
    with open(out_fasta, "w") as f:
        SeqIO.write(rc_record, f, "fasta")

if __name__ == "__main__":
    in_gb = sys.argv[1]
    out_gb = sys.argv[2]
    reverse_complement(in_gb, out_gb)
