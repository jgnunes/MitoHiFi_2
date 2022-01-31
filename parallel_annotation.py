import concurrent.futures
import sys
import functools
import subprocess
import getMitoLength 

def annotate_mito(max_contig_size, related_gbk, gen_code, contig_id):
    #mitofinder = ["/software/team311/jf18/MitoFinder/mitofinder", "--max-contig-size", str(max_contig_size), "-j", contig_id + ".annotation", "-a", contig_id + ".mitogenome.fa", "-r", related_gbk, "-o", gen_code, "-p", str(threads)]
    print(f"Annotating contig {contig_id}...")
    mitofinder = ["mitofinder", "--max-contig-size", str(max_contig_size), "-j", contig_id + ".annotation", "-a", contig_id + ".mitogenome.fa", "-r", related_gbk, "-o", gen_code, "-p", "1"]
    subprocess.run(mitofinder)

if __name__ == "__main__":
    
    if sys.argv[1] == "-h":
        print("""Usage:
        arg1 = contigs
        arg2 = related_mito_fasta
        arg3 = related_mito_gbk
        arg4 = genetic code
        """)
        sys.exit()

    contigs = sys.argv[1].split()
    print(f"Contigs: {contigs}")


    rel_mito_len = getMitoLength.get_mito_length(sys.argv[2])
    print("Length of related mitogenome is: {} bp".format(rel_mito_len))

    # calculate maximum contig size accepted by mitofinder when annotating the contigs
    max_contig_size = 5*rel_mito_len

    related_gbk = sys.argv[3]
    gen_code = sys.argv[4]

    partial_annotate_mito = functools.partial(annotate_mito, max_contig_size, related_gbk, gen_code)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_annotate_mito, contigs)
    print("Finished annotation!")
