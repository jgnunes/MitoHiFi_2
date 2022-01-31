import subprocess
import pandas as pd
import argparse
from datetime import date
import parse_blast
from Bio import SeqIO
import sys
import os
import filterfasta
from circularizationCheck import circularizationCheck

def get_contigs_ids(blast_output):
    """
    args:
    blast_line is a blast output, that can be either parsed_blast.txt or parsed_blast_all.txt
    
    returns:
    The ID from each contig from blast_output, i.e., the BLAST queries
    """
    contigs_ids = set()

    num_lines = sum(1 for line in open(blast_output, "r") if line.strip()) # counts the number of lines in the input file
    if num_lines >= 2:
        with open(blast_output, "r") as f:
            next(f) # skips blast header
            for line in f:
                contigs_ids.add(line.split()[0])    

    return contigs_ids

def get_circo_mito(contig_id, circular_size, circular_offset):
    """
    It gets a contig ID and circularizes it, registering the circularization points
    args:
    the ID from the contig that we want to circularize
    """

    def cut_coords(contig_fasta, circularization_position, fasta_out):

        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        get= record[circularization_position:]
        with open(mitogenome_fasta_out, "w") as f:
            f.write(get.format('fasta'))

    # find circularization point
    circularization_history = []
    contig_fasta = "".join([contig_id, ".mito.fa"])
    mitogenome_fasta_out = "".join([contig_id, ".mitogenome.fa"])
    circularization_checks = "".join([contig_id, ".circularisationCheck.txt"])

    circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset)

    # writes circularization information to '[contig_id].circularisationCheck.txt' file
    with open(circularization_checks, "w") as f:
        f.write(str(circularization_info))
        circularization_history.append(str(circularization_info))
    
    is_circularizable = circularization_info[0]
    circularization_position = int(circularization_info[2])

    # if the contig is not circularizable, then create the "mitogenome.fasta" directly from it
    if not is_circularizable:
        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        with open(mitogenome_fasta_out, "w") as f:  
            f.write(record.format('fasta'))        
    else:
        # if the contig is circularizable, then run circularization iteratively until the "mitogenome.fasta"
        # is no longer circularizable    
        while is_circularizable:
            cut_coords(contig_fasta, circularization_position, mitogenome_fasta_out)
            contig_fasta = mitogenome_fasta_out
            circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset)
            is_circularizable = circularization_info[0]
            circularization_position = int(circularization_info[2])
            with open(circularization_checks, "a") as f:
                f.write("\n" + str(circularization_info))
                circularization_history.append(str(circularization_info))

    return circularization_history

def main():

    today = date.today()

    parser= argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "blablabli") 
    parser.add_argument("-r", help= "-r: Pacbio Hifi Reads from your species", required = "True")
    parser.add_argument("-f", help= "-f: Close-related Mitogenome is fasta format", required = "True")
    parser.add_argument("-g", help= "-k: Close-related species Mitogenome in genebank format", required = "True")
    parser.add_argument("-t", help= "-t: Number of threads for (i) hifiasm and (ii) the blast search", required = "True", type=int)
    parser.add_argument("-l", help="-l: Length of the close-related mitogenome (to filter out possible NUMTs)", required = "True", type=int)
    parser.add_argument("-p", help="-p: Percentage of query in the blast match with close-related mito", type=int, default=50)
    parser.add_argument('--circular-size', help='Size to consider when checking for circularization', type=int, default=220)
    parser.add_argument('--circular-offset', help='Offset from start and finish to consider when looking for circularization', type=int, default=40)
    parser.add_argument("-o", help="""-o: Organism genetic code following NCBI table (for mitogenome annotation):
    1. The Standard Code 2. The Vertebrate MitochondrialCode 3. The Yeast Mitochondrial Code 
    4. The Mold,Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code 5. The Invertebrate Mitochondrial Code 
    6. The Ciliate, Dasycladacean and Hexamita Nuclear Code 9. The Echinoderm and Flatworm Mitochondrial Code 10. The Euplotid Nuclear Code 
    11. The Bacterial, Archaeal and Plant Plastid Code 12. The Alternative Yeast Nuclear Code 13. The Ascidian Mitochondrial Code 
    14. The Alternative Flatworm Mitochondrial Code 16. Chlorophycean Mitochondrial Code 21. Trematode Mitochondrial Code 
    22. Scenedesmus obliquus Mitochondrial Code 23. Thraustochytrium Mitochondrial Code 24. Pterobranchia Mitochondrial Code 
    25. Candidate Division SR1 and Gracilibacteria Code 
        """, type=str, default='1')
    args = parser.parse_args()

    print("MitoHifi v2.2" + "/n")
    print("Started at:" , today)
    print("\n" + "1-)First we map your PacbioHiFi reads to the close-related mitogenome+ ""\n")

    minimap = "minimap2 "+ "-t " + str(args.t) + " " + "--secondary=no -ax map-pb " + args.f + " " + args.r + " | " + "samtools view -@ " + str(args.t) + "-S -b -F4 -F 0x800 > reads.HiFiMapped.bam"
    subprocess.run(minimap, shell=True)

    print("\n" + "2-)Now we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS." + "\n")
    samtools = "samtools fasta " + "reads.HiFiMapped.bam >  gbk.HiFiMapped.bam.fasta"
    print(samtools)
    subprocess.run(samtools, shell=True)

    filterfasta.filterFasta(minLength=args.l, neg=True, inStream="gbk.HiFiMapped.bam.fasta", outPath="gbk.HiFiMapped.bam.filtered.fasta")

    print("\n" + "3-)Now let's run hifiasm to assemble the mapped and filtered reads!" + "\n")

    hifiasm = "hifiasm -t " + str(args.t) + " -o gbk.HiFiMapped.bam.filtered.assembled gbk.HiFiMapped.bam.filtered.fasta 2>hifiasm.log"
    print(hifiasm)
    subprocess.run(hifiasm, shell=True)

    get_p = "scripts/gfa2fa gbk.HiFiMapped.bam.filtered.assembled.p_ctg.gfa >  gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa"
    get_a = "scripts/gfa2fa gbk.HiFiMapped.bam.filtered.assembled.a_ctg.gfa >  gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa"
    cat = "cat gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa > hifiasm.contigs.fasta"

    subprocess.run(get_p, shell=True)
    subprocess.run(get_a, shell=True)
    subprocess.run(cat, shell=True)

    print("\n" + "4-)Now let's run the blast of the assembled contigs with the close-related mitogenome" + "\n")

    makeblastdb = "makeblastdb -in " + args.f + " -dbtype nucl"
    print(makeblastdb)
    subprocess.run(makeblastdb, shell=True)
    blast = "blastn -query hifiasm.contigs.fasta -db " + args.f + " -num_threads " + str(args.t) + " -out contigs.blastn -outfmt '6 std qlen slen'"
    print(blast)
    subprocess.run(blast, shell=True)
    print("Blast done!")

    #the next script parses a series of conditions to exclude blast with NUMTs. 
    parse_blast.parse_blast(query_perc=args.p)

    #We check for circularisation

    # select contigs to be circularized
    # first look for contigs in parsed_blast.txt
    contigs_ids = get_contigs_ids("parsed_blast.txt")

    # if we don't find contigs in parse_blast.txt 
    # look for contigs in parsed_blast_all.txt
    if len(contigs_ids) == 0:
        contigs_ids = get_contigs_ids("parsed_blast_all.txt")

    # if we can't find any contigs even in parsed_blast_all.txt, then we exit the pipeline
    if len(contigs_ids) == 0:
        sys.exit("""\n Attention! \n The 'parsed_blast.txt' and 'parsed_blast_all.txt' files are empty. The pipeline has stopped !! \n You need to run further scripts to check if you have mito reads pulled to a large NUMT!!""")

    # records all contigs kept for the downstream steps in a file called 'contigs_ids.txt'
    with open("contigs_ids.txt", "w") as f:
        for contig_id in contigs_ids:
            f.write(contig_id + "\n")

    # removes file that contains history of circularization of it already exists
    try:
        os.remove('all_contigs.circularisationCheck.txt')
    except OSError:
        pass

    #retrieves the FASTA files for each contig
    for contig_id in contigs_ids:
        # retrieves the FASTA files for each contig
        filterfasta.filterFasta(idList=[contig_id], inStream="hifiasm.contigs.fasta", outPath="".join([contig_id, ".mito.fa"]))
        # circularizes each contig and saves circularization history to a file
        circularization_history = get_circo_mito(contig_id, args.circular_size, args.circular_offset)
        for circularization_event in circularization_history: 
            with open('all_contigs.circularisationCheck.txt', 'a') as f:
                f.write("\t".join([contig_id, circularization_event, "\n"]))    
        # annotates mitogenome(s) using mitofinder
        #annotate_mito = "mitofinder -j {} -a {}.mitogenome.fa -r {} -o {}".format(contig_id, contig_id, args.g, args.o)

if __name__ == '__main__':
    main()
