import concurrent.futures
from concurrent.futures import wait
import logging
import shutil
import subprocess
import time
import warnings
import pandas as pd
import argparse
import parse_blast
from Bio import SeqIO
import sys
import os
import cleanUpCWD
import filterfasta
import findFrameShifts
import fixContigHeaders
import functools
import rotation
import getMitoLength
import getReprContig
import shlex
from circularizationCheck import circularizationCheck
import alignContigs

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
    logging.debug("Defining cut_coords function...\n") # for debugging
    def cut_coords(contig_fasta, circularization_position, fasta_out):
        print(f"cut_coords function called!\ncontig_fasta: {contig_fasta}; circularization_position: {circularization_position}; fasta_out: {fasta_out}") # for debugging
        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        get= record[circularization_position:]
        with open(fasta_out, "w") as f:
            f.write(get.format('fasta'))

    # find circularization point
    circularization_history = []
    contig_fasta = "".join([contig_id, ".mito.fa"])
    mitogenome_fasta_out = "".join([contig_id, ".mitogenome.fa"])
    circularization_checks = "".join([contig_id, ".circularisationCheck.txt"])

    circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset, contig_id=contig_id)

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
            circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset, contig_id=contig_id)
            is_circularizable = circularization_info[0]
            circularization_position = int(circularization_info[2])
            with open(circularization_checks, "a") as f:
                f.write("\n" + str(circularization_info))
                circularization_history.append(str(circularization_info))
    
    print(f"circularization_history: {circularization_history}") # for debugging
    return circularization_history

def get_most_freq_tRNA():
    tRNAs = {}
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.trnas'):
            with open(curr_file, "r") as infile:
                for line in infile:
                    tRNA = line.split("\t")[0]
                    if tRNA not in tRNAs:
                        tRNAs[tRNA] = 1
                    else:
                        tRNAs[tRNA] += 1
    print(f"tRNAS: {tRNAs}") #debug
    # if any contig has a tRNA-Phe, use it as the reference gene for rotation
    if 'tRNA-Phe' in tRNAs:
        most_freq_tRNA = 'tRNA-Phe'
    else:
        most_freq_tRNA = max(tRNAs, key=tRNAs.get)
    return most_freq_tRNA
    
def process_contig(threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, contig_id): 
    print(f"\nWorking with contig {contig_id}. Annotation will be done with {threads_per_contig} threads\n")
    print(f"contigs: {contigs} | args.circular_size: {circular_size}") # for debugging
    # retrieves the FASTA files for each contig
    filterfasta.filterFasta(idList=[contig_id], inStream=contigs, outPath="".join([contig_id, ".mito.fa"]))
    # circularizes each contig and saves circularization history to a file
    circularization_history = get_circo_mito(contig_id, circular_size, circular_offset)
    #for circularization_event in circularization_history: 
        #with open('all_contigs.circularisationCheck.txt', 'a') as f:
            #f.write("\t".join([contig_id, circularization_event, "\n"])) 
    # annotates mitogenome(s) using mitofinder
    print("Running mitofinder with maximum contig size of {} bp".format(max_contig_size))
    subprocess.run(["mitofinder", "--new-genes", "--max-contig-size", str(max_contig_size), "-j", contig_id+".annotation", "-a", contig_id+".mitogenome.fa", "-r", rel_gbk, "-o", gen_code, "-p", str(threads_per_contig)], stderr=subprocess.STDOUT)
    # rotates the mitogenome
    mitogenome_gb = os.path.join(contig_id + ".annotation", contig_id + ".annotation_MitoFinder_mitfi_Final_Results", contig_id + ".annotation_mtDNA_contig.gb") 
    if not os.path.isfile(mitogenome_gb):
        warnings.warn("Contig "+ contig_id + " does not have an annotation file, check MitoFinder's output")
        return
    
    trnas = rotation.get_trna_pos(mitogenome_gb)
    print(f"trnas: {trnas}") #debug
    if trnas:
        with open(f"{contig_id}.trnas", "w") as outfile:
            for trna in trnas:
                outfile.write("\t".join([trna, trnas[trna][0], trnas[trna][1], "\n"]))
        return
    else:
        warnings.warn('No tRNA gene found in ' + mitogenome_gb + '... Skipping contig ' + contig_id + '\n')
        return


def process_contig_02(ref_tRNA, threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, contig_id): 
    print("Started process_contig_02...") #debug 
    if not os.path.isfile(f"{contig_id}.trnas"):
        warnings.warn(f"Contig {contig_id} does not have annotated tRNAs, skipping it...")
        return
    
    with open(f"{contig_id}.trnas", "r") as infile:
        for line in infile:
            if line.strip().split("\t")[0] == ref_tRNA:
                start = int(line.strip().split("\t")[1])
                strand = int(line.strip().split("\t")[2])
    
    if start == None:
        warnings.warn(f"Reference gene {ref_tRNA} is not present in contig {contig_id}. Skipping contig...")
        return

    print(f"start: {start} | strang: {strand}") #debug
    mitogenome_gb = os.path.join(contig_id + ".annotation", contig_id + ".annotation_MitoFinder_mitfi_Final_Results", contig_id + ".annotation_mtDNA_contig.gb")

    genome = contig_id + ".mitogenome.fa"
    new_gb = None
    if strand == -1:
        genome_rc = contig_id + "_RC.mitogenome.fa"
        rc = os.path.join(os.path.dirname(genome), genome_rc)
        rotation.make_rc(genome, rc)
        new_gb = rotation.annotate(os.path.dirname(genome), os.path.abspath(genome_rc), os.path.abspath(rel_gbk), contig_id, gen_code, max_contig_size, str(threads_per_contig))
        #start, strand = rotation.get_phe_pos(new_gb)
        genome = rc
    rotation.rotate(genome, start, contig_id)
    print(' '.join([f'Rotated to {ref_tRNA} genome is at', os.path.join(os.path.dirname(genome), contig_id + '.mitogenome.rotated.fa')]))
    if new_gb:
        print(f"new_gb: {new_gb}")
        print('Mitogenome annotation is at ', new_gb)
    # check frameshifts in genes from contig and append save findings to 
    # `{contig_id}.individual.stats` intermediate file
    frameshifts = findFrameShifts.find_frameshifts(mitogenome_gb)
    gb_len, num_genes = findFrameShifts.get_gb_stats(mitogenome_gb)
    contig_dir = os.path.join("potential_contigs", contig_id)
    mitogenome_location = os.path.join(contig_dir, mitogenome_gb)
    if not frameshifts:
        all_frameshifts = "No frameshift found"
    elif len(frameshifts)==1:
        all_frameshifts = "".join(frameshifts)
    elif len(frameshifts)>1:
        all_frameshifts = ";".join(frameshifts)
    with open(f"{contig_id}.individual.stats", "w") as outfile:
        outfile.write("\t".join([contig_id, all_frameshifts, mitogenome_location, gb_len, num_genes+"\n"]))

def main():
    
    start_time = time.time()

    parser= argparse.ArgumentParser(add_help=False)
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-r", help= "-r: Pacbio Hifi Reads from your species")
    mutually_exclusive_group.add_argument("-c", help= "-c: Assemnbled fasta contigs/scaffolds to be searched to find mitogenome")
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "Print this help message.")
    parser.add_argument("-f", help= "-f: Close-related Mitogenome is fasta format", required = "True")
    parser.add_argument("-d", help="-d: debug mode to output additional info on log", default=False)
    parser.add_argument("-g", help= "-k: Close-related species Mitogenome in genebank format", required = "True")
    parser.add_argument("-a", help="-a: Choose between animal (default) or plant", default="animal", choices=["animal", "plant"])
    parser.add_argument("-t", help= "-t: Number of threads for (i) hifiasm and (ii) the blast search", required = "True", type=int)
    parser.add_argument("-p", help="-p: Percentage of query in the blast match with close-related mito", type=int, default=50)
    parser.add_argument("-m", help="-m: Number of bits for HiFiasm bloom filter [it maps to -f in HiFiasm] (default = 0)", type=int, default=0)
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
    

    # Set log message format
    FORMAT='[%(asctime)s %(levelname)s] %(message)s'

    if args.d:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                            format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                            format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
    
    # Welcome message
    logging.info("Welcome to MitoHifi v2. Starting pipeline...")
    logging.debug("Running MitoHiFi on debug mode.")
    
    # Measure the length of the related mitogenome 
    rel_mito_len = getMitoLength.get_mito_length(args.f)
    rel_mito_num_genes = getMitoLength.get_mito_genes(args.g)
    logging.info("Length of related mitogenome is: {} bp".format(rel_mito_len))
    logging.info("Number of genes on related mitogenome: {}".format(rel_mito_num_genes))
    
    # Calculate maximum contig size accepted by mitofinder when annotating the contigs
    max_contig_size = 5*rel_mito_len

    # If input are reads, map them to the related mitogenome and assemble the mapped ones
    if args.r:
        logging.info("Running MitoHifi pipeline in reads mode\n")
       
        logging.info("First we map your PacbioHiFi reads to the close-related mitogenome")

        #print(shlex.split(args.r))
        minimap_cmd = ["minimap2", "-t", str(args.t), "--secondary=no", "-ax", "map-pb", args.f] + shlex.split(args.r) 
        samtools_cmd = ["samtools", "view", "-@", str(args.t), "-S", "-b", "-F4", "-F", "0x800"] 
        minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        mapped_reads_f = open("reads.HiFiMapped.bam", "w")
        subprocess.run(samtools_cmd, stderr=subprocess.STDOUT, stdin=minimap.stdout, stdout=mapped_reads_f)
        minimap.wait()
        minimap.stdout.close()

        logging.info("Now we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS")
        mapped_fasta_f = open("gbk.HiFiMapped.bam.fasta", "w")
        subprocess.run(["samtools", "fasta", "reads.HiFiMapped.bam"], stdout=mapped_fasta_f)

        filterfasta.filterFasta(minLength=rel_mito_len, neg=True, inStream="gbk.HiFiMapped.bam.fasta", outPath="gbk.HiFiMapped.bam.filtered.fasta")

        logging.info("Now let's run hifiasm to assemble the mapped and filtered reads!")
        
        with open("hifiasm.log", "w") as hifiasm_log_f:
            subprocess.run(["hifiasm", "-t", str(args.t), "-f", str(args.m), "-o", "gbk.HiFiMapped.bam.filtered.assembled", "gbk.HiFiMapped.bam.filtered.fasta", ], stderr=subprocess.STDOUT, stdout=hifiasm_log_f)
        
        gfa2fa_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),"gfa2fa") # gets path to gfa2fa script
        
        with open("gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa", "w") as p_ctg_f:
            subprocess.run([gfa2fa_script, "gbk.HiFiMapped.bam.filtered.assembled.p_ctg.gfa"], stdout=p_ctg_f)
        with open("gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa", "w") as a_ctg_f:
            subprocess.run([gfa2fa_script, "gbk.HiFiMapped.bam.filtered.assembled.a_ctg.gfa"], stdout=a_ctg_f)
        
        with open("hifiasm.contigs.fasta", "w") as hifiasm_f:
            subprocess.run(["cat", "gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa", "gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa"], stdout=hifiasm_f)
        
        contigs = "hifiasm.contigs.fasta"
    
    else:
        logging.info("Running MitoHifi pipeline in contigs mode")
       
        logging.info("Fixing potentially conflicting FASTA headers...")
        original_contigs = args.c
        fixContigHeaders.fix_headers(original_contigs, "fixed_header_contigs.fasta")
        
        os.remove(original_contigs) # remove original contig file  
        shutil.move("fixed_header_contigs.fasta", original_contigs) # replace original contigs file by the version that has the headers fixed
        
        contigs = original_contigs
        
    
    logging.info("Let's run the blast of the contigs versus the close-related mitogenome")

    makeblastdb = "makeblastdb -in " + args.f + " -dbtype nucl"
    logging.debug(makeblastdb)
    subprocess.run(["makeblastdb", "-in", args.f, "-dbtype", "nucl"], stderr=subprocess.STDOUT)
    logging.debug("makeblastdb done. Running blast with the contigs")
    subprocess.run(["blastn", "-query", contigs, "-db", args.f, "-num_threads", str(args.t), "-out", "contigs.blastn", "-outfmt", "6 std qlen slen"], stderr=subprocess.STDOUT)
    logging.debug("Blast done!" + "\n")

    #the next script parses a series of conditions to exclude blast with NUMTs. 
    if args.a == "plant":
        ## if species is a plant, set minimum query percentage equal to 0% of related mito 
        ## and maximum query lenght 10 times the lenght of the related
        parse_blast.parse_blast(query_perc=args.p, min_query_len=0, max_query_len=10)
    else:
        ## if species is an animal, set minimum query percentage equal to 80% of related mito
        ## and maximum query length 5 times the length of the related (default values from 
        ## parse_blast function
        parse_blast.parse_blast(query_perc=args.p)

    # select contigs to be circularized
    # first look for contigs in parsed_blast.txt
    contigs_ids = get_contigs_ids("parsed_blast.txt")

    # if we don't find contigs in parse_blast.txt 
    # look for contigs in parsed_blast_all.txt
    if len(contigs_ids) == 0:
        contigs_ids = get_contigs_ids("parsed_blast_all.txt")

    # if we can't find any contigs even in parsed_blast_all.txt, then we exit the pipeline
    if len(contigs_ids) == 0:
        sys.exit("""Attention!
'parsed_blast.txt' and 'parsed_blast_all.txt' files are empty.
The pipeline has stopped !! You need to run further scripts to check if you have mito reads pulled to a large NUMT!!""")

    # records all contigs kept for the downstream steps in a file called 'contigs_ids.txt'
    with open("contigs_ids.txt", "w") as f:
        for contig_id in contigs_ids:
            f.write(contig_id + "\n")

    # removes file that contains history of circularization of it already exists
    try:
        os.remove('all_contigs.circularisationCheck.txt')
    except OSError:
        pass
      
    logging.info("Now we are going to circularize, annotate and rotate each contig which is a potential mitogenome")
    
    # creates a dictionary that will save frameshifts information for each contig
    contig_shifts = {}

    threads_per_contig = 1
    if args.t // len(contigs_ids) > 1:
        threads_per_contig = args.t // len(contigs_ids)

    partial_process_contig = functools.partial(process_contig, threads_per_contig, args.circular_size, args.circular_offset, contigs, max_contig_size, args.g, args.o)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_process_contig, contigs_ids)
    
    tRNA_ref = get_most_freq_tRNA() 
    logging.debug(f"tRNA_ref: {tRNA_ref}") #debug 
    
    
    partial_process_contig_02 = functools.partial(process_contig_02, tRNA_ref, threads_per_contig, args.circular_size, args.circular_offset, contigs, max_contig_size, args.g, args.o)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_process_contig_02, contigs_ids)

    # concatenates circularization histories from all contigs into 
    # a single concatenate history file 
    circularization_history_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.circularisationCheck.txt'):
            circularization_history_files.append(curr_file)
    logging.debug(f"circularization_history_files: {circularization_history_files}") # for debugging
    with open("all_contigs.circularisationCheck.txt", "a") as outfile:
        for circ_hist_file in circularization_history_files:
            contig_id = circ_hist_file.split('.')[0]
            with open(circ_hist_file, "r") as infile:
                for line in infile:
                    outfile.write("\t".join([contig_id, line.strip()+"\n"]))
    
    logging.debug(f"contig_shifts after parallel processing: {contig_shifts}") # for debugging 
    #align final mitogenome rotated contigs
    logging.info("Now the final rotated contigs will be aligned")
    # list all final rotated contigs 
    contigs_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('mitogenome.rotated.fa'):
            contigs_files.append(curr_file)
    # first concatenate all rotated contigs into a single multifasta file
    logging.info("List of contigs that will be aligned: " + str(contigs_files) + "\n")
    concat_fasta = alignContigs.concatenate_contigs(contigs_files)
    # then run MAFFT alignment between the rotated contigs using the multifasta as input and clustal as output format
    alignContigs.mafft_align(multifasta_file=concat_fasta, threads=args.t, clustal_format=True)
  
    logging.info("\n" + "8-) Now we will choose the most representative contig" + "\n")
    repr_contig_id, repr_contig_cluster = getReprContig.get_repr_contig("all_mitogenomes.rotated.fa", args.t)
    logging.info("Representative contig is {} that belongs to {}. This contig will be our final mitogenome. See all contigs and clusters in cdhit.out.clstr".format(repr_contig_id, repr_contig_cluster))
    
    repr_contig_fasta = repr_contig_id + ".mitogenome.rotated.fa"
    repr_contig_get_gb = ["mitofinder", "--new-genes", "--max-contig-size", str(max_contig_size), "-j", "final_mitogenome.annotation", "-a", repr_contig_fasta, "-r", args.g, "-o", args.o, "-p", str(args.p)]
    subprocess.run(repr_contig_get_gb, stderr=subprocess.STDOUT)

    final_fasta = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.fasta")
    final_gbk = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.gb")
    
    # Generate contigs stats
    ## Print first three lines (comment and header)
    frameshifts = findFrameShifts.find_frameshifts(final_gbk)  
    contig_len, num_genes = findFrameShifts.get_gb_stats(final_gbk)
    if not frameshifts:
        all_frameshifts = "No frameshift found"
    elif len(frameshifts)==1:
        all_frameshifts = "".join(frameshifts)
    elif len(frameshifts)>1:
        all_frameshifts = ";".join(frameshifts)
    with open("contigs_stats.tsv", "w") as f: 
        f.write("# Related mitogenome is {} bp long and has {} genes\n".format(rel_mito_len, rel_mito_num_genes))
        f.write("\t".join(["contig_id", "frameshifts_found", "genbank_file", "length(bp)", "number_of_genes\n"]))
        f.write("\t".join(["final_mitogenome", all_frameshifts, "final_mitogenome.gb", contig_len, num_genes+"\n"]))
    ## Iterate over each contig and print its info (ID, framshifts and genbank file used 
    ## to search for the frameshifts)
    
    # get list with all individual contig stats files 
    contigs_stats_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.individual.stats'):
            # skips addition of representative contig, which is the 
            # same as the final_mitogenome
            if curr_file.split('.')[0] != repr_contig_id: 
                contigs_stats_files.append(curr_file)
    logging.debug(f"contigs_stats_files: {contigs_stats_files}") # for debugging
    
    with open("contigs_stats.tsv", "a") as outfile:
        for contig_stats in contigs_stats_files:
            with open(contig_stats, "r") as infile:
                shutil.copyfileobj(infile, outfile)

    # copying final FASTA and GBK to working directory
    shutil.copy(final_fasta, "final_mitogenome.fasta")
    shutil.copy(final_gbk, "final_mitogenome.gb")

    # cleaning up working directory 
    cleanUpCWD.clean_up_work_dir(contigs_ids)

    logging.info("Pipeline finished!" )
    runtime= time.time() - start_time
    logging.info(f"Runtime: {runtime} seconds")

if __name__ == '__main__':
    main()
