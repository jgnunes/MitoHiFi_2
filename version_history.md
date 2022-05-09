mitohifi-v3-fromCirc_v02.py	runs until (included) circularization
mitohifi-v3-fromCirc_v02.1.py	runs until (included) annotation
mitohifi-v3-fromCirc_v02.2.py	runs until (included) rotation 
mitohifi-v3-fromCirc_v02.3.py 	runs until (included) alignment of rotated contigs 
mitohifi-v3-fromCirc_v02.4.py	offers reads or contigs mode for running the pipeline
mitohifi-v3-fromCirc_v02.5.py	don't raise exception when contig has not tRNA-PHe. Instead, just raise a warning and skip contig
mitohifi-v3-fromCirc_v02.6.py	removes the -l flag and calculates the length of the related mito
mitohifi-v3-fromCirc_v02.7.py	changes MAFFT output to clustalW format
mitohifi-v3-fromCirc_v02.8.py	changes maximum contig size of MitoFinder to 3 times the length of the related mito
mitohifi-v3-fromCirc_v02.9.py	includes -p PROCESSORSTOUSE option in MitoFinder annotation
mitohifi-v3-fromCirc_v02.10.2.py	includes CDHIT step to pick most representative contig	
mitohifi-v3-fromCirc_v02.11.py	cleans up working directory and copies final (FASTA and GBK) mitogenome to it
mitohifi-v3-fromCirc_v02.11.5.3.py	calculates frameshifts stats for processed contigs
mitohifi-v3-fromCirc_v02.11.5.4.py	deals with contigs that don't have tRNA-Phe
mitohifi-v3-fromCirc_v02.11.5.5.py	pass more informative message about tRNA/annotation missing contigs
mitohifi-v3-fromCirc_v02.11.5.6.py	remove shell=True from subprocesses 		
mitohifi-v3-fromCirc_v02.11.5.7.py	integrate 02.11.5.5 and 02.11.5.6
mitohifi-v3-fromCirc_v02.11.5.8.py	fixes bug on contig mode when input contigs have a header with _rc_ pattern	
mitohifi-v3-fromCirc_v02.11.6.py	includes HiFiasm bloom filter option (as recommended by Cristo Gallardo on Github) and replaces os.rename by shutil.move on line 177
mitohifi-v3-fromCirc_v02.11.7.py	includes --new-genes option in MitoFinder call to annotate non-standard animal mitochondrial genes (e.g. rps3 in fungi) (in progress)
mitohifi-v3-fromCirc_v02.11.8.py	includes new stats in contigs_stats.tsv file 
mitohifi-v3-fromCirc_v02.11.9.py	supports running the script using symlinks (besides rel/abs paths)
mitohifi-v3-fromCirc_v02.12.3.py	parallelizes annotation step
mitohifi-v3-fromCirc_v02.13.1.py	includes related mitogenome info on contigs_stats.tsv
mitohifi-v3-fromCirc_v02.13.2.py	fixes bugs
mitohifi-v3-fromCirc_v02.13.3.py	creates flags for plants/animal (-a animal/plant) 
mitohifi-v3-fromCirc_v02.13.4.py	deal with cases where none contig has tRNA-Phe gene
mitohifi-v3-fromCirc_v02.13.5.py	improves logging
mitohifi-v2.2	changed the version numbering to match the version number we're using withing tola group
mitohifi-v2.3|4	coverage plots for final_mitogenome; annotation plots for all contigs; update minimap call to use hifi flag 
