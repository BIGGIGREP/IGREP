"""
This will be a pipeline for annotating RAW Fastq files from MISEQ 2x300 sequencer using GGLAB IgFFT program
Steps of pipeline:
	1) Pass in R1, R2 files
	2) Unzip R1, R2 files
	3) (Optionally)Run trimmomatic on sequences
	4) Stitch R1-R2 reads
	5) Quality filter sequences
	6) Run annotation
	7) Run pairing
"""

import sys
import os
# LOAD IGREP SCRIPTS INTO PYTHON PATH
script_path = os.path.dirname(os.path.abspath(__file__))
igrep_scripts = os.path.join(os.path.dirname(script_path), 'common_tools')
sys.path.insert(0, igrep_scripts)

# IMPORT SCRIPTS USED IN ANALYSIS
import immunogrep_ngs_pair_tools as processing
import immunogrep_fastx_toolkit as fastx
import immunogrep_igfft_functions as igfft
import immunogrep_useful_functions as useful

# Set up some default settings to run program
# How many threads should we use when multithreading is possible
number_threads = 3

# Trimmomatic settings: Trimmoamtic will remove low quality bases from ends of sequences using a sliding window approach
# Size of window in trimmomatic
trim_seqs = False
window_trim = 5
# average quality of window
quality_cutoff_trim = 20
# minimum lenght of sequence after trim
min_read_len_post_trim = 100
# encoding of quality scores from sequencer
phred_encode = 33

# PEARING settings: Pear will stitch together R1-R2 reads
min_overlap_length = 10
max_assembly_length = 500
min_assembly_length = 50
max_fraction_uncalled = 1
pear_memory = '1G'

# Quality filtering settings: fastx quality filter will remove whole sequences whose percent of bases are below the quality cutoff
quality_cutoff = 20
percent_bases = 50

# AB Annotation settings
# Use these sequences to determine isotype
isotyping_barcodes = [
	['CCTCCACCAAGGGCCCATCGCAG', 'IGHG'],
	['GGAGTGCATCCGCCCCAACC', 'IGHM'],
	['CATCCCCGACCAGCCCCAAGC', 'IGHA'],
	['GAACTGTGGCTGCACCATCT', 'IGK'],
	['GTCACTCTGTTCCCGCCCTC', 'IGL']
]
remove_insertions = 0  # Never remove insertions from alignment
#remove_insertions = 1  # Aloways remove insertions from alignment
#remove_insertions = 2  # Only remove insertions if there is a stop codon in final sequence

# Pairing settings
# How are CDRH3 sequences clustered (min cluster, max cluster, step size)
cluster_setting = [0.85, 0.98, 0.01]
# This is the cluster we usually want and will store in database if possible
annotation_cluster_cutoff = 0.9


def run_gglab_pipeline(input_files, species, loci, group_name=''):
	# Unzip files
	print('Processing raw fastq files')
	processed_files = []
	
	for pair_of_files in input_files:		
		folder_path = os.path.dirname(pair_of_files[0])
		for i, f in enumerate(pair_of_files):		
			if f.endswith('.gz'):
				print('Unzipping: ', f)
				pair_of_files[i] = useful.gunzip_python(f)

		# Run trimmomatic
		if trim_seqs:
			print('Trimming low quality bases')
			trimming_parameters = {
				'SLIDINGWINDOW': str(window_trim) + ':' + str(quality_cutoff_trim),				
				'MINLEN': min_read_len_post_trim
			}
			method = 'PE'		
			input_files = processing.run_trimmomatic(pair_of_files, folder_path, method, phred_encode, trimming_parameters)

		# Stitch R1-R2 files
		pairing_parameters = {
			'v': min_overlap_length,
			'm': max_assembly_length,
			'n': min_assembly_length,
			'u': max_fraction_uncalled,					
		}
		print('Stitching R1-R2 reads')
		pear_results = processing.run_pear(pair_of_files[0], pair_of_files[1], working_directory=folder_path, parameters=pairing_parameters, num_threads=number_threads, memory=pear_memory)[0]		
		# Run quality filtering
		filtered_file = fastx.Run_Quality_Filter(pear_results, output_dir=folder_path, quality=quality_cutoff, percent=percent_bases)		
		os.remove(pear_results)
		processed_files.append(filtered_file)
		processed_files.append('/home/costas/Documents/testpairing/HD150K.fasta')
	
	print('Annotating processed fastq files')
	annotated_files = []
	for i, f in enumerate(processed_files):
		annotated_f = igfft.igfft_multiprocess(f, species=species, locus=loci, parsing_settings={'isotype': isotyping_barcodes, 'remove_insertions': remove_insertions}, num_processes=number_threads, delete_alignment_file=True)			
		annotated_files.append(annotated_f[0])
	print('Pipeline complete')
		
if __name__ == "__main__":
	arguments = sys.argv[1:]
	argnum = 0
	fastq_files = []
	while True:
		if argnum >= len(arguments):
			break
		if arguments[argnum] == '-species':
			argnum += 1
			species = arguments[argnum]
		elif arguments[argnum] == '-locus': 
			argnum += 1
			locus = arguments[argnum]
		else:
			fastq_files.append(arguments[argnum].split(','))
			if len(fastq_files[-1]) > 2:
				raise Exception('When providing input files, on provide two R1-R2 files seperated by commas: i.e file1_r1.fastq,file1_r2.fastq file2_r1.fastq,file2_r2.fastq')
			for i, f in enumerate(fastq_files[-1]):
				f = os.path.expanduser(f)
				fastq_files[-1][i] = os.path.relpath(f)
				if not os.path.isfile(fastq_files[-1][i]):
					raise Exception('The provided file {0} does not exist'.format(fastq_files[-1][i]))
		argnum += 1
	
	run_gglab_pipeline(fastq_files, species, locus.split(','))
