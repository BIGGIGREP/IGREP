"""
This will be a pipeline for annotating RAW Fastq files from MISEQ 2x300 sequencer using GGLAB IgFFT program and then subsequntly running GGLAB pairing.
Steps of pipeline:
	1) Pass in R1, R2 files
	2) Unzip R1, R2 files
	3) Run trimmomatic on sequences
	4) Quality filter sequences
	5) Run annotation
	6) Run pairing
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
import immunogrep_gglab_pairing as pairing
import immunogrep_useful_functions as useful

# Set up some default settings to run program
# How many threads should we use when multithreading is possible
number_threads = 6

# Trimmomatic settings: Trimmoamtic will remove low quality bases from ends of sequences using a sliding window approach
# Size of window in trimmomatic
window_trim = 5
# average quality of window
quality_cutoff_trim = 20
# minimum lenght of sequence after trim
min_read_len_post_trim = 100
# encoding of quality scores from sequencer
phred_encode = 33

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
	for i, f in enumerate(input_files):
		folder_path = os.path.dirname(f)
		if f.endswith('.gz'):
			print('Unzipping: ', f)
			f = useful.gunzip_python(f)
		annotated_f = igfft.igfft_multiprocess(f, file_type='FASTQ', species=species, locus=loci, parsing_settings={'isotype': isotyping_barcodes, 'remove_insertions': remove_insertions}, num_processes=number_threads, delete_alignment_file=True)			
		annotated_files.append(annotated_f[0])
	output_file_list = ','.join(annotated_files)
	print output_file_list
	return output_file_list

		
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
			fastq_files.append(arguments[argnum])
			f = os.path.expanduser(fastq_files[-1])			
			if not os.path.isfile(f):
				raise Exception('The provided file {0} does not exist'.format(f))
			fastq_files[-1] = os.path.abspath(f)			
		argnum += 1	
	run_gglab_pipeline(fastq_files, species, locus.split(','))
