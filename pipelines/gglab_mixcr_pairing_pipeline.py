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
import immunogrep_mixcr_tools as mixcr
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
# AB Annotation settings
exportPrettyAlignment = False

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

		# Run trimmomatic
		trimming_parameters = {
			'SLIDINGWINDOW': str(window_trim) + ':' + str(quality_cutoff_trim),				
			'MINLEN': min_read_len_post_trim
		}
		method = 'SE'		
		trimmedf = processing.run_trimmomatic(f, folder_path, method, phred_encode, trimming_parameters)[0]		
		# Run quality filtering
		filtered_trimmed_file = fastx.Run_Quality_Filter(trimmedf, output_dir=folder_path, quality=quality_cutoff, percent=percent_bases)		
		os.remove(trimmedf)
		processed_files.append(filtered_trimmed_file)
	
	print('Annotating processed fastq files')
	annotated_files = []
	for i, f in enumerate(processed_files):
		output_file = useful.removeFileExtension(f) + '.mixcr.alignment'
		output_file_annotation = useful.removeFileExtension(f) + '.mixcr.annotation'
		# Run MIXCR file
		print('Running MIXCR')
		[annotated_f, command_val] = mixcr.RunMixcr(f, output_file, filetype='FASTQ', loci=[], species='', exportPrettyAlignment=False, num_threads=number_threads)
		# Parse MIXCR file
		print('Parsing MIXCR')
		annotated_file = mixcr.parseMIXCR(f, output_file, 'FASTQ', output_file_annotation, command_val=command_val)  # again, annotated_file should be equal to outfile_annotation
		annotated_files.append(annotated_file)	
	print('Pairing sequences')	
	output_dir = os.path.dirname(annotated_files[0])
	pairing.RunPairing(annotated_files, annotated_file_formats='TAB', analysis_method='MIXCR', output_folder_path=output_dir, prefix_output_files=group_name, cluster_cutoff=cluster_setting, annotation_cluster_setting=annotation_cluster_cutoff)
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
			fastq_files.append(arguments[argnum])
			f = os.path.expanduser(fastq_files[-1])			
			if not os.path.isfile(f):
				raise Exception('The provided file {0} does not exist'.format(f))
			fastq_files[-1] = os.path.abspath(f)			
		argnum += 1	
	run_gglab_pipeline(fastq_files, species, locus.split(','))
