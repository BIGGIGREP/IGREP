import immunogrep_read_file as readfile
import pandas as pd 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy
from collections import defaultdict
from collections import OrderedDict
import immunogrep_useful_functions as useful
import os 
import subprocess
import datetime


import gc
import time

#for making folders in igrep
import immunogrep_file_system_tools as filesystem

try:
	import appsoma_api
	appsoma_api.resource_pull("https://www.appsoma.com/programs/get/cchrysostomou/default/make_unique_aa.bash",'process_unique_ab_aa_file.bash')
except:
	pass

delim = '_:_'

#this is a simple variable for mapping the fields we require to run the analysis to the field names in different annotation files 
#so when reading files we know that the field name 'Functionality_1' in IMGT files correponds to our functionality variable we use in the analysis 
fields_for_analysis = {
	'IMGT':{		
		'functionality':'Functionality_1',
		'vgene':'V-GENE and allele_1',
		'jgene':'J-GENE and allele_1',
		'dgene':'D-GENE and allele_1',		
		'cdr3':'CDR3-IMGT_5',			
		'cdr1':'CDR1-IMGT_5',			
		'cdr2':'CDR2-IMGT_5',			
		'shm':'V-REGION Nb of mutations_8',					
		'recomb':'R_TYPE',
		'locus':'',
		'full_len_ab':'PREDICTED_AB_SEQ.AA'#THIS FIELD IS NOT PRESENT IN THE IMGT FILE, BUT IT DOES GET ADDED TO THE INFORMATION WHEN READING IMGT FILES USING OUR READER CLASS					
	},
	'IGBLAST':{		
		'functionality':'PRODUCTIVE',
		'vgene':'VREGION.VGENES',
		'jgene':'JREGION.JGENES',
		'dgene':'DREGION.DGENES',		
		'cdr3':'CDR3.AA',		
		'cdr2':'CDR2.AA',		
		'cdr1':'CDR1.AA',		
		'shm':'VREGION.SHM_NT',			
		'recomb':'RECOMBINATION_TYPE',
		'locus':"LOCUS_NAME",
		'full_len_ab':'PREDICTED_AB_SEQ.AA',							
	},
	'IGFFT':{
		'functionality':'Productive',
		'vgene':'Top_V-Gene_Hits',
		'dgene':'',
		'jgene':'Top_J-Gene_Hits',				
		'cdr3':'CDR3_Sequence.AA',
		'cdr2':'CDR2_Sequence.AA',
		'cdr1':'CDR1_Sequence.AA',		
		'shm':'VRegion.SHM.NT',
		'full_len_ab':'Full_Length_Sequence.AA',
		'recomb':'Recombination_Type',
		'locus':'Locus'
	},
	'MIXCR':{		
		'functionality':'Productivity',
		'vgene':'FirstVgene',
		'jgene':'FirstJgene',
		'dgene':'FirstDgene',		
		'cdr3':'AA. seq. CDR3',
		'cdr2':'AA. seq. CDR2',
		'cdr1':'AA. seq. CDR1',
		'seq_header':'Seqheader',
		'shm':'VGENE: Shm.nt',
		'full_len_ab':'Full AA',	
		'locus':'',
		'recomb':'Recombination Type'
	},
	'DATABASE':{		
		'functionality':'PRODUCTIVE',
		'vgene':'VREGION.VGENES',
		'jgene':'JREGION.JGENES',
		'dgene':'DREGION.DGENES',		
		'cdr3':'CDR3.AA',
		'cdr2':'CDR2.AA',
		'cdr1':'CDR1.AA',
		'seq_header':'SEQUENCE_HEADER',
		'shm':'VREGION.SHM_NT',
		'full_len_ab':'PREDICTED_AB_SEQ.AA',
		'locus':'LOCUS',
		'recomb':'RECOMBINATION_TYPE'
	}
	
}

#variable for determining whether the antibody is heavy or light chain 
recomb_call = {
	'TRB':'VDJ',
	'TB':'VDJ',
	'IGH':'VDJ',
	'TRD':'VDJ',
	'IGK':'VJ',
	'TRA':'VJ',
	'TA':'VJ',
	'IGL':'VJ',
	'TRG':'VJ',
}

def ProcessGene(gene):
	"""
		Removes the allele calls from genes
		
		Assumption: 
			
			Multiple genes are seperated by ',' 
			
			Alleles are seperated by '*' 
			
			If a gene is seperated by multiple spaces, then the gene should be identified by a gene that contains either - or '*' 
			For example: 

				Imgt genes may be: 

					Homo sapiens IGHV1-3*01

					We only want to isolate the word IGHV1-3
	"""
	
	if not gene:
		return ''	
	#extract all genes in field 
	gene_array = gene.split(',')
	for gene_num,each_gene in enumerate(gene_array):
		#split each gene by spaces. Go through each word
		split_words  = each_gene.split(' ')
		if len(split_words)==1: #there is only one word 
			gene_array[gene_num] = split_words[0].split('*')[0]
		else:
			g = ''
			for subv in split_words:
				#if we find a word with gene characters in it, it must be our gene 
				if '*' in subv or '-' in subv:						
					#extract everything before '*'
					g = subv.split('*')[0]
					break
			gene_array[gene_num] = g	
	return ','.join(gene_array)
	

def defaultdictcdr3(num_exp):
	return [0]*(num_exp)
			
def defaultdictgenes(num_exp):
	return [0]*(num_exp)

def FloatRange(min,max,step,num_dec=None):
	"""
		Uses a generator to return a series of floating numbers starting from min, increasing in increments of step, and continuing until max 
	"""
	min = float(min)
	max = float(max)
	x = min
	if max<min:
		min2=min
		min=max
		max=min2
	if min==max:
		yield min-1
		yield min+1
	else:
		while x<=max:
			if num_dec:
				yield round(x,num_dec)
			else:
				yield x
			x+=float(step)
			

def PlotGeneDist(dataframe,figure_path,plot_title='',xlabel='',ylabel='',max_val=None,min_val=None,step=None):	
	#replace nan and infinity values with 0
	g = dataframe.replace([np.inf, -np.inf,np.nan], 0)
		
	if len(g)==0:
		return
	if not max_val:
		max_val = max(g.max())
	if not min_val:
		min_val = min(g.min())
	if not step:
		step = 10 #10 (max_val-min_val)/float(len(g)) #10
	
	# Create a figure of given size
	fig = plt.figure(figsize=(16,24))
	# Add a subplot
	ax = fig.add_subplot(111)
	# Set title		
	
	g.plot(kind='barh',ax=ax, title=plot_title,width=0.6)
	# Remove grid lines (dotted lines inside plot)
	ax.grid(False)
	
	# Remove plot frame
	ax.set_frame_on(False)
	# Pandas trick: remove weird dotted line on axis
	if len(ax.lines)>0:
		ax.lines[0].set_visible(False)
	
	# Customize title, set position, allow space on top of plot for title
	ax.set_title(ax.get_title(), fontsize=20, alpha=0.7, ha='left')
	plt.subplots_adjust(top=0.9)
	ax.title.set_position((0,1))
	
	# Set x axis label on bottom of plot, set label text
	#ax.xaxis.set_label_position('bottom')	
	if xlabel:
		ax.set_xlabel(xlabel, fontsize=14, alpha=0.7, ha='left')
	if ylabel:
		ax.set_ylabel(ylabel, fontsize=14,alpha=0.7)
	
	#ax.xaxis.set_label_coords(0, 1.04)

	# Position x tick labels on bottom
	ax.xaxis.tick_bottom()
	# Remove tick lines in x and y axes
	ax.yaxis.set_ticks_position('none')
	ax.xaxis.set_ticks_position('none')
	
	# Customize x tick lables
	
	xticks = [r for r in FloatRange(min_val,max_val,(float(max_val)-float(min_val))/step,4)]
	
	ax.xaxis.set_ticks(xticks)
	#xticks = [item.get_text() for item in ax.get_xticklabels()]
	
	ax.set_xticklabels(xticks, fontsize=10, alpha=1)
	
	# Customize y tick labels
	yticks = [item.get_text() for item in ax.get_yticklabels()]
	ax.set_yticklabels(yticks, fontsize=10, alpha=1)
	ax.yaxis.set_tick_params(pad=12)
	
	# Shrink current axis by 20%
	box = ax.get_position()
	#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	
	# Put a legend to the right of the current axis
	#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
	plt.savefig(figure_path+'.png',bbox_inches='tight',dpi=100)
	plt.savefig(figure_path+'.svg',format='svg',bbox_inches='tight',dpi=100)
	

def PlotVJGeneHeatMap(dataframe,figure_path,max_val=None,min_val=None):

	subtable = dataframe	
	if len(subtable)==0:
		return
	#list of columns in current table. each column should correspond to an experiment
	c = subtable.columns 

	

	#set the figure size 
	fig = plt.figure(figsize=(30,30))				

	
	#cbar = fig.colorbar(ticks=[-1, 0, 1], orientation='horizontal')    
	if not max_val:
		max_val = max(subtable.max())
	if not min_val:
		min_val = min(subtable.min())
	
	#loop through all of the experiments. We want to make a VJ gene heatmap for each one 
	for i,name in enumerate(c):			
		#use PIVOT to change the dataframe table. Basically make a new table whose rows are JGENSE and whose columns are VGENES. the value of each talbe is the frequency for that experiment 
		exp_table = subtable.reset_index().pivot(index='JGENE',columns='VGENE',values=name)			
		exp_table = exp_table.fillna(0)
								
		#subplot => (numrows,numcolumns,current plot)
		ax = fig.add_subplot(len(subtable.columns),1,i+1)
		
		#norm => matplotlib.colors.LogNorm() => maek a log scale for color bar 
		#create a heatmap use the table 
		cax = ax.imshow(exp_table,cmap=matplotlib.cm.Blues,interpolation='nearest')
		#set the colorlimit on the figure
		cax.set_clim(vmin=min_val,vmax=max_val)
		
		row_labels = list(exp_table.index)    
		column_labels = list(exp_table.columns)
		n = len(column_labels)
		m = len(row_labels)
		
		#put the colored boxes/heatmap in the middle of the ticks
		ax.set_xticks(np.linspace(0,n-1,n))
		ax.set_yticks(np.linspace(0,m-1,m))
	
		#set up the xtic labels to match the labels
		labels = ax.set_xticklabels(column_labels, minor=False,fontsize=10, alpha=1)
		ax.set_yticklabels(row_labels, minor=False,fontsize=10, alpha=1)
	
		#make the xaxis rotated
		for label in labels:
			label.set_rotation(70) 			
		
		ax.title.set_position((0,1))
		ax.set_title(name, fontsize=14, alpha=0.7, ha='left')
					
		ax.set_ylabel('JGENES', fontsize=14, alpha=0.7,)
	
	plt.tight_layout() #=> took out tight_layout because of warning message, used subplot adjust instead
	#add vgene label to the last plot		
	ax.set_xlabel('VGENES', fontsize=14, alpha=0.7,)
	#create an axis to draw the colormap	
	colormap_box = fig.add_axes([1.05,0.5 , 0.02, 0.35])
	cb = fig.colorbar(cax,cax=colormap_box)
	cb.set_label("Frequency", size=14)
	
	cb.ax.tick_params(axis='x', labelsize=10)
	
	#left  = 0.125  # the left side of the subplots of the figure
	#right = 0.9    # the right side of the subplots of the figure
	#bottom = 0.1   # the bottom of the subplots of the figure
	#top = 0.9      # the top of the subplots of the figure
	#wspace = 0.2   # the amount of width reserved for blank space between subplots
	#hspace = 0.5   # the amount of height reserved for white space between subplots
	#plt.subplots_adjust(hspace=10)
	plt.savefig(figure_path+'.png',bbox_inches='tight',dpi=100)
	plt.savefig(figure_path+'.svg',format='svg',bbox_inches='tight',dpi=100)
	plt.clf()


#calculate cdr3 length information 
def PlotCDR3Histogram(cdr3,figure_path):					
	#step1=> filter out all CDR3lengths<=2
	cdr3=cdr3[cdr3['CDR3_LENGTH']>2]
			
	#setup figure
	fig = plt.figure(figsize=(20,10))
		
	#step3=>Make histograme plot for each VDJ and VJ 
	plotnum = 1
	xlabels=['CDRH3 Length','CDRL3 Length']
	cdr3_length_stats = {'VDJ':{'mean':{},'std':{}},'VJ':{'mean':{},'std':{}}}
	
	#step2=> split dataframe into VDJ adn VJ 	
	rtype = ['VDJ','VJ']
	for l,r in enumerate(rtype):                    				
		if r not in list(cdr3.index.levels[cdr3.index.names.index('recomb')]):
			continue
		df = cdr3.loc[r]
			
		
		#these are a column of cdr3 lengths observed in ALL experiments
		total_count = df['CDR3_LENGTH']
			
		if len(total_count)==0:
			#no experiments have any sequences with this recombination type 
			continue
		
		cdr3_length_stats[rtype[l]]['mean']['TOTAL_COUNTS'] = total_count.mean()
		cdr3_length_stats[rtype[l]]['std']['TOTAL_COUNTS'] = total_count.std()
		for cnames in df.columns:
			#if cnames=='TOTAL_COUNTS':
			#	continue
			cdr3_length_stats[rtype[l]]['mean'][cnames] = df.loc[df[cnames]>0,'CDR3_LENGTH'].mean()
			cdr3_length_stats[rtype[l]]['std'][cnames] = df.loc[df[cnames]>0,'CDR3_LENGTH'].std()
				
		#remove total counts and cdr3 column
		#IGNORE the number of counts we observe each CDR3  , we just want to analyze unique cdr3s
		sub_df = df.reset_index().drop(['TOTAL_COUNTS','CDR3'],1).set_index('CDR3_LENGTH').applymap(lambda x:1 if x>0 else 0)        
		#group by cdr3 length and sum all columns in group (this will work b/c of the unique 1/0 in apply map)
		sub_df = sub_df.groupby(level=0).sum()    
		#setup axes
		ax = fig.add_subplot(1,2,plotnum)                                
		max_cdr3_len = df['CDR3_LENGTH'].max()        
		index_vals = sub_df.index        
		#add 0 to missing cdr3 lengths...its a hack to force bar plot to shift this many x -coordinates
		for j in [r for r in range(min(index_vals)) ]: 
			sub_df.loc[j] = 0                         
		sub_df.sort().plot(ax=ax, kind='bar',alpha=0.7)                					
		
		max_val = max_cdr3_len+3
		xticks = [r for r in FloatRange(0,max_val,5)]

		ax.set_xlim([0, max_val])             
		ax.set_xlabel(xlabels[l], fontsize=14)
		
		[n,bins] = np.histogram(total_count)            
		max_hist_height = max(n)    
		max_all_exp_counts = max(sub_df.max())    
		normalize_histogram_height = max_all_exp_counts/float(max_hist_height)
			
		# Plot the resulting histogram
		#center = (bins[:-1]+bins[1:])/2
		
		width = 1*(bins[1]-bins[0])
		heights = [bn*normalize_histogram_height for bn in n]
		plt.bar(bins[:-1],heights,width,color='black',alpha=0.15)
		
		ax.xaxis.set_ticks(xticks)
		ax.xaxis.set_ticklabels([str(int(l)) for l in xticks])
			
		for tick in ax.xaxis.get_major_ticks():
			tick.label.set_fontsize(10)
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(10)    
		if plotnum==1:
			#only set the y axis label for the first plot
			ax.set_ylabel('Counts', fontsize=14)
		
		# Put a legend to the right of the current axis
		ax.legend(loc='center left', bbox_to_anchor=(0, 1.1))
		plt.grid('off')
		ax.set_xlim([0, df['CDR3_LENGTH'].max()+3])
		plotnum+=1
	
		del df
	
	plt.savefig(figure_path+'.png',bbox_inches='tight',dpi=100)
	plt.savefig(figure_path+'.svg',format='svg',bbox_inches='tight',dpi=100)
	return cdr3_length_stats
	
def CalculateDiversities(cdr3_df,figure_path):
	#make a total of 6 plots
	#column 1= > VDJ results 
	#column 2= > VJ results 
	#row 1 = > DIVERSITY calls 
	#row 2 = > plot of unique cdr3 and cum sum 
	#row 3 => plt of # cdr3s > x count 
	
	diversity_results = {}
	r_types = ['VDJ','VJ']
	fig = plt.figure(figsize=(40,40))
	#NOW  plot everythiong
	plot_nums = [[1,3,5],[2,4,6]]
	
	#loop through recombination types 
	for num,r in enumerate(r_types):
		if r not in list(cdr3_df.index.levels[cdr3_df.index.names.index('recomb')]):
			#no results for that recombination type
			continue
			
		#take slice of cdr3 (isolate VDJ and VJ plots separately)
		temp = cdr3_df.loc[r]
		
		try:
			temp = temp.drop('CDR3_LENGTH',axis=1)
		except:
			pass
		
		#ANALYSIS 1
		#make a dictionary whose keys refer to counts and whose values refers to => nubmer of unique sequences in each experiment with counts above this 
		counts_above = [1,2,8,16,50,100,500,1000]
		total_seq_count = (temp>0).sum()	
		
		#loop through the total_seq_count:
		#if there is essentially only one experiment, (all experiments are 0 counts or only one epxerimenrt was analyzed) then exclude TOTAL_COUNTS
		num_above_0 = (total_seq_count>0).sum()
		
		if num_above_0 <= 2 and 'TOTAL_COUNTS' in temp.columns:
			temp = temp.drop('TOTAL_COUNTS',axis=1)
		
			
				
		counter_dict = {np.log2(c):((temp>=c).sum())/total_seq_count for c in counts_above}
		#counter_dict[0] = ((temp>0).sum())/((temp>0).sum())
					
		#convert table into frequencies
		temp = temp/(temp.fillna(0).sum()) 
		
		#ANALYSIS 2
		#calculate various diversity measurements for cdr3s
		richness = temp.applymap(lambda x:1 if x>0 and np.abs(x)!=np.inf else 0).sum()
		shannon_index = temp.applymap(lambda x:-1*np.log(x)*x).sum()
		ginni_simpsons_index = 1-temp.applymap(lambda x:x*x).sum()
		#Normalize diversity indices using the "true diversity" calls 
		shannon_diversity = (np.exp(shannon_index)).fillna(0).replace([np.inf, -np.inf],0)
		simpsons_diversity = (1/(1-ginni_simpsons_index)).fillna(0).replace([np.inf, -np.inf],0)
		richness_diversity = (richness).fillna(0).replace([np.inf, -np.inf],0)
		#remove NA numbers from diversity indices
		df_diversities = OrderedDict()
		df_diversities['Species Richness']=richness_diversity    
		df_diversities['Shannon entropy']=shannon_diversity
		df_diversities['Ginni Simpsons Index']=simpsons_diversity				
		
		#ANALYSIS 3 
		#sort each column by frequency and then calculate/plot the cumulative sum
		sub_dfs = OrderedDict()
		for exp in temp.columns:        
			sub_dfs[exp] = temp.loc[temp[exp]>0,exp].reset_index().drop('CDR3',axis=1).sort(exp,ascending=0,inplace=False).cumsum()							
		
		#1) MAKE A PLOT FOR DIVERSITY 
		ax = fig.add_subplot(3,2,plot_nums[num][0])
		#provide a title ot the first column of subplot
		if num == 0:
			ax.set_title('CDRH3 Results', fontsize=20, alpha=0.7)
		else:
			ax.set_title('CDRL3 Results', fontsize=20, alpha=0.7)				
		pd.DataFrame(df_diversities).transpose().plot(ax=ax,kind='bar',sharex=False)		
		ax.set_ylabel('Normalized Diversity', fontsize=14, alpha=1,)
		ax.set_xlabel('Diversity Index',fontsize=14,alpha=1)	
		for tick in ax.yaxis.get_major_ticks():        
			tick.label.set_fontsize(10)				
		labels = ax.set_xticklabels([k for k in df_diversities.keys()], fontsize=10, alpha=1)	
		for label in labels:
			label.set_rotation(0)        				
		plt.grid('off')	
		plt.tight_layout()	
		
		#2) MAKE A PLOT FOR FREQUENCY OF EACH UNIQUE CDR3 
		ax2 = fig.add_subplot(3,2,plot_nums[num][1])
		#each key in dictionary is an experiment
		#each value is a dataframe defining that experiments' CDF for unique CDR3s
		#for exp_name,df in sub_dfs.iteritems():
		#	df.plot(ax=ax2,sharex=False,sharey=False)    		
		for exp_name,ddf in sub_dfs.iteritems():
			if exp_name == 'TOTAL_COUNTS':			
				ax2.plot(ddf[exp_name],label = exp_name,linewidth=4,alpha=0.5)
			else:
				ax2.plot(ddf[exp_name],label = exp_name,linewidth=4)								
		#fix x and y axes
		for tick in ax2.xaxis.get_major_ticks():        
			tick.label.set_fontsize(10)
		for tick in ax2.yaxis.get_major_ticks():        
			tick.label.set_fontsize(10)
		ax2.legend(loc='lower right')
		ax2.set_ylabel('Cumulative Frequency of total CDR3\ncount in sample\n\n\n\n',fontsize=14, alpha=1,horizontalalignment='center')
		ax2.set_xlabel('Relative rank of unique CDR3 in each sample',fontsize=14, alpha=1)
		plt.grid('off')		
		plt.tight_layout()	
		
		#3) MAKE  third plot that shows a CDF of how many unique CDR3s have a count above x    
		ax3 = fig.add_subplot(3,2,plot_nums[num][2])    
		pd.DataFrame(counter_dict).transpose().plot(marker='o',ax=ax3,sharex=False,linewidth=4,markersize=12)    		
		max_tick = int(np.log2(max(counts_above)))+1    
		ticks = [tt for tt in range(1,max_tick)]    
		tick_labels = [np.power(2,tt) for tt in ticks]                    
		ax3.set_xticks(ticks)           
		ax3.set_xticklabels(tick_labels, minor=False,fontsize=10, alpha=1)        
		for tick in ax3.yaxis.get_major_ticks():        
			tick.label.set_fontsize(10)
		ax3.legend(loc='upper right')
		ax3.set_ylabel('CDF\n(Frequency of unique CDR3 sequences\n whose observed counts are >= x)\n\n\n\n',fontsize=14, alpha=1,horizontalalignment='center')
		ax3.tick_params(axis='y', pad=15)
		ax3.set_xlabel('CDR3 Count',fontsize=14, alpha=1)	
		plt.grid('off')
		plt.tight_layout()	
	
		diversity_results[r]  = {
		'unique_cdr3s':total_seq_count.to_dict(),
		'shannon_entropy':{
			'index':shannon_index.replace([np.inf, -np.inf],0).to_dict(),
			'true_diversity':shannon_diversity.to_dict()
		},
		'ginni_simpsons':{
			'index':ginni_simpsons_index.replace([np.inf, -np.inf],0).to_dict(),
			'true_diversity':simpsons_diversity.to_dict()
		},
		'num_above_2':(counter_dict[1]*total_seq_count).to_dict()
		
		}					
		
		del shannon_diversity
		del shannon_index
		del simpsons_diversity
	
	plt.tight_layout()	
	plt.savefig(figure_path+'.png',bbox_inches='tight',dpi=100)
	plt.savefig(figure_path+'.svg',format='svg',bbox_inches='tight',dpi=100)
							
	return diversity_results			
			
	
def Descriptive_Statistics(list_of_files,input_file_type,analysis_name='',exp_names = [],output_file_prefix='',fields={},statistics_to_run=['ab_aa','cdr3','vgene','jgene','vjgene']):
	analysis_name = analysis_name.upper()
	if input_file_type=='IMGT' and not isinstance(list_of_files[0],list):
		list_of_files = [list_of_files]
	elif not isinstance(list_of_files,list):
		list_of_files = [list_of_files]
		
	if len(exp_names)!=len(list_of_files):
		exp_names = []
	
	#by default, save results to the same folder as the input file
	if not output_file_prefix:		
		output_file_prefix = useful.removeFileExtension(list_of_files[0])
	
	analysis_name = analysis_name.upper()
	supported_analyses = fields_for_analysis.keys()
	if (not analysis_name or analysis_name=='CUSTOM' or analysis_name not in supported_analyses) and not fields:
		raise Exception('The required fields for the provided analysis, {0}, is not currently automated. Please explicity provide the fields names'.format(str(analysis_name)))
	
	#first we use default fields defined ehere
	if analysis_name in supported_analyses:
		fields_to_use = copy.deepcopy(fields_for_analysis[analysis_name])
	else:
		fields_to_use = {}
	#next we add in user defined fields just in case there are any changes/mistakes
	for f,name in fields.iteritems():
		fields_to_use[f] = name 
	
	
	
	filenames_to_use = [f[0] if isinstance(f,list) else f for f in list_of_files]
	print('Performing descriptive statistics at {0}.'.format(str(datetime.datetime.now())))
	print('Analyzing the following files:\n\t {0}'.format('\n\t'.join(filenames_to_use)))
	unique_aa_file = None 
	unique_cdr3_file = None 	
	v_gene_analysis = None
	j_gene_analysis = None
	vj_gene_analysis = None
	gene_analysis_plot = output_file_prefix
	plots_created = []
	gene_summary_file = output_file_prefix+'.summary_of_stats.txt'
	
	
	output_file_names = {}
	
	aa_files = ['AB AA SEQUENCE','RECOMBINATION_TYPE','LOCUS','CDR1','CDR2','CDR3','STOP CODONS','PRODUCTIVE','VGENES','DGENES','JGENES','TOTAL COUNTS']
	fields_order = ['full_len_ab','recomb','locus','cdr1','cdr2','cdr3','stopc','functionality','vgene','dgene','jgene']
	num_exp= len(list_of_files)
	if not exp_names:
		if input_file_type=='IMGT':
			pass
		else:			
			exp_names = []
			for file in list_of_files:
				count = 1
				str_file = os.path.basename(file)
				while True:					
					if str_file in exp_names:
						str_file = os.path.basename(file)+'_'+str(count)
						count+=1
					else:
						exp_names.append(str_file)
						break			
		
	if 'ab_aa' in statistics_to_run:
		intermediate_file =  output_file_prefix+'.unique_aa_file_temp'
		#first we will use a temp file/intermeidate file 
		output_file_names['ab_aa'] = open(intermediate_file,'w')
		#output_file_names['ab_aa'].write('\t'.join(aa_files)+'\n')
	
	
	cdr3analysis = True if 'cdr3' in statistics_to_run else False
	aaanalysis = True if 'ab_aa' in statistics_to_run else False
	
	vjgene_dict=defaultdict(lambda:defaultdictgenes(num_exp))
	
	#cdr3_dict=defaultdict(lambda:defaultdictcdr3(num_exp))	
	cdr3_dict_vdj = defaultdict(lambda:defaultdictcdr3(num_exp))	
	cdr3_dict_vj = defaultdict(lambda:defaultdictcdr3(num_exp))	
	cdr3_dict_unk = defaultdict(lambda:defaultdictcdr3(num_exp))	
	
	use_these_fields = fields_to_use.values()
	fields_to_use['stopc'] = 'stopc'
	num_results = [0]*(num_exp)
	num_cdr3 = [0]*(num_exp)
	num_stop_codon = [0]*(num_exp)
	num_vdj = [0]*(num_exp)
	num_vj = [0]*(num_exp)
	num_sequences = [0]*(num_exp)
	
	
	if not fields_to_use['recomb']:
		#maybe the user never defined a feild for recombinoation type..that coudl be a problem because we will have to guess it using the variable at the top of the script: recomb_call
		recomb_not_defined = True	
		fields_to_use['recomb'] = 'recomb'
	else:
		recomb_not_defined = False
	
	
	print('Reading through sequences in file(s)')
	seqnum=1
	#go through all of the files and report the relevant fields 
	#if we are creating a unique amino acid file, then report thiese fields to temp file
	for fnum,each_file in enumerate(list_of_files):				
		annotated_file = readfile.immunogrepFile(each_file,input_file_type,field_names = use_these_fields)
		#loop through each file 
		for seq_lines in annotated_file.read():						
			if not seq_lines:
				continue
			if seqnum%500000==0:
				print('Read {0} sequences'.format(str(seqnum)))
			seqnum+=1
			num_sequences[fnum]+=1			
			seq_lines = defaultdict(str,seq_lines)
			if seq_lines[fields_to_use['full_len_ab']]:
				#full length antibody sequence not found
				num_results[fnum]+=1				
												
			#only select the first gene in the list. alos remove allelic name ('*')
			seq_lines[fields_to_use['vgene']] = seq_lines[fields_to_use['vgene']].split(',')[0].split('*')[0]
			seq_lines[fields_to_use['dgene']] = seq_lines[fields_to_use['dgene']].split(',')[0].split('*')[0]
						
			#IF NO RECOMBINATION TYPE IS FOUND or provided, THEN guess it using the vgene or jgene call
			if recomb_not_defined or not seq_lines[fields_to_use['recomb']]:
				r = '' #not sure what the recombation type is yet
				#try to guess the recombination type 				
				if seq_lines[fields_to_use['vgene']]:
					#use vgene if present
					# look at the first three characters in vgene to predict recombioation type
					gn = ProcessGene(seq_lines[fields_to_use['vgene']])
					if gn[:3] in recomb_call:
						r = recomb_call[gn[:3]]
					elif gn[:2] in recomb_call: #next check the first two letters (IGBLAST REPORTS TA RATHER THAN TRA for example 
						r = recomb_call[gn[:2]]
				if not r and seq_lines[fields_to_use['jgene']]:
					#still not r found, so use jgene 
					gn = ProcessGene(seq_lines[fields_to_use['jgene']])
					if gn[:3] in recomb_call:
						r = recomb_call[gn[:3]]
					elif gn[:2] in recomb_call: #next check the first two letters (IGBLAST REPORTS TA RATHER THAN TRA for example 
						r = recomb_call[gn[:2]]					
				
				#update recomb result 
				seq_lines[fields_to_use['recomb']] = r								
			
			if not seq_lines[fields_to_use['recomb']]:
				continue
								
			if seq_lines[fields_to_use['recomb']] == 'VDJ':
				num_vdj[fnum]+=1								
			elif seq_lines[fields_to_use['recomb']] == 'VJ':				
				num_vj[fnum]+=1
				
			
			seq_lines[fields_to_use['jgene']] = seq_lines[fields_to_use['jgene']].split(',')[0].split('*')[0]
			seq_lines['stopc'] = 'YES' if '*' in seq_lines[fields_to_use['full_len_ab']] else 'NO'			
			if seq_lines['stopc'] == 'YES':
				num_stop_codon[fnum]+=1
			if aaanalysis:
				exp_str = str(fnum+1)
				#make an intermediate file where we only put the fields we want in the proper order from any file 
				#we will use this field for sorting afterwards
				#also output exp_num to account for which sequence came from which experiment 
				output_file_names['ab_aa'].write('\t'.join([seq_lines[fields_to_use[f]] for f in fields_order])+'\t'+str(exp_str)+'\n')
			if seq_lines[fields_to_use['vgene']] or seq_lines[fields_to_use['jgene']]:							
				key_v =delim.join([seq_lines[fields_to_use['vgene']],seq_lines[fields_to_use['jgene']],seq_lines[fields_to_use['recomb']]])
				vjgene_dict[key_v][fnum]+=1
			
			if not seq_lines[fields_to_use['cdr3']]:
				#no cdr3 found 	
				continue
			
			#add unique cdr3_recomb and vjgene info to dictionaires
			num_cdr3[fnum]+=1
			
			if cdr3analysis:
				key = seq_lines[fields_to_use['cdr3']]
				#key_cdr3 = delim.join([],seq_lines[fields_to_use['recomb']]])
				if seq_lines[fields_to_use['recomb']]=='VDJ':
					cdr3_dict_vdj[key][fnum]+=1
				elif seq_lines[fields_to_use['recomb']]=='VJ':
					cdr3_dict_vj[key][fnum]+=1					
				else:
					print('unknown recombination types: ',seq_lines[fields_to_use['recomb']])
					cdr3_dict_unk[key][fnum]+=1 
									  
			if seqnum>10000:
				break
				
					
	if aaanalysis:
		output_file_names['ab_aa'].close()
		print('Generating a file of unique AB amino acid sequences')
		unique_aa_file = output_file_prefix+'.unique_aa_file.txt'
		#Use some bash to make a unique amino acid file using sorting and then some awk 
		GenerateAAFile(intermediate_file,unique_aa_file,aa_files,exp_names)
		#number of amino acid sequences observed
		if not os.path.isfile(unique_aa_file):
			num_unique_aa = 0 
		else:
			num_unique_aa = useful.file_line_count(unique_aa_file)-1 #-1 => remove header row count
	
	#Now have some fun with pandas 	
	if set(['vgene','jgene','vjgene']) & set(statistics_to_run):
		#vjgene_dict format = {
			#'key' = 'vgene',_,'jgene',_,'recombtype'
			#value = [count,count] => a list of counts for presence of that key in EACH providced file/experiment. Length of list = number of experiments
		#}
		gene_df = pd.DataFrame(vjgene_dict).transpose()
		if 'VGENE' not in gene_df.columns:
			gene_df['VGENE'] = ''
		if 'JGENE' not in gene_df.columns:
			gene_df['JGENE'] = ''
		if 'recomb' not in gene_df.columns:
			gene_df['recomb'] = ''
		gene_df['TOTAL_COUNTS'] = gene_df.sum(axis=1)		
		gene_df = gene_df.reset_index()				
		gene_df = gene_df.apply(ModifyPDTable,axis=1,args=(['VGENE','JGENE','recomb'],delim))
		
		
		new_names = {}
		for f,v in enumerate(exp_names):
			new_names[f]=v
			#key = experiment index number
			#value = new name

		#rename the columns 0,1,...num experiments to match the experiment names 
		gene_df = gene_df.rename(columns=new_names)
		
		#format of gene_df:
			#index => no index set, just use default numbers
			#columns => start with column for each experiment, then add the following columns: VGENE, JGENE, recomb, TOTAL_COUNTS

		if 'vgene' in statistics_to_run:
			print('Performing V gene analysis')
			
			v_gene_analysis = output_file_prefix+'.vgenes.txt'
			#group elements by VH GENE CALLS and VL gene calls 
			sorted_v_counts =  gene_df.groupby(['recomb','VGENE']).sum()#.count()#.sort('VGENE',ascending=1)						
			
			#find out which level in multilevel index corresponds to 'VGENE' => looking at above code , it should be level 1 (recomb should be level 0)
			vgene_level = sorted_v_counts.index.names.index('VGENE')			
			
			#remove results where vGENE is empty
			if '' in list(sorted_v_counts.index.levels[vgene_level]):
				sorted_v_counts = sorted_v_counts.drop('',level='VGENE')			
			
			ignore_counts = ['TOTAL_COUNTS','JGENE']
			keep_col = [n for n in sorted_v_counts.columns if n not in ignore_counts]
			g = sorted_v_counts[keep_col]			
			
			#NOW PLOT the FREQUENCY for every exeprement 
			if 'VDJ' in list(g.index.levels[g.index.names.index('recomb')]):
				
				vdj_g = g.xs('VDJ',level='recomb')
				
				PlotGeneDist(vdj_g/vdj_g.sum(),gene_analysis_plot+'.vdj.vgenes','VH Gene Distribution','Frequency','V Gene',max_val=None,min_val=0)
				
				plots_created.append(gene_analysis_plot+'.vdj.vgenes.png') #.png extension is added in the function plotgenedist
				
			if 'VJ' in list(g.index.levels[g.index.names.index('recomb')]):
				
				vj_g = g.xs('VJ',level='recomb')
				PlotGeneDist(vj_g/vj_g.sum(),gene_analysis_plot+'.vj.vgenes','VL Gene Distribution','Frequency','V Gene',max_val=None,min_val=0)			
				plots_created.append(gene_analysis_plot+'.vj.vgenes.png') #.png extension is added in the function plotgenedist
			sorted_v_counts.reset_index().sort(['recomb','TOTAL_COUNTS'],ascending=[1,0]).iloc[:,:-1].to_csv(v_gene_analysis,sep='\t',index=False)			
			
		
		#do the same as above, except for J genes this time 
		if 'jgene' in statistics_to_run:
			print('Performing J gene analysis')
			j_gene_analysis = output_file_prefix+'.jgenes.txt'
			sorted_j_counts =  gene_df.groupby(['recomb','JGENE']).sum()#.sort('VGENE',ascending=1)						
			jgene_level = sorted_j_counts.index.names.index('JGENE')			
			if '' in list(sorted_j_counts.index.levels[jgene_level]):
				sorted_j_counts.drop('',level='JGENE',inplace=True)			
			ignore_counts = ['TOTAL_COUNTS','VGENE']
			keep_col = [n for n in sorted_j_counts.columns if n not in ignore_counts]
			g = sorted_j_counts[keep_col]			
			sorted_j_counts.reset_index().sort(['recomb','TOTAL_COUNTS'],ascending=[1,0]).iloc[:,:-1].to_csv(j_gene_analysis,sep='\t',index=False)
			
			#NOW CALCULATE FREQUENCY for every exeprement 						
			if 'VDJ' in list(g.index.levels[g.index.names.index('recomb')]):
				vdj_g = g.xs('VDJ',level='recomb')
				PlotGeneDist(vdj_g/vdj_g.sum(),gene_analysis_plot+'.vdj.jgenes','JH Gene Distribution','Frequency','J Gene',max_val=None,min_val=0,step=5)
				plots_created.append(gene_analysis_plot+'.vdj.jgenes.png') #.png extension is added in the function plotgenedist
			if 'VJ' in list(g.index.levels[g.index.names.index('recomb')]):
				vj_g = g.xs('VJ',level='recomb')			
				PlotGeneDist(vj_g/vj_g.sum(),gene_analysis_plot+'.vj.jgenes','JL Gene Distribution','Frequency','J Gene',max_val=None,min_val=0,step=5)
				plots_created.append(gene_analysis_plot+'.vj.jgenes.png') #.png extension is added in the function plotgenedist		
		
		#now perform a V-J gene analysis (heat map) for each experiment 
		if 'vjgene' in statistics_to_run:
			print('Performing V-J gene analysis')
			vj_gene_analysis = output_file_prefix+'.v_and_jgene_analysis.txt'
			#group datafraom by recombination, vgene, and jgene 
			#first rename all V and J gnees that are empyt as No call 						
			#Then Group H / L results by  by v and j gnees and take the sum of each column in the group 
			vj_df =  gene_df.replace([''],[' No call']).groupby(['recomb','VGENE','JGENE']).sum()
			vj_df.to_csv(vj_gene_analysis,sep='\t')			
			
			#remove TOTAL_COUNTS			
			vj_df.drop('TOTAL_COUNTS', axis=1, inplace=True)
			
			#calculate frequency for each recomb type 
			if 'VDJ' in list(vj_df.index.levels[g.index.names.index('recomb')]):						
				v1 =  vj_df.loc['VDJ',:]/vj_df.loc['VDJ',:].sum()
				PlotVJGeneHeatMap(v1,gene_analysis_plot+'.vdj.v_and_jgene_analysis',max_val=None,min_val=None)
				plots_created.append(gene_analysis_plot+'.vdj.v_and_jgene_analysis.png') #.png extension is added in the function plotgenedist
			if 'VJ' in list(vj_df.index.levels[g.index.names.index('recomb')]):
				v2 =  vj_df.loc['VJ',:]/vj_df.loc['VJ',:].sum()
				PlotVJGeneHeatMap(v2,gene_analysis_plot+'.vj.v_and_jgene_analysis',max_val=None,min_val=None)																							
				plots_created.append(gene_analysis_plot+'.vj.v_and_jgene_analysis.png') #.png extension is added in the function plotgenedist
			del vj_df
		del gene_df
		
	#lets do some cdr3 analysis 									
	cdr3_length_stats = {}
	diversity_measurements = {}
	if cdr3analysis:	
		unique_cdr3_file = output_file_prefix+'.unique_cdr3_counts.txt' 
		print('Performing CDR3 analyisis')
		if sum(num_cdr3)>0:
			#again create a pandas dataframe but this time using the unique cdr3 calls 
			print('Loading CDR3s into a dataframe')
			cdr3_df_list = [pd.DataFrame.from_dict(c,orient='index') for c in [cdr3_dict_vdj,cdr3_dict_vj,cdr3_dict_unk]]
			#merge all dftogether
			keys=['VDJ','VJ','UNK']
			cdr3_df = pd.concat(cdr3_df_list,keys=keys)
			#cdr3_df = pd.DataFrame(cdr3_dict).transpose()			
			cdr3_df['TOTAL_COUNTS'] = cdr3_df.sum(axis=1)
			print('Dataframe created')
			
			cdr3_df.index.names = ['recomb','CDR3']
			cdr3_df = cdr3_df.reset_index()				
			#cdr3_df['CDR3'] = ''
			#cdr3_df['recomb'] = ''
			#cdr3_df = cdr3_df.apply(ModifyPDTable,axis=1,raw=True,reduce=True,args=(['CDR3','recomb'],delim))			
			
			new_names = {}
			#performm 			
			cdr3_df['CDR3_LENGTH'] = cdr3_df.CDR3.map(len) 
			for f,v in enumerate(exp_names):
				new_names[f]=v
			#rename the columns to match the experiment names 
			
			cdr3_df = cdr3_df.rename(columns=new_names)
			cdr3_df.sort(['recomb','TOTAL_COUNTS'],ascending=[1,0],inplace=True)
			cdr3_df.set_index(['recomb','CDR3'],inplace=True)					
			
			#save dataframe as tab dleim file 						
			cdr3_df.to_csv(unique_cdr3_file,sep='\t')									
			
			cdr3_length_stats = PlotCDR3Histogram(cdr3_df,gene_analysis_plot+'.cdr3_length_histogram')
			plots_created.append(gene_analysis_plot+'.cdr3_length_histogram.png')
			
			diversity_measurements = CalculateDiversities(cdr3_df,gene_analysis_plot+'.cdr3_diversity_plots')
			plots_created.append(gene_analysis_plot+'.cdr3_diversity_plots.png')			
		del cdr3_df
	
	print('Writing summary to file')
	#finally make a results text file that summarizes all the information	
	GenerateResultsSummaryFile(gene_summary_file,statistics_to_run,list_of_files,exp_names,unique_aa_file,unique_cdr3_file,v_gene_analysis,j_gene_analysis,vj_gene_analysis,plots_created,num_sequences,num_results,num_vdj,num_vj,num_cdr3,num_stop_codon,cdr3_length_stats,diversity_measurements)	
	
	files_generated = [gene_summary_file]
	if unique_aa_file:
		files_generated.append(unique_aa_file)
	if unique_cdr3_file:
		files_generated.append(unique_cdr3_file)
	if v_gene_analysis:
		files_generated.append(v_gene_analysis)
	if j_gene_analysis:
		files_generated.append(j_gene_analysis)
	if vj_gene_analysis:
		files_generated.append(vj_gene_analysis)
	
	print('Descriptive statistics completed at {0}.'.format(str(datetime.datetime.now())))
	
	gc.collect()

	
	return {'files':files_generated,'figures':plots_created}

def GenerateResultsSummaryFile(output_file,analysis_requests,input_file_paths,exp_names,aa_file_location,cdr3_file_location,v_gene_file_location,j_gene_file_location,vj_gene_file_location,plots_created,num_sequences,num_results,num_vdj,num_vj,num_cdr3,num_stop_codon,cdr3_length_stats,diversity_measurements):
	
	with open(output_file,'w') as summary:
		summary.write('Data generated on: %s\n' %(datetime.datetime.now()))		
		summary.write("******************************************************************************\n")
		summary.write("******************Summary Report for Descriptive Statistics*******************\n")
		summary.write("******************************************************************************\n\n\n")
		summary.write('The following experiments were used in the analysis:\n')
		summary.write('\tExperiment Name\tFile path\n')
		exp_string = ""
		for num,exp in enumerate(exp_names):
			exp_string += '\t'+exp
			if isinstance(input_file_paths[num],list):
				exp_string+='\t'+input_file_paths[num][0]+'\n'
				for k in range(1,len(input_file_paths[num][0])):
					exp_string+='\t\t'+input_file_paths[num][k]+'\n'
			else:
				exp_string+='\t'+input_file_paths[num]+'\n'
		summary.write(exp_string+'\n\n')
		summary.write('The following analyses were requested: {0}'.format(','.join(analysis_requests))+'\n\n')
		
		total_unique_aa_count = 0
		if aa_file_location:
			total_unique_aa_count = useful.file_line_count(aa_file_location)-1				
				
		summary.write('General File Info\n')
		summary.write('\t'+'\t'.join(['Experiment name','Number of results found','Number sequences containing antibody sequence','Number sequences with stop codon','Num VDJ','Num VJ','Number sequences with CDR3 amino acid sequence'])+'\n')
		for ind,each_exp in enumerate(exp_names):
			num_total_seqs = num_sequences[ind]			
			res_string = str(num_results[ind])+' ('+str(round(100*float(num_results[ind])/num_total_seqs,3))+'%)'
			stop_string = str(num_stop_codon[ind])+' ('+str(round(100*float(num_stop_codon[ind])/num_total_seqs,3))+'%)'
			summary.write('\t'+'\t'.join([each_exp,str(num_total_seqs),res_string,stop_string,str(num_vdj[ind]),str(num_vj[ind]),str(num_cdr3[ind])])+'\n')
		if len(exp_names)>1:
			num_total_seqs = sum(num_sequences)
			res_string = str(sum(num_results))+' ('+str(round(100*float(sum(num_results))/num_total_seqs,3))+'%)'
			stop_string = str(sum(num_stop_codon))+' ('+str(round(100*float(sum(num_stop_codon))/num_total_seqs,3))+'%)'
			summary.write('\t'+'\t'.join(['Total',str(num_total_seqs),res_string,stop_string,str(sum(num_vdj)),str(sum(num_vj)),str(sum(num_cdr3))])+'\n')
		if total_unique_aa_count:
			summary.write('A total of {0} unique amino acid antibody sequences were found across all files\n'.format(str(total_unique_aa_count)))
		
		if 'cdr3' in analysis_requests:								
			summary.write('\nCDR3 Results Summary\n')			
			summary.write('*The average CDR3 length only considers CDR3 amino acid sequences with more than 2 amino acids\n')
			use_names = exp_names
			if len(exp_names)>1:
				use_names.append('TOTAL_COUNTS')					
			for rtype in ['VDJ','VJ']:				
				if rtype in cdr3_length_stats or rtype in diversity_measurements:
					summary.write('CDR3 - '+rtype+'\n')			
					summary.write('\t'+'\t'.join(['Experiment name','Average CDR3 length*','Standard deviation CDR3 length','Number unique CDR3 AA','Number unique CDR3 AA identified more than two times','Shannon diversity index', 'Ginni-simpsons diversity index','Normalized shannon diversity','Normalized ginni-simpsons diversity'])+'\n')									
					
					for exp_num,each_exp in enumerate(use_names):															
						if each_exp == 'TOTAL_COUNTS':
							results = ['All experiments']
						else:
							results = [each_exp]																						
						if each_exp in cdr3_length_stats[rtype]['mean']:
							results.extend([str(cdr3_length_stats[rtype]['mean'][each_exp]),str(cdr3_length_stats[rtype]['std'][each_exp])])
						else:
							results.extend(['',''])																		
						try:
							results.append(str(diversity_measurements[rtype]['unique_cdr3s'][each_exp]))
						except:
							results.append('')							
						try:
							results.append(str(diversity_measurements[rtype]['num_above_2'][each_exp]))
						except:
							results.append('')						
						try:
							results.append(str(diversity_measurements[rtype]['shannon_entropy']['index'][each_exp]))
						except:
							results.append('')						
						try:
							results.append(str(diversity_measurements[rtype]['ginni_simpsons']['index'][each_exp]))
						except:
							results.append('')						
						try:
							results.append(str(diversity_measurements[rtype]['shannon_entropy']['true_diversity'][each_exp]))
						except:
							results.append('')						
						try:
							results.append(str(diversity_measurements[rtype]['ginni_simpsons']['true_diversity'][each_exp]))
						except:
							results.append('')																		
						summary.write('\t'+'\t'.join(results)+'\n')																																		
				summary.write('\n')
		
		summary.write('Result files created from program\n')
		if aa_file_location:
			summary.write('\tA file containing the unique antibody amino acid sequences and counts can be found in the following location:\n\t\t{0}\n'.format(aa_file_location)) 			
		if cdr3_file_location:
			summary.write('\tA file containing a list of unique cdr3 amino acid sequences and counts can be found in the following location:\n\t\t{0}\n'.format(cdr3_file_location)) 
		if v_gene_file_location:
			summary.write('\tA file summarizing VGENE usage can be found in the following location:\n\t\t{0}\n'.format(v_gene_file_location))
		if j_gene_file_location:
			summary.write('\tA file summarizing JGENE usage can be found in the following location:\n\t\t{0}\n'.format(j_gene_file_location))
		if vj_gene_file_location:
			summary.write('\tA file summarizing V and JGENE usage can be found in the following location:\n\t\t{0}\n'.format(vj_gene_file_location))
		if len(plots_created)>0:
			summary.write('\tThe following figures in both PNG and SVG format were generated during the analysis:\n')
			for fig in plots_created:
				summary.write('\t\t'+fig+'\n')		
				
#if 'vgene','jgene','vjgene','cdr3','diversity' 
def ModifyPDTable(row,new_col_names,delim):			
	for i,c in enumerate(row.pop('index').split(delim)):				
		row[new_col_names[i]] = c
	return row
	

def GenerateAAFile(intermediate_file,unique_aa_file,header_row,experiment_names):
	#we use bash for calling efficient functions (sort and awk) for processing of a TAB delimited file 
	#assume the last column of the tab file corresponds to COUNTS for a SPECIFIC EPXERIMENT NUMBER (it has to be a number corersponding to an experiment number)
	#these experiment numbers should correlate to indexes in experiment_names variable
	
	
	total_count_column = len(header_row) #the column for total counts will always be the last column in the header_row variable. 
	#THEN we have to add counts for each experiment 
	for e in experiment_names:
		header_row.append(e+' Counts')
	num_exp_counts = len(experiment_names)	
	hfile = intermediate_file+'.header.txt'
	parent = useful.get_parent_dir(intermediate_file) # os.path.dirname(os.path.dirname(os.path.abspath(intermediate_file)))# '/'.join(intermediate_file.split('/')[:-1])+'/'
	#write a header row 
	with open(hfile,'w') as w:
		w.write('\t'.join(header_row)+'\n')
	
	#the bash script process_aa_file will sort sequences in file then count their occurrences and collapse	
	this_script_folder = useful.get_parent_dir(__file__)
	#the bash script,process_unique_ab_aa_file.bash, MUST be in the same folder as this script (immunogrep_descriptive_statistics)
	bash_script_path = os.path.join(this_script_folder,'process_unique_ab_aa_file.bash')
	bash_command = '''bash '{5}' '{0}' '{1}' '{2}' {3} {4}'''.format(intermediate_file,unique_aa_file,hfile,str(num_exp_counts),str(total_count_column),bash_script_path)
	
	#run the bash 
	#subprocess.call(bash_command,shell=True)
	subprocess.call(bash_command,shell=True)
	
