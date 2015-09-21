from immunogrep_global_variables import idIdentifier

	
def grab_translator_variable(analysis_name):
	"""
		This is a simple function for returning the proper dictionary for defining a translator for inserting annotated data to the database. 
		The translator defines how fields from specific annotation files should be stored in the database. It simply maps the fieldnames 
		in a specific file to the field names int he database. This function will return the hardcoded translator for our currently supported
		analysis names in the pipeline 
		
		We currently define the fields for: IMGT, IGBLAST, IGFFT, MIXCR

	"""

	all_translators = {
		'MIXCR':{},
		'IGFFT':{},
		'IMGT':{#IMGT FILE 1 WILL ALWAYS BE REQUIRED                                                 
			"ANALYSIS_NAME":'IMGT',
			"RECOMBINATION_FIELD":{		
				#LAMBDA FUNCTION FOR RECOMBINATION FIELD IS A HACK TO ENSURE THAT, AT A MINIMUM, ONLY IMGT FILE ONE IS REQUIRED (However fieldnames from other files will not be added to the database)
				"LAMBDA_FUNCTION":"lambda x: 'VDJ' if x['V-D-J-REGION_3'] else ('VJ' if x['V-J-REGION_3'] else ('' if not(x['V-GENE and allele_1']) else ('VDJ' if ('IGHV' in x['V-GENE and allele_1'] or 'TRBV' in x['V-GENE and allele_1']) else 'VJ')))",
			},
			"FIELDS_UPDATE":{ #updates these field into database. key = > database fields, value = > IMGT field names
				idIdentifier:idIdentifier,
				'COMMAND':'COMMAND',
				'SEQUENCE_HEADER':'Sequence ID_1',		
				'SEQUENCE':'Sequence_1',
				'STRAND':'Orientation_1',
				'PRODUCTIVE':'Functionality_1',
				'VREGION.VGENE_SCORES':'V-REGION score_1',		
				'JREGION.JGENE_SCORES':'J-REGION score_1',		
				'NOTES':'Functionality comment_1',
		
				'VREGION.CDR1.NT':'CDR1-IMGT_3',
				'GAPPED.VREGION.CDR1.NT':'CDR1-IMGT_2',
		
				'VREGION.CDR1.AA':'CDR1-IMGT_5',
				'GAPPED.VREGION.CDR1.AA':'CDR1-IMGT_4',
		
				'VREGION.CDR2.NT':'CDR2-IMGT_3',
				'GAPPED.VREGION.CDR2.NT':'CDR2-IMGT_2',
		
				'VREGION.CDR2.AA':'CDR2-IMGT_5',
				'GAPPED.VREGION.CDR2.AA':'CDR2-IMGT_4',
		
				'CDR3.NT':'CDR3-IMGT_3',		
				'GAPPED.CDR3.NT':'CDR3-IMGT_2',				
		
				'CDR3.AA':'CDR3-IMGT_5',
				'GAPPED.CDR3.AA':'CDR3-IMGT_4',		
		
				'VREGION.FR1.NT':'FR1-IMGT_3',
				'GAPPED.VREGION.FR1.NT':'FR1-IMGT_2',
		
				'VREGION.FR1.AA':'FR1-IMGT_5',
				'GAPPED.VREGION.FR1.AA':'FR1-IMGT_4',
		
				'VREGION.FR2.NT':'FR2-IMGT_3',
				'GAPPED.VREGION.FR2.NT':'FR2-IMGT_2',
		
				'VREGION.FR2.AA':'FR2-IMGT_5',
				'GAPPED.VREGION.FR2.AA':'FR2-IMGT_4',
				
				'VREGION.FR3.NT':'FR3-IMGT_3',
				'GAPPED.VREGION.FR3.NT':'FR3-IMGT_2',
		
				'VREGION.FR3.AA':'FR3-IMGT_5',	
				'GAPPED.VREGION.FR3.AA':'FR3-IMGT_4',
		
				'JREGION.FR4.NT':'FR4-IMGT_3',		
				'GAPPED.JREGION.FR4.NT':'FR4-IMGT_2',
		
				'JREGION.FR4.AA':'FR4-IMGT_5',								
				'GAPPED.JREGION.FR4.AA':'FR4-IMGT_4',								
		
				'VREGION.VGENE_QUERY_START':'V-REGION start_3',
				'VREGION.VGENE_QUERY_END':'V-REGION end_3',
				'JREGION.JGENE_QUERY_START':'J-REGION start_3',
				'JREGION.JGENE_QUERY_END':'J-REGION end_3',
				'VREGION.VGENES':'V-GENE and allele_1',
				'JREGION.JGENES':'J-GENE and allele_1',
				'DREGION.DGENES':'D-GENE and allele_1',
		
				'VREGION.SHM.NT':'VREGION.SHM.NT', #created in filereader class (uses the field from fithe file summary V_REGION IDENTITY)
				'VREGION.SHM.NT_PER':'VREGION.SHM.NT_PER', #created in filereader class 
				'VREGION.SHM.AA':'VREGION.SHM.AA', #created in filereader class
				'VREGION.SHM.AA_PER':'VREGION.SHM.AA_PER', #created in filereader class
		
				'JREGION.SHM.NT':'JREGION.SHM.NT', #created in filereader class (uses the field from fithe file summary J_REGION IDENTITY)
				'JREGION.SHM.NT_PER':'JREGION.SHM.NT_PER', #created in filereader class 
		
				"PREDICTED_AB_SEQ.NT":"PREDICTED_AB_SEQ.NT",#created in filereader class
				"PREDICTED_AB_SEQ.AA":"PREDICTED_AB_SEQ.AA"#created in filereader class				
			}
		},		
		'IGBLAST':{}
	}

	return all_translators[analysis_name]