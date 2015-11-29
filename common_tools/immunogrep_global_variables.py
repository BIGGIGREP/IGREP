#########################################################################################################
# Please place the name of the class/script that imports this script and uses these global variables:    
#  immunogrep_immunogrepfile.py
#
#
#
#

# +++++++++++++++++++
# +Global variables +
# +++++++++++++++++++

import os

program_folder = os.path.join(os.path.dirname(os.path.relpath(__file__)), "program_binaries")

descriptor_symbol='#@'   #To signify the line has information needed to transact in immunogrep pipeline (PANDAS does not like comments > 1 char)
fasta_file_delimiter=' <' #To signify additional info in sequence header
idIdentifier='SEQ_ID'        #for use with MongoDB
expIdentifier='EXP_ID'
seqRawData = '@SEQ' #this is the ANALYSIS_NAME provided to mongo documents corresponding to raw nucleotide sequeences 

# Note, an example of the decoratorline would be as follows:
# #@D{'DESCRIPTION':[list of label headers; order matters here],'TRANSLATON':{dictionary for communicating with the database}}
description_var='DESCRIPTION'
translation_var='TRANSLATION'
filetype_var = 'FILETYPE'
tracktion_var='TRACKTION'
type_var='TYPE' # This is the analysis type key (i.e. 'TYPE':'IMGT' or 'TYPE':'INHOUSE', etc.) 

# File extensions in immunogrep pipeline
annotationExt = 'annotation'
analysisExt = 'analysis'
queryExt = 'query'
IFExt = 'IFFile'
listofextension = ['fasta', 'fastq', 'fq', 'fna', 'txt', 'csv', 'fas', 'fa', 'json', 'query', 'annotation', 'analysis', 'IFFile']		


imgt_files_info = [ #this summarizes each of the column names associated with each imgt file from IMGT V-QUEST
    {
        'Name':  '1_Summary',                                    
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION score','V-REGION identity %',
            'V-REGION identity nt','V-REGION identity % (with ins/del events)','V-REGION identity nt (with ins/del events)',
            'J-GENE and allele','J-REGION score','J-REGION identity %','J-REGION identity nt','D-GENE and allele',
            'D-REGION reading frame','CDR1-IMGT length','CDR2-IMGT length','CDR3-IMGT length','CDR-IMGT lengths','FR-IMGT lengths',
            'AA JUNCTION','JUNCTION frame','Orientation','Functionality comment','V-REGION potential ins/del','J-GENE and allele comment',
            'V-REGION insertions','V-REGION deletions','Sequence'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION score','V-REGION identity %',
            'V-REGION identity nt','V-REGION identity % (with ins/del events)','V-REGION identity nt (with ins/del events)',
            'J-GENE and allele','J-REGION score','J-REGION identity %','J-REGION identity nt','D-GENE and allele',
            'D-REGION reading frame','CDR1-IMGT length','CDR2-IMGT length','CDR3-IMGT length','CDR-IMGT lengths','FR-IMGT lengths',
            'AA JUNCTION','JUNCTION frame','Orientation','Functionality comment','V-REGION potential ins/del','J-GENE and allele comment',
            'V-REGION insertions','V-REGION deletions','Sequence'
        ]             
    },
    {
        'Name':  '2_IMGT-gapped-nt-sequences',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ]                    
    },            
    {                         
        'Name':  '3_Nt-sequences',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele','V-D-J-REGION',
            'V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT','JUNCTION',"3'V-REGION",
            '(N-D)-J-REGION','(N-D)-REGION',"P3'V",'N-REGION','N1-REGION',"P5'D",'D-REGION',"P3'D","P5'D1",'D1-REGION',"P3'D1",
            'N2-REGION',"P5'D2",'D2-REGION',"P3'D2",'N3-REGION',"P5'D3",'D3-REGION',"P3'D3",'N4-REGION',"P5'J","5'J-REGION",'D-J-REGION',
            'J-REGION','FR4-IMGT','V-D-J-REGION start','V-D-J-REGION end','V-J-REGION start','V-J-REGION end','V-REGION start','V-REGION end',
            'FR1-IMGT start','FR1-IMGT end','CDR1-IMGT start','CDR1-IMGT end','FR2-IMGT start','FR2-IMGT end','CDR2-IMGT start','CDR2-IMGT end',
            'FR3-IMGT start','FR3-IMGT end','CDR3-IMGT start','CDR3-IMGT end','JUNCTION start','JUNCTION end',"3'V-REGION start","3'V-REGION end",
            '(N-D)-J-REGION start','(N-D)-J-REGION end','(N-D)-REGION start','(N-D)-REGION end',"P3'V start","P3'V end",'N-REGION start',
            'N-REGION end','N1-REGION start','N1-REGION end',"P5'D start","P5'D end",'D-REGION start','D-REGION end',"P3'D start","P3'D end",
            "P5'D1 start","P5'D1 end",'D1-REGION start','D1-REGION end',"P3'D1 start","P3'D1 end",'N2-REGION start','N2-REGION end',
            "P5'D2 start","P5'D2 end",'D2-REGION start','D2-REGION end',"P3'D2 start","P3'D2 end",'N3-REGION start','N3-REGION end',
            "P5'D3 start","P5'D3 end",'D3-REGION start','D3-REGION end',"P3'D3 start","P3'D3 end",'N4-REGION start','N4-REGION end',
            "P5'J start","P5'J end","5'J-REGION start","5'J-REGION end",'D-J-REGION start','D-J-REGION end','J-REGION start','J-REGION end',
            'FR4-IMGT start','FR4-IMGT end'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele','V-D-J-REGION',
            'V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT','JUNCTION', 'FR4-IMGT','V-D-J-REGION start',
            'V-D-J-REGION end','V-J-REGION start','V-J-REGION end','V-REGION start','V-REGION end','FR1-IMGT start','FR1-IMGT end','CDR1-IMGT start',
            'CDR1-IMGT end','FR2-IMGT start','FR2-IMGT end','CDR2-IMGT start','CDR2-IMGT end','FR3-IMGT start','FR3-IMGT end','CDR3-IMGT start',
            'CDR3-IMGT end','JUNCTION start', 'JUNCTION end','D-J-REGION start','D-J-REGION end','J-REGION start','J-REGION end','FR4-IMGT start','FR4-IMGT end'
        ]                           
    },
    {
        'Name':  '4_IMGT-gapped-AA-sequences',                          
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ]
    },             
    {
        'Name':  '5_AA-sequences',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        ]
    },      
    {
        'Name':  '6_Junction',        
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'JUNCTION frame','JUNCTION','JUNCTION (with frameshift)',"3'V-REGION","P3'V",'N-REGION','N1-REGION',"P5'D",
            'D-REGION',"P3'D","P5'D1",'D1-REGION',"P3'D1",'N2-REGION',"P5'D2",'D2-REGION',"P3'D2",'N3-REGION',"P5'D3",
            'D3-REGION',"P3'D3",'N4-REGION',"P5'J","5'J-REGION",'JUNCTION-nt nb',"3'V-REGION-nt nb","P3'V-nt nb",'N-REGION-nt nb',
            'N1-REGION-nt nb',"P5'D-nt nb",'D-REGION-nt nb',"P3'D-nt nb","P5'D1-nt nb",'D1-REGION-nt nb',"P3'D1-nt nb",
            'N2-REGION-nt nb',"P5'D2-nt nb",'D2-REGION-nt nb',"P3'D2-nt nb",'N3-REGION-nt nb',"P5'D3-nt nb",'D3-REGION-nt nb',
            "P3'D2-nt nb",'N4-REGION-nt nb',"P5'J-nt nb","5'J-REGION-nt nb","3'V-REGION trimmed-nt nb","5'D-REGION trimmed-nt nb",
            "3'D-REGION trimmed-nt nb","5'D1-REGION trimmed-nt nb","3'D1-REGION trimmed-nt nb","5'D2-REGION trimmed-nt nb",
            "3'D2-REGION trimmed-nt nb","5'D3-REGION trimmed-nt nb","3'D3-REGION trimmed-nt nb","5'J-REGION trimmed-nt nb",
            "3'V-REGION mut-nt nb",'D-REGION mut-nt nb','D1-REGION mut-nt nb','D2-REGION mut-nt nb','D3-REGION mut-nt nb',
            "5'J-REGION mut-nt nb",'D-REGION reading frame','Ngc','CDR3-IMGT length','Molecular mass','pI',"3'V-REGION accepted mut nb",
            'D-REGION accepted mut nb',"5'J-REGION accepted mut nb",'Accepted D-GENE nb','CDR3-IMGT','CDR3-IMGT-nt nb','CDR3-IMGT (with frameshift)',
            'CDR3-IMGT (AA)','CDR3-IMGT (AA) (with frameshift)','JUNCTION (AA)','JUNCTION (AA) (with frameshift)'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'JUNCTION frame','JUNCTION','JUNCTION (with frameshift)','CDR3-IMGT','CDR3-IMGT-nt nb','CDR3-IMGT (with frameshift)',
            'CDR3-IMGT (AA)','CDR3-IMGT (AA) (with frameshift)','JUNCTION (AA)','JUNCTION (AA) (with frameshift)'
        ]
    },
    {
        'Name':  '7_V-REGION-mutation-and-AA-change-table',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT',
            'CDR2-IMGT','FR3-IMGT','CDR3-IMGT'
        ],
        'Important_Headers':[
            'Sequence number','Sequence ID','Functionality','V-GENE and allele'
        ],
    },
    {
        'Name':  '8_V-REGION-nt-mutation-statistics',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION Nb of positions including IMGT gaps (nt)',
            'V-REGION Nb of nucleotides','V-REGION Nb of identical nucleotides','V-REGION Nb of mutations','V-REGION Nb of silent mutations',
            'V-REGION Nb of nonsilent mutations','V-REGION a>g','V-REGION g>a','V-REGION c>t','V-REGION t>c','V-REGION a>c',
            'V-REGION c>a','V-REGION a>t','V-REGION t>a','V-REGION g>c','V-REGION c>g','V-REGION g>t','V-REGION t>g',
            'FR1-IMGT Nb of positions including IMGT gaps (nt)','FR1-IMGT Nb of nucleotides','FR1-IMGT Nb of identical nucleotides',
            'FR1-IMGT Nb of mutations','FR1-IMGT Nb of silent mutations','FR1-IMGT Nb of nonsilent mutations','FR1-IMGT a>g',
            'FR1-IMGT g>a','FR1-IMGT c>t','FR1-IMGT t>c','FR1-IMGT a>c','FR1-IMGT c>a','FR1-IMGT a>t','FR1-IMGT t>a',
            'FR1-IMGT g>c','FR1-IMGT c>g','FR1-IMGT g>t','FR1-IMGT t>g','CDR1-IMGT Nb of positions including IMGT gaps (nt)',
            'CDR1-IMGT Nb of nucleotides','CDR1-IMGT Nb of identical nucleotides','CDR1-IMGT Nb of mutations','CDR1-IMGT Nb of silent mutations',
            'CDR1-IMGT Nb of nonsilent mutations','CDR1-IMGT a>g','CDR1-IMGT g>a','CDR1-IMGT c>t','CDR1-IMGT t>c','CDR1-IMGT a>c',
            'CDR1-IMGT c>a','CDR1-IMGT a>t','CDR1-IMGT t>a','CDR1-IMGT g>c','CDR1-IMGT c>g','CDR1-IMGT g>t','CDR1-IMGT t>g',
            'FR2-IMGT Nb of positions including IMGT gaps (nt)','FR2-IMGT Nb of nucleotides','FR2-IMGT Nb of identical nucleotides',
            'FR2-IMGT Nb of mutations','FR2-IMGT Nb of silent mutations','FR2-IMGT Nb of nonsilent mutations','FR2-IMGT a>g',
            'FR2-IMGT g>a','FR2-IMGT c>t','FR2-IMGT t>c','FR2-IMGT a>c','FR2-IMGT c>a','FR2-IMGT a>t','FR2-IMGT t>a','FR2-IMGT g>c',
            'FR2-IMGT c>g','FR2-IMGT g>t','FR2-IMGT t>g','CDR2-IMGT Nb of positions including IMGT gaps (nt)','CDR2-IMGT Nb of nucleotides',
            'CDR2-IMGT Nb of identical nucleotides','CDR2-IMGT Nb of mutations','CDR2-IMGT Nb of silent mutations','CDR2-IMGT Nb of nonsilent mutations',
            'CDR2-IMGT a>g','CDR2-IMGT g>a','CDR2-IMGT c>t','CDR2-IMGT t>c','CDR2-IMGT a>c','CDR2-IMGT c>a','CDR2-IMGT a>t',
            'CDR2-IMGT t>a','CDR2-IMGT g>c','CDR2-IMGT c>g','CDR2-IMGT g>t','CDR2-IMGT t>g','FR3-IMGT Nb of positions including IMGT gaps (nt)',
            'FR3-IMGT Nb of nucleotides','FR3-IMGT Nb of identical nucleotides','FR3-IMGT Nb of mutations','FR3-IMGT Nb of silent mutations',
            'FR3-IMGT Nb of nonsilent mutations','FR3-IMGT a>g','FR3-IMGT g>a','FR3-IMGT c>t','FR3-IMGT t>c','FR3-IMGT a>c',
            'FR3-IMGT c>a','FR3-IMGT a>t','FR3-IMGT t>a','FR3-IMGT g>c','FR3-IMGT c>g','FR3-IMGT g>t','FR3-IMGT t>g','CDR3-IMGT Nb of positions including IMGT gaps (nt)',
            'CDR3-IMGT Nb of nucleotides','CDR3-IMGT Nb of identical nucleotides','CDR3-IMGT Nb of mutations','CDR3-IMGT Nb of silent mutations',
            'CDR3-IMGT Nb of nonsilent mutations','CDR3-IMGT a>g','CDR3-IMGT g>a','CDR3-IMGT c>t','CDR3-IMGT t>c','CDR3-IMGT a>c',
            'CDR3-IMGT c>a', 'CDR3-IMGT a>t','CDR3-IMGT t>a','CDR3-IMGT g>c','CDR3-IMGT c>g','CDR3-IMGT g>t','CDR3-IMGT t>g'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION Nb of positions including IMGT gaps (nt)',
            'V-REGION Nb of nucleotides','V-REGION Nb of identical nucleotides','V-REGION Nb of mutations','V-REGION Nb of silent mutations',
            'V-REGION Nb of nonsilent mutations','FR1-IMGT Nb of positions including IMGT gaps (nt)','FR1-IMGT Nb of nucleotides','FR1-IMGT Nb of identical nucleotides',
            'FR1-IMGT Nb of mutations','FR1-IMGT Nb of silent mutations','FR1-IMGT Nb of nonsilent mutations','CDR1-IMGT Nb of nucleotides','CDR1-IMGT Nb of identical nucleotides','CDR1-IMGT Nb of mutations','CDR1-IMGT Nb of silent mutations',
            'CDR1-IMGT Nb of nonsilent mutations','FR2-IMGT Nb of positions including IMGT gaps (nt)','FR2-IMGT Nb of nucleotides','FR2-IMGT Nb of identical nucleotides',
            'FR2-IMGT Nb of mutations','FR2-IMGT Nb of silent mutations','FR2-IMGT Nb of nonsilent mutations','CDR2-IMGT Nb of positions including IMGT gaps (nt)','CDR2-IMGT Nb of nucleotides',
            'CDR2-IMGT Nb of identical nucleotides','CDR2-IMGT Nb of mutations','CDR2-IMGT Nb of silent mutations','CDR2-IMGT Nb of nonsilent mutations','FR3-IMGT Nb of positions including IMGT gaps (nt)',
            'FR3-IMGT Nb of nucleotides','FR3-IMGT Nb of identical nucleotides','FR3-IMGT Nb of mutations','FR3-IMGT Nb of silent mutations',
            'FR3-IMGT Nb of nonsilent mutations','CDR3-IMGT Nb of positions including IMGT gaps (nt)',
            'CDR3-IMGT Nb of nucleotides','CDR3-IMGT Nb of identical nucleotides','CDR3-IMGT Nb of mutations','CDR3-IMGT Nb of silent mutations',
            'CDR3-IMGT Nb of nonsilent mutations'
        ],
    },
    {
        'Name':  '9_V-REGION-AA-change-statistics',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION Nb of positions including IMGT gaps (AA)',
            'V-REGION Nb of AA','V-REGION Nb of identical AA','V-REGION Nb of AA changes','V-REGION +++','V-REGION ++-',
            'V-REGION +-+','V-REGION +--','V-REGION -+-','V-REGION --+','V-REGION ---','V-REGION Very similar','V-REGION Similar',
            'V-REGION Dissimilar','V-REGION Very dissimilar','FR1-IMGT Nb of positions including IMGT gaps (AA)','FR1-IMGT Nb of AA',
            'FR1-IMGT Nb of identical AA','FR1-IMGT Nb of AA changes','FR1-IMGT +++','FR1-IMGT ++-','FR1-IMGT +-+','FR1-IMGT +--',
            'FR1-IMGT -+-','FR1-IMGT --+','FR1-IMGT ---','FR1-IMGT Very similar','FR1-IMGT Similar','FR1-IMGT Dissimilar','FR1-IMGT Very dissimilar',
            'CDR1-IMGT Nb of positions including IMGT gaps (AA)','CDR1-IMGT Nb of AA','CDR1-IMGT Nb of identical AA','CDR1-IMGT Nb of AA changes',
            'CDR1-IMGT +++','CDR1-IMGT ++-','CDR1-IMGT +-+','CDR1-IMGT +--','CDR1-IMGT -+-','CDR1-IMGT --+','CDR1-IMGT ---',
            'CDR1-IMGT Very similar','CDR1-IMGT Similar','CDR1-IMGT Dissimilar','CDR1-IMGT Very dissimilar','FR2-IMGT Nb of positions including IMGT gaps (AA)',
            'FR2-IMGT Nb of AA','FR2-IMGT Nb of identical AA','FR2-IMGT Nb of AA changes','FR2-IMGT +++','FR2-IMGT ++-','FR2-IMGT +-+',
            'FR2-IMGT +--','FR2-IMGT -+-','FR2-IMGT --+','FR2-IMGT ---','FR2-IMGT Very similar','FR2-IMGT Similar','FR2-IMGT Dissimilar',
            'FR2-IMGT Very dissimilar','CDR2-IMGT Nb of positions including IMGT gaps (AA)','CDR2-IMGT Nb of AA','CDR2-IMGT Nb of identical AA',
            'CDR2-IMGT Nb of AA changes','CDR2-IMGT +++','CDR2-IMGT ++-','CDR2-IMGT +-+','CDR2-IMGT +--','CDR2-IMGT -+-','CDR2-IMGT --+',
            'CDR2-IMGT ---','CDR2-IMGT Very similar','CDR2-IMGT Similar','CDR2-IMGT Dissimilar','CDR2-IMGT Very dissimilar',
            'FR3-IMGT Nb of positions including IMGT gaps (AA)','FR3-IMGT Nb of AA','FR3-IMGT Nb of identical AA','FR3-IMGT Nb of AA changes',
            'FR3-IMGT +++','FR3-IMGT ++-','FR3-IMGT +-+','FR3-IMGT +--','FR3-IMGT -+-','FR3-IMGT --+','FR3-IMGT ---','FR3-IMGT Very similar',
            'FR3-IMGT Similar','FR3-IMGT Dissimilar','FR3-IMGT Very dissimilar','CDR3-IMGT Nb of positions including IMGT gaps (AA)',
            'CDR3-IMGT Nb of AA','CDR3-IMGT Nb of identical AA','CDR3-IMGT Nb of AA changes','CDR3-IMGT +++','CDR3-IMGT ++-', 'CDR3-IMGT +-+',
            'CDR3-IMGT +--','CDR3-IMGT -+-','CDR3-IMGT --+','CDR3-IMGT ---','CDR3-IMGT Very similar','CDR3-IMGT Similar','CDR3-IMGT Dissimilar',
            'CDR3-IMGT Very dissimilar'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION Nb of positions including IMGT gaps (AA)',
            'V-REGION Nb of AA','V-REGION Nb of identical AA','V-REGION Nb of AA changes','FR1-IMGT Nb of positions including IMGT gaps (AA)','FR1-IMGT Nb of AA',
            'FR1-IMGT Nb of identical AA','FR1-IMGT Nb of AA changes','CDR1-IMGT Nb of positions including IMGT gaps (AA)','CDR1-IMGT Nb of AA',
            'CDR1-IMGT Nb of identical AA','CDR1-IMGT Nb of AA changes','FR2-IMGT Nb of positions including IMGT gaps (AA)',
            'FR2-IMGT Nb of AA','FR2-IMGT Nb of identical AA','FR2-IMGT Nb of AA changes','CDR2-IMGT Nb of AA','CDR2-IMGT Nb of identical AA',
            'CDR2-IMGT Nb of AA changes','FR3-IMGT Nb of positions including IMGT gaps (AA)','FR3-IMGT Nb of AA','FR3-IMGT Nb of identical AA','FR3-IMGT Nb of AA changes',
            'CDR3-IMGT Nb of positions including IMGT gaps (AA)','CDR3-IMGT Nb of AA','CDR3-IMGT Nb of identical AA','CDR3-IMGT Nb of AA changes'
        ] 
    },
    {
        'Name':  '10_V-REGION-mutation-hotspots',
        'Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele','(a/t)a','t(a/t)','(a/g)g(c/t)(a/t)','(a/t)(a/g)c(c/t)'
        ],
        'Important_Headers':  [
            'Sequence number','Sequence ID','Functionality','V-GENE and allele'
        ]
    },
    {
        'Name':  '11_Parameters',
        'Headers':  [
            'Date','IMGT/V-QUEST programme version','IMGT/V-QUEST reference directory release','Species','Receptor type or locus','IMGT/V-QUEST reference directory set', 'Search for insertions and deletions', "Nb of nucleotides to add (or exclude) in 3' of the V-REGION for the evaluation of the alignment score", "Nb of nucleotides to exclude in 5' of the V-REGION for the evaluation of the nb of mutations", 'Number of submitted sequences'
        ],
        'Important_Headers':  [
            'Date','IMGT/V-QUEST programme version','IMGT/V-QUEST reference directory release','Species','Receptor type or locus','IMGT/V-QUEST reference directory set', 'Search for insertions and deletions', "Nb of nucleotides to add (or exclude) in 3' of the V-REGION for the evaluation of the alignment score", "Nb of nucleotides to exclude in 5' of the V-REGION for the evaluation of the nb of mutations", 'Number of submitted sequences'
        ],                            
        'Notes':  'PERFORM_TRANSPOSE_BEFORE_READING_DATA'
    }
]
