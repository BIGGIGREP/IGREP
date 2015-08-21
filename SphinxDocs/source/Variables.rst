Global Variables
================
.. _Variables:

This covers the global variables used in all of the other files in Immunogrep.

List of variables
-----------------

    ====================   =====================   ==================

      **Variables**            **Value**             **Description**
    --------------------   ---------------------   ------------------
    descriptor_symbol      #@                      To signify the line has information needed to transact in immunogrep pipeline
    fasta_file_delimiter    <                      To signify additional info in sequence header
    idIdentifier            SEQ_ID                 For use with mongoDB
    expIdentifier           EXP_ID
    seqRawData              @SEQ
    description_var         DESCRIPTION
    translation_var         TRANSLATION
    tracktion_var           TRACKTION
    type_var                TYPE                   This is the analysis type key (i.e. 'TYPE':'IMGT' or 'TYPE':'INHOUSE', etc.)
    annotationExt           annotation
    analysisExt             analysis
    queryExt                Query
    IFExt                   IFFile
    listofextension                                 Consists of various extensions
    ====================   =====================   ==================

    .. note:: An example of the decorator label

        #@D{'DESCRIPTION':[list of label headers; order matters here],'TRANSLATON':{dictionary for communicating with the database}}

imgt_files_info
'''''''''''''''
    .. todo:: This part is ugly and I haven't figured out a better way of looking at it. Also not all 'Name's are in.

    this summarizes each of the column names associated with each imgt file from IMGT V-QUEST
        :Name: **1_Summary**
        :Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION score','V-REGION identity %',
            'V-REGION identity nt','V-REGION identity % (with ins/del events)','V-REGION identity nt (with ins/del events)',
            'J-GENE and allele','J-REGION score','J-REGION identity %','J-REGION identity nt','D-GENE and allele',
            'D-REGION reading frame','CDR1-IMGT length','CDR2-IMGT length','CDR3-IMGT length','CDR-IMGT lengths','FR-IMGT lengths',
            'AA JUNCTION','JUNCTION frame','Orientation','Functionality comment','V-REGION potential ins/del','J-GENE and allele comment',
            'V-REGION insertions','V-REGION deletions','Sequence'
        :Important_Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','V-REGION score','V-REGION identity %',
            'V-REGION identity nt','V-REGION identity % (with ins/del events)','V-REGION identity nt (with ins/del events)',
            'J-GENE and allele','J-REGION score','J-REGION identity %','J-REGION identity nt','D-GENE and allele',
            'D-REGION reading frame','CDR1-IMGT length','CDR2-IMGT length','CDR3-IMGT length','CDR-IMGT lengths','FR-IMGT lengths',
            'AA JUNCTION','JUNCTION frame','Orientation','Functionality comment','V-REGION potential ins/del','J-GENE and allele comment',
            'V-REGION insertions','V-REGION deletions','Sequence'

|
|

        :Name: **2_IMGT-gapped-nt-sequences**
        :Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        :Important_Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'

|
|

        :Name: **3_Nt-sequences**
        :Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele','V-D-J-REGION',
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
        :Important_Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele','V-D-J-REGION',
            'V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT','JUNCTION', 'FR4-IMGT','V-D-J-REGION start',
            'V-D-J-REGION end','V-J-REGION start','V-J-REGION end','V-REGION start','V-REGION end','FR1-IMGT start','FR1-IMGT end','CDR1-IMGT start',
            'CDR1-IMGT end','FR2-IMGT start','FR2-IMGT end','CDR2-IMGT start','CDR2-IMGT end','FR3-IMGT start','FR3-IMGT end','CDR3-IMGT start',
            'CDR3-IMGT end','JUNCTION start', 'JUNCTION end','D-J-REGION start','D-J-REGION end','J-REGION start','J-REGION end','FR4-IMGT start','FR4-IMGT end'

|
|

        :Name: **4_IMGT-gapped-AA-sequences**
        :Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'
        :Important_Headers: 'Sequence number','Sequence ID','Functionality','V-GENE and allele','J-GENE and allele','D-GENE and allele',
            'V-D-J-REGION','V-J-REGION','V-REGION','FR1-IMGT','CDR1-IMGT','FR2-IMGT','CDR2-IMGT','FR3-IMGT','CDR3-IMGT',
            'JUNCTION','J-REGION','FR4-IMGT'

|
|

        :Name:
        :Headers:
        :Important_Headers:

        :Name:
        :Headers:
        :Important_Headers:
