
//global variable that will store metadata for each experiment returned from query 
var exp_metadata_info = {}
var exp_metadata_list = []

var analysis_schema_datatypes = {}
var global_annotation = {}
var all_analyses = []
var all_recomb = []
var exp_user
var expdT
var global_table_header
var global_use_qtip
var current_selection = []
var current_exp_set = []

//var->map experimentid+recombination type: schema and experimentid+annotation: schema
var analysis_schema_annotations = {}

//simillar to above var, but  more info 
//->map experimentid+recombination type: schema and experimentid+annotation: schema
//replace field names containining '.' with '___'

var analysis_schema_annotations_html = {}
var analysis_schema_annotations_html_by_analysis = {}



//this variable will control how the 'OUTPUT FIELD ORDER'/'MODIFY FIELD OUTPUT' will appear
//the order of these fields will match the order they will be returned as outputs to TAB/CSV files (or FASTA/FASTQ)
//any fields NOT in this list will appear at the end of the list 
//the prefix DATA. is intentially left off from this list 
var default_field_sort_order = [
	'SEQ_ID',
	'SEQUENCE_HEADER',
	'SEQUENCE',
	'QUALITY_SCORE',
	'PREDICTED_AB_SEQ.AA',
	'PREDICTED_AB_SEQ.NT',
	'STRAND',	
	'PRODUCTIVE',
	'RECOMBINATION_TYPE',
	'PREDICTED_CHAIN_TYPE',
	'LOCUS_NAME',

	'VREGION.VGENES',
	'DREGION.DGENES',
	'JREGION.JGENES',	
	
	'FULL_LENGTH',
	'STOP_CODONS',
	
	'VREGION.FR1.AA',
	'VREGION.CDR1.AA',
	'VREGION.FR2.AA',
	'VREGION.CDR2.AA',
	'VREGION.FR3.AA',	
	'CDR3.AA',
	'JREGION.FR4.AA',
	
	'VREGION.FR1.NT',	
	'VREGION.CDR1.NT',
	'VREGION.FR2.NT',
	'VREGION.CDR2.NT',
	'VREGION.FR3.NT',	
	'CDR3.NT',
	'JREGION.FR4.NT',
	
	'VREGION.CDR1.AA_LENGTH',
	'VREGION.CDR1.NT_LENGTH',
	'VREGION.CDR2.AA_LENGTH',
	'VREGION.CDR2.NT_LENGTH',
	'CDR3.AA_LENGTH',
	'CDR3.NT_LENGTH',
	'VREGION.VGENE_SCORES',
	'JREGION.JGENE_SCORES',
	'DREGION.DGENE_SCORES',
	'VREGION.SHM.NT',
	'VREGION.SHM.AA',
	
	'GAPPED.VREGION.FR1.AA',
	'GAPPED.VREGION.CDR1.AA',
	'GAPPED.VREGION.FR2.AA',
	'GAPPED.VREGION.CDR2.AA',
	'GAPPED.VREGION.FR3.AA',	
	'GAPPED.CDR3.AA',
	'GAPPED.JREGION.FR4.AA',
	
	'GAPPED.VREGION.FR1.NT',
	'GAPPED.VREGION.CDR1.NT',
	'GAPPED.VREGION.FR2.NT',
	'GAPPED.VREGION.CDR2.NT',
	'GAPPED.VREGION.FR3.NT',	
	'GAPPED.CDR3.NT',
	'GAPPED.JREGION.FR4.NT',
	
	'VREGION.VGENE_QUERY_START',
	'VREGION.VGENE_QUERY_END',
	'JREGION.JGENE_QUERY_START',
	'JREGION.JGENE_QUERY_END'
	
]

//similar to varible above 
//EXCEPT..REGARDLES OF ANY NEW FIELDS FOUND IN THE SCHEMA, THESE FIELDS WILL APPEAR LAST/AT THE END AND IN THIS ORDER
var default_field_sort_order_always_appear_last = [
	'NOTES',
	'EXP_ID',
	'EXPERIMENT_NAME',
	'PROJECT_NAME',
	'ANALYSIS_NAME',
	'DATE_UPDATED',
	'FILENAME',
	'SETINGS'	
]

//FOR SORTING THE TREE DATA (CHECKBOX TREE)
//MAKING THE CHECKBOX TREE IS A RECURSIVE FUNCTION THAT PROGRESSIVELY REMOVES STRINGS PREPENDED WITH X.Y (SO FIRST WE SORT BY THE FIRST LAYER X, THEN WE DEFINE SORT ORDER FOR THE SECOND LAYER => Y)
var sort_tree_in_this_order = [
	'SEQ_ID',
	'EXP_ID',
	'ANALYSIS_NAME',
	'RECOMBINATION_TYPE',
	'EXPERIMENT_NAME',	
	'PROJECT_NAME',
	'DATE_UPDATED',
	'DATA',	
	'SEQUENCE_HEADER',
	'SEQUENCE',
	'QUALITY_SCORE',
	'STRAND',
	'PRODUCTIVE',
	'PREDICTED_AB_SEQ',
	'PREDICTED_AB_SEQ',		
	'PREDICTED_CHAIN_TYPE',
	'LOCUS_NAME',
	'FULL_LENGTH',
	'STOP_CODONS',
	'VREGION',
	'CDR3',
	'JREGION',
	'DREGION',				
	'GAPPED',
	'NOTES',
	'VGENES',
	'VGENE_SCORES',	
	'JGENES',
	'JGENE_SCORES',
	'DGENES',
	'DGENE_SCORES',
	'FR1',
	'CDR1',
	'FR2',
	'CDR2',
	'FR3',	
	'VGENE_QUERY_START',
	'VGENE_QUERY_END',
	'JGENE_QUERY_START',
	'JGENE_QUERY_END',
	'NT',
	'AA',
	'NT_LENGTH',
	'AA_LENGTH'
]
	


//sorts an 2D array by its first column 
function SortArray(row1,row2){	
	if (row1[0] === row2[0]) {
		return 0;
	}
	else {
		return (row1[0] < row2[0]) ? -1 : 1;
	}
}

//sorts an 2D array by its second column 
function SortArrayP2(row1,row2){	
	if (row1[1] === row2[1]) {
		return 0;
	}
	else {
		return (row1[1] < row2[1]) ? -1 : 1;
	}
}

//SORT THE POTENTIAL TREE FOR SELECTING OUTPUT FILES 
//THIS WILL SHOW THE ORDER WE WANT THE CHECKBOX TREE TO APPEAR 
function SortSchemaTree(row1,row2){
	//take the schema names from both elements (row1, row2)
	//loop through the var sort_tree_in_this_order
	//if either elemnet row1/row2 starts with a field defined in the array, 
	//then that is its sort order 
	//if it doesnt exist, then have it return by default the order of their fields  
	var r1pos = sort_tree_in_this_order.length+1
	var r2pos = sort_tree_in_this_order.length+1 		
		
	
	for (var i = 0; i<sort_tree_in_this_order.length;i++){			
		if (row1 == sort_tree_in_this_order[i] || row1.startsWith(sort_tree_in_this_order[i]+'.'))
			r1pos = i
		if (row2 == sort_tree_in_this_order[i] || row2.startsWith(sort_tree_in_this_order[i]+'.'))
			r2pos = i
	}
	
	//since they are both equally weighted, return the sort order allphabetically by their value
	if (r1pos==r2pos){
		if (row1 === row2) {
			return 0;
		}
		else {
			return (row1< row2) ? -1 : 1;
		}	
	}
	else{
		//return by their sort position 
		return (r1pos < r2pos) ? -1 : 1;
	}
	
	
}


function GetSortedFields(parent_div ){
	//find all the checkbox in this current level /nodes 		
	//find all the checkbox in this current level /nodes 	
	if (typeof parent_div == 'undefined') {
        parent_div = 'schema-fields';
    }
	
	
	var me = $('#'+parent_div)	
	var field_array = []
	var name
	if(me.children('div.check_box_parent').length==0){
		if(me.children('input:checkbox').is(':checked')){
			name = me.attr('name').replace(/___/g,'.')
			if (name.startsWith('DATA.'))
				name = name.substr(5)
			field_array.push(name)
		}
	}
	else{
		console.log('a')
		me.children('div.check_box_parent.foundfield').each(function(){ //only choose non hidden objects that are sub sections of this 				
			field_array = field_array.concat(GetSortedFields($(this).attr('id')))
			console.log($(this))
		})
	}		
	return field_array
	
}

//This will take ALL of the data from exp metadata and generate a dictionaty that contains all possible fields from all experiments
//var global_schema = {}
//input a list of experiment metadata
//From all experiment metadata, it will determine the overall schema based on current settings
function GetGlobalSchema(exp_data){
	var combined_schema = {}
	for (each_exp in exp_data){
		for(all_schema_fields in exp_data[each_exp])
			combined_schema[all_schema_fields] = 1
	}
	//THESE FIELDS DO NOT APPEAR IN THE 'SCHEMA' UNDER EXPERIMENT METADATA BUT MAY BE USEFUL TO ALLOW USER TO OUTPUT THEM 
	var fields_in_schema = ['SEQ_ID','EXP_ID','FILENAME','RECOMBINATION_TYPE','ANALYSIS_NAME','EXPERIMENT_NAME','PROJECT_NAME','DATE_UPDATED','SETTINGS','DATA.SEQUENCE_HEADER','DATA.SEQUENCE','DATA.QUALITY_SCORE']
	fields_in_schema = fields_in_schema.concat(Object.keys(combined_schema))
	
	
	//NOW LETS CONTROL THE ORDER OF THE FIELDS IN THE SCHEMA 
	/*var sorted_fields_in_schema = []
	var sorted_fields_last = []
	var missing_fields = []
	var index
	//loop through each of the fields found in the global schema 
	for (var i = 0; i<fields_in_schema.length;i++){
		//search for its position in the defaultsortorder arrays described above 
		//if a field is not found, then add it to missing_fields array 
		index = default_field_sort_order.indexOf(fields_in_schema[i])
		if(index>-1){ //found in the main sort_order_array 
			sorted_fields_in_schema.push([fields_in_schema[i],index])
		}
		else{
			index = default_field_sort_order_always_appear_last.indexOf(fields_in_schema[i])
			if(index>-1)//found in the sort order array for elements that always appear last 
				sorted_fields_last.push([fields_in_schema[i],index])
			else
				missing_fields.push(fields_in_schema[i])//just store this as a new field found 
		}
	}*/	
	//now sort the fields by their order defined in sort_order_array 
	//sorted_fields_in_schema.sort(SortArrayP2)
	
	//NOW sort the global schema fields by their prefereed tree order 
	//fields_in_schema.sort(SortSchemaTree)
	
	
	
	//based on the fields in schema, create an HTML checkbox tree using the scortschematree sort fucntion
	MakeCheckBoxFamily_dictionary_flattened('schema-fields',fields_in_schema,SortSchemaTree)
	
	
	//Now by default set the following filds as checked 
	var check_default_fields = ['SEQ_ID','EXP_ID','ANALYSIS_NAME','RECOMBINATION_TYPE','DATA']
	
	for (var  i = 0; i<check_default_fields.length;i++){		
		$("input[name='"+check_default_fields[i]+"']").trigger("click")
	}
	$("input[name='SEQ_ID']").prop('disabled',true)//always force SEQ_ID to be selected
	
	return fields_in_schema
}


//returns a unique array 
function MakeUnique(input_array,sort){
	if (typeof sort == 'undefined') {
        sort = false;
    }
	var uniqueNames = [];
	$.each(input_array, function(i, el){
		if($.inArray(el, uniqueNames) === -1) uniqueNames.push(el);		
	});
	if (sort)
		uniqueNames.sort()
	return uniqueNames
}


//main function for generating the app. 
//It will accept two fields: 
	//a list of dictionaries describing metadata for each experiment 
	//a string corresponding to the users name in the database (want to avoid use appsoma_api.get_user_info)
function InitializeAPP(info_from_db){
	if ('metadata' in info_from_db)
		metadata = info_from_db.metadata;
	else
		metadata = []
	if('user_data' in info_from_db)
		db_username = info_from_db.user_data.user
	else
		db_username = ''
	
	var IgnoreFields = ['DUPLICATED_FIELDS','ANALYSES_SETTINGS']
	var temp_dict 
	exp_user = db_username
	//parse metadata provided 
	//remove keys we dont want to show 
	//also set global variables for the app (schema fields, annotation methods) 
	for (var i = 0; i<metadata.length;i++){
		analysis_schema_annotations[metadata[i]['_id']] = {}
		analysis_schema_datatypes[metadata[i]['_id']] = {}
		analysis_schema_annotations_html[metadata[i]['_id']] = {}
		exp_metadata_info[metadata[i]['_id']] = {}
		for (each_key in metadata[i]){
			if (each_key == 'ANALYSIS_SCHEMA'){				
				temp_dict = FlattenDictionary(metadata[i][each_key])
				for (schema_fields in temp_dict){
					if(schema_fields.endsWith('.DATATYPE'))
						analysis_schema_datatypes[metadata[i]['_id']]['DATA.'+schema_fields.substr(0,schema_fields.length-'.DATATYPE'.length)] = temp_dict[schema_fields]					
					else if(schema_fields.endsWith('.ANALYSES')){
						//important! prepend fields with data!
						analysis_schema_annotations[metadata[i]['_id']]['DATA.'+schema_fields.substr(0,schema_fields.length-'.ANALYSES'.length)] = temp_dict[schema_fields]	
						analysis_schema_annotations_html[metadata[i]['_id']]['DATA___'+schema_fields.substr(0,schema_fields.length-'.ANALYSES'.length).replace(/\./g,'___')] = temp_dict[schema_fields]	
						for (var j = 0; j<temp_dict[schema_fields].length;j++){//now add an extra variable where we add in the analysis name (we use this for later filtering by annotation types)
							if(!(metadata[i]['_id']+'_'+temp_dict[schema_fields][j] in analysis_schema_annotations_html_by_analysis)){
								analysis_schema_annotations_html_by_analysis[metadata[i]['_id']+'_'+temp_dict[schema_fields][j]] = {}																
							}
							analysis_schema_annotations_html_by_analysis[metadata[i]['_id']+'_'+temp_dict[schema_fields][j]]['DATA___'+schema_fields.substr(0,schema_fields.length-'.ANALYSES'.length).replace(/\./g,'___')] = temp_dict[schema_fields]	
						}	
					}
				}
			}
			else if (each_key == 'ANALYSES_COUNT'){					
				exp_metadata_info[metadata[i]['_id']]['ANALYSES'] = Object.keys(metadata[i][each_key])
				all_analyses = all_analyses.concat(exp_metadata_info[metadata[i]['_id']]['ANALYSES']) //keep track of all observed analyses
				exp_metadata_info[metadata[i]['_id']]['ANALYSES_COUNT'] = metadata[i][each_key]		
				for (each_analysis in metadata[i]['ANALYSES_COUNT']){
					all_recomb=MakeUnique(all_recomb.concat(all_recomb,Object.keys(metadata[i]['ANALYSES_COUNT'][each_analysis])))
				}
			}
			else if(IgnoreFields.indexOf(each_key)==-1){
				exp_metadata_info[metadata[i]['_id']][each_key] = metadata[i][each_key]	
			}		
		}
		exp_metadata_list.push(exp_metadata_info[metadata[i]['_id']])
	}	
		
	
	//make all analyses unique and then sort 
	all_analyses = MakeUnique(all_analyses,true)
	all_recomb = MakeUnique(all_recomb,true)
	current_exp_set = exp_metadata_list
	expdT = MakeHTMLTable('exp-table-goes-here','name',[{'colname':'PROJECT_NAME'},{'colname':'EXPERIMENT_NAME'},{'colname':'READ_ACCESS'},{'colname':'OWNERS_OF_EXPERIMENT'},{'colname':'DESCRIPTION'},{'colname':'SEQ_COUNT'},{'colname':'ANALYSES'},{'colname':'_id','hidden':true}],exp_metadata_list,true,false,true,'Experiment Metadata',[],'','300px','200px', '',false)	
	
	//adds annotation section to the html page
	AddAnnotationFilter()
	AddAnalysesRecombResults()
	
	//figure out ALL possible annotation fields based on users experiments. //create checkboxes for selecting specific fields 
	global_annotation= GetGlobalSchema(analysis_schema_annotations)
	
	TableListener()
	FindAnnotatedFields([])		
}


//will add html code inot the page. Basically undername the section for 'selected annotated experiments' it will add fields for specific annotation types 
function AddAnnotationFilter(){
	//The global variable 'all_analyses' stores all annotation types found in the provided experiments (this could change by users running app)
	var subsection_html = ''
	
	for (var i = 0; i<all_analyses.length;i++){
		//loop through eachannotation type and add html below the tag for ffilterby annotation 
		subsection_html += "<div class='exp-filter div-btn' name='"+all_analyses[i]+"'>"+
								all_analyses[i]+
								"<input type='checkbox' class='hidden'>"+
							"</div>"
	}
	
	//now find the parent div and replace its child (subsection) wiht this html 	
	$('.filter-annotation').children('.subsection').html(subsection_html)
}


//will add html code inot the page. Basically undername the section for 'selected annotated experiments' it will add fields for specific annotation types and recombination types 
function AddAnalysesRecombResults(){
	//The global variable 'all_analyses' stores all annotation types found in the provided experiments (this could change by users running app)
	var subsection_html = ''
	
	for (var i = 0; i<all_analyses.length;i++){
		//loop through eachannotation type and add html below the tag for ffilterby annotation 
		subsection_html += "<div name='"+all_analyses[i]+"' class='report div-btn analyses_program annotation_"+all_analyses[i]+"'>"+
								"<span name='"+all_analyses[i]+"'>"+all_analyses[i]+"</span>"+
								"<input type='checkbox' class='hidden'>"+
							"</div>"
	}
	subsection_html+="<div class='space_up_down_small'>-----------</div>"
	for (var i = 0; i<all_recomb.length;i++){
		//loop through eachannotation type and add html below the tag for ffilterby annotation 
		subsection_html += "<div name='"+all_recomb[i]+"' class='report div-btn recombination recombination_"+all_recomb[i]+"'>"+
								"<span name='"+all_recomb[i]+"'>"+all_recomb[i]+"</span>"+
								"<input type='checkbox' class='hidden'>"+
								"</div>"
	}
	
	//now find the parent div and replace its child (subsection) wiht this html 	
	$('.annotateddata').children('.subsection').html(subsection_html)
}



//input -> nested dictionary
//output -> flattens contents such that all nested objects are replaced by '.' character
	//{'vregion':{'cdr1':{'nt':'a'}}} => becomes => vretion.cdr1.nt:'a'
function FlattenDictionary(nested_json,parent_string,flattened_json){

	if (typeof parent_string == 'undefined') {
        parent_string='';
    }
	
	
	if (typeof flattened_json== 'undefined') {
        flattened_json = {};
    }

	var me = nested_json	
	var prefix
	if (parent_string == '')
		prefix = ''
	else
		prefix = parent_string+'.'
		
	for (each_key in me){		
		if(me[each_key].constructor==Object){					
			FlattenDictionary(me[each_key],prefix+each_key,flattened_json)
		}
		else{
			flattened_json[prefix+each_key] = me[each_key]
			
		}
	}	
	return flattened_json	
}



//input location of where to make the interface
//a NESTED (not flat) json object 
//output => will make html checkbox items that are organized in a tree like structure
//clickling/selecting a parent will select all children below it 
function MakeCheckBoxFamily(divname,json,sort_fxn,parent_node){
	if (typeof sort_fxn == 'undefined') {
        sort_fxn = '';
    }
	var me = json	
	var child_node = 0
	html_string = ''
	var fieldnames=[]
	for (each_key in me){
		fieldnames.push(each_key)
	}
	
	if(sort_fxn=='')	
		fieldnames.sort()		
	else
		fieldnames.sort(sort_fxn)		
		
	for (var i = 0;i<fieldnames.length;i++){
		child_node+=1
		html_string+="<div class='space_up_down_small indent1 check_box_parent' id='Leaf"+parent_node+'_'+String(child_node)+"'><input type='checkbox' name='"+fieldnames[i]+"'><span class='checkbox-label'>"+fieldnames[i]+"</span></div>"
	}		
	$('#'+divname).append(html_string);	
	ListenCheckboxAction(divname);
	child_node=0;
	for (var i = 0;i<fieldnames.length;i++){
		child_node+=1;				
		if(me[fieldnames[i]].constructor==Object){		
			MakeCheckBoxFamily('Leaf'+parent_node+'_'+String(child_node),me[fieldnames[i]],sort_fxn,String(parseInt(parent_node)+1))
		}
	}	
}


//takes a LIST of KEY names from a FLATTENED DICTIONARY and using '.' to designated nested objects will generate a checkbox tree
function MakeCheckBoxFamily_dictionary_flattened(divname,key_list,sort_fxn,prefix){
	if (typeof sort_fxn == 'undefined') {
        sort_fxn = '';
    }
	
	if (typeof prefix == 'undefined') {
        prefix = '';
    }

	if (prefix!='')
		prefix = prefix+'___' //WE HAVE TO USE THIS CHARACTER BECAUSE '.' DOES NOT WOKR CORERCTLY MUST BE INCOMPITABLE WITH INTERNAL JAVASCRIPT/JQUERY 

	//var me = json		
	html_string = ''
	var fieldnames=[]
	var outerkeys = []
	var key_string=''
	var current_key=''
	
	var nested_keys = {}
	
	//sort all keys in dictionary 
	//for (each_key in me){
	//	fieldnames.push(each_key)		
	//}		
	if (sort_fxn=='')
		fieldnames = key_list.sort()		
	else
		fieldnames = key_list.sort(sort_fxn)		
			
	for (var i = 0;i<fieldnames.length;i++){
		key_string = fieldnames[i].split('.')
		if(key_string[0]!=current_key){
			current_key=key_string[0]
			outerkeys.push(current_key)						
			if (key_string.length>1)
				nested_keys[current_key] = [key_string.slice(1).join('.')]						
		}
		else{
			if (key_string.length>1)
				nested_keys[current_key].push(key_string.slice(1).join('.'))
		}				
	}			
	
	for (var i =0; i<outerkeys.length;i++)	
		html_string+="<div class='space_up_down_small indent1 check_box_parent' id='"+prefix+outerkeys[i]+"' name='"+prefix+outerkeys[i]+"'><input type='checkbox' name='"+prefix+outerkeys[i]+"'><span class='checkbox-label'>"+outerkeys[i]+"</span></div>"			
		
	$('#'+divname).append(html_string);	
	ListenCheckboxAction(divname);
	
	for (var key in nested_keys)		
		//divname of parent key = > prefix+key
		MakeCheckBoxFamily_dictionary_flattened(prefix+key, nested_keys[key],sort_fxn,prefix+key )
	
}





//Add javascript functionality to checkboxes 
//when a parent is selected, the select all of its childrent
//when a div is unselected, then deselect its parent 
function ListenCheckboxAction(divname){	
	$('#'+divname).on('click','input:checkbox',function(){		
		var check_value = $(this).is(':checked');
		
		//find all of the div within div name (siblings to the checkbox)
		var element = $(this).siblings('.check_box_parent')
		element.each(function(){
			//for each of these divs add the checkbox value 
			var subelement = $(this)
			subelement.find('input:checkbox').each(function(){				
				$(this).attr('checked',check_value)
			})
		})
		
		//now compare whether all of the sibling check boxes in same levle are checked
		
		//console.log($(this).parent().siblings('.check_box_parent:visible').children('input:checked').length+check_value) //include the check_value of current element 
		//console.log($(this).parent().siblings('.check_box_parent:visible').children('input').length+1) // include current element 
		
		var sibling_checkboxes = $(this).parent().siblings('.check_box_parent:visible')
		
		var bool_all_selected =sibling_checkboxes.children('input:checked').length+check_value == sibling_checkboxes.children('input').length+1
		
				
		if(bool_all_selected){
			//if they are all checke,d then the parent shoudl be checked
			$(this).parent().parent().children('input:checkbox').prop('checked',true)
		}
		else{
			//if not then the parent should be unchecked
			$(this).parent().parent().children('input:checkbox').prop('checked',false)
		}
	});
}

//this will loop through all of the chekcboxes under schema 
//based on the selected fields it will create a mongo projection query (basically a dictionary of values to project)
function GetMongoProjectionFromCheckboxes(parent_div,projection){
	//find all the checkbox in this current level /nodes 	
	if (typeof parent_div == 'undefined') {
        parent_div = 'schema-fields';
    }
	
	if (typeof projection == 'undefined') {
        projection = {}
    }
	
	var me = $('#'+parent_div)
	
	if(me.children('input:checkbox').is(':checked')){
		//add to the projection 
		projection[me.attr('name').replace(/___/g,'.')] = 1
		console.log('checked')
	}
	else{ //it is not selected, but may some of its leafs are selected 		
		me.children('div.check_box_parent').each(function(){ //only choose non hidden objects that are sub sections of this 							
			if($(this).is(':visible'))
				GetMongoProjectionFromCheckboxes($(this).attr('id') ,projection)
		})
	}
		
	return projection
}


//tool_tip_text_array => a list of text to associate as qtips to each of the rows 
function AddQTipTitle(datatable_var,tool_tip_text_array,tooltipheader){
	
	if (typeof tooltipheader == 'undefined') {
        tooltipheader = '';
    }
	
	//$('#'+tableid+' tbody tr').each(function(){{	--> DO NOT USE THIS, USE DATATABLE ROWS. TABLEID WILL ONLY LIST THE FIRST TEN ROWS
	datatable_var.$('tr').each(function(){
		
		
		$(this).attr('title',tool_tip_text_array[datatable_var.fnGetPosition(this)]);																		
	});	
	//finally associate a qtip jquery to all of the rows in the datatable
	datatable_var.$('tr').qtip({
		content:{	
			title: tooltipheader,
			text: false
		},
		show:{
			delay:1000
		},
		hide: {
			fixed: true, // Make it fixed so it can be hovered over
			delay: 100//add a delay from hiding tooltip
		},
		position:{
			target: 'mouse',
			adjust: {
			  mouse: false  
			}
		},
		style:{
		
		}
	});		
}

function RenderTable(datatable_var,data){	
	var dT = datatable_var;
	dT.fnClearTable();
	dT.fnAddData(data);		
	$(window).trigger( "resize" );//ADJUST TABLE HEIGHT AND COLUMNS SEE window resize function below 	
	return dT
}

function RefreshTableWidths(){
	$(window).trigger( "resize" );//ADJUST TABLE HEIGHT AND COLUMNS SEE window resize function below 	
}

//data => list of json docs. key corresponds to headr/colum name
//header=>list of dictioary containing three keys: 
	//colname
	//hidden: true/false
	//width: column width 	
/*,minColWidth='200px',maxColWidth='',wrap_text=false*/
//class_table_name => optional if user wants to add a new class to the table 
function MakeHTMLTable(divid,name,header,data,multiple_row_selection,use_cell_tooltips,use_row_qtips,qTipHeader,row_tooltips,class_table_name,maxTableHeight,minColWidth, maxColWidth,wrap_text){	
	
	
	
	if (typeof header == 'undefined') {
        header = [];
    }
	
	if (typeof data == 'undefined') {
        data = [];
    }
	
	if (typeof multiple_row_selection == 'undefined') {
        multiple_row_selection=true
    }
	
	if (typeof use_cell_tooltips == 'undefined') {
        use_cell_tooltips=false
    }
	
	if (typeof use_row_qtips == 'undefined') {
        use_row_qtips = false
    }
	
	if (typeof qTipHeader == 'undefined') {
        qTipHeader = ''
    }
	
	if (typeof row_tooltips == 'undefined') {
        row_tooltips=[]
    }
	
	if (typeof class_table_name == 'undefined') {
        class_table_name=''
    }
	
	if (typeof maxTableHeight == 'undefined') {
        maxTableHeight='500px'
    }
	
	if (typeof minColWidth == 'undefined') {
        minColWidth=''
    }
	if (typeof maxColWidth == 'undefined') {
        maxColWidth=''
    }
	
	if (typeof wrap_text == 'undefined') {
        wrap_text=false
    }
	
	
	var myDataTable;
	global_use_qtip=use_row_qtips
	//this variable will control visibilty to datatable
	var datatable_column_visibilty = []
	
	//initialize string for making header row 
	table_hd = "<tr>"
	//parse throguh header variable
		//based on keys, update datacolumn visibility field 
	for (var i = 0; i< header.length;i++){
		var headername = header[i]['colname'];		
		if( ('hidden' in header[i]) && header[i]['hidden']==true){
			
			if('width' in header[i])
				datatable_column_visibilty.push({'bVisible':false,'sWidth':String(header[i]['width'])+'px'})
			else
				datatable_column_visibilty.push({'bVisible':false})
		}
		else{
			if ('width' in header[i])
				datatable_column_visibilty.push({'sWidth':String(header[i]['width'])+'px'})
			else
				datatable_column_visibilty.push('null')					
		}
		table_hd+="<td class='"+headername+"'>"+headername+"</td>"		
	}
	table_hd+="</tr>"
	global_table_header = header
	
	//make html code for the necessary elements we will want 
	//we will inject html TABLE code and inject inputs/buttons 	
	html_code = "<table class='igrep-datatable "+class_table_name+"'>" +
	"<thead>"+table_hd+"</thead>"+
	"<tbody></tbody>"+
	"</table>"+
	"<div align='right'>"+
		"<input type='button' class='all-row-select multirowselect not-multiple' value='Select all filtered rows'>"+
		"<input type='button' class='all-row-deselect multirowselect not-multiple' value='Deselect all rows'>"+
		"<div class='table-errors' style='color:red;'></div>"
	"</div>"
		
	$('#'+divid).html(html_code);
	
	tableid = divid+' table.igrep-datatable'	
			
	//settings used for making datatable
	datatable_parameters = {
		'bSort':true,                                                                            		
		"sScrollX": "100%",
		"bScrollCollapse": true,
		"bJQueryUI": true,
		"sPaginationType": "full_numbers",
	//	"bAutoWidth": false, //does weird things to the scrollbar...it moves the horizontal scrollbar back after sorting..
		"autowidth":false,					
	//	"sScrollXInner": "100%",
	//	"bStateSave": true,
	//	"bDeferRender": true,							
		"bRetrieve": true,
		"aoColumns": datatable_column_visibilty,					 		
	}
	
	if (maxTableHeight!=''){
		if(typeof maxTableHeight!='string')
			maxTableHeight = String(maxTableHeight)+'px'				
		datatable_parameters["sScrollY"] = maxTableHeight
	}
			
	//add this key to the variable above if we want to use cell tool tips
	//basically each of the values in the table will appear as tooltips also 
	if (use_cell_tooltips==true){
		datatable_parameters['fnCreatedRow']=function(nRow,aData,iDataIndex){
					for (var i =0 ;i<aData.length;i++){
						$('td:eq('+i.toString()+')', nRow).attr('title',aData[i]);
					}				
				}
	}
	
	myDataTable = $('#'+tableid).dataTable(datatable_parameters)
				
		
	if (multiple_row_selection){
		$('.multirowselect').removeClass('not-multiple')
	
		//Allow the HTML table to be selectable using the jquery selectable	
		$('#'+tableid+' tbody').selectable({			
			items:'tr',
			stop: function( event, ui ) {					
				var current_rows
				if(event.ctrlKey){//the user is pressing down the ctrl key for selection (so perhaps we want to save all hihglihged rows in all pages 
					//find all selected rows that are ONLY in the filtered table 
					current_rows = myDataTable.$('tr.ui-selected', {"filter": "applied"}); 
				}
				else{
					//only find selected rows that IN THE CURRENT PAGE!
					current_rows = myDataTable.$('tr.ui-selected', {"page":"current"});
				}				
				//remove all rows in the entire datatable (including rows not filtered)
				myDataTable.$('tr').removeClass('ui-selected');
				//now re-select the rows from first statement (selected rows only in filtered table) 				
				current_rows.addClass('ui-selected');		
				//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
				myDataTable.$('tr').first().trigger('click')
				
			}
		}).disableSelection()
		
		//add javascript/jquery to allow multiple selection of fields
		$('#'+divid+' .all-row-deselect').on('click',function(){					
			//FIND ALL ROWS in the datatable var (must use datatable() function to include non -filtered rows
			myDataTable.$('tr').removeClass('ui-selected');	
			//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
			myDataTable.$('tr').first().trigger('click')
		});				
		
		$('#'+divid+' .all-row-select').on('click',function(){
			//FIND ROWS THAT ARE FILTERED
			var filteredrows = myDataTable.$('tr', {"filter": "applied"});					
			myDataTable.$('tr').removeClass('ui-selected'); //remove class from ALL rows in table
			filteredrows.addClass('ui-selected');//add class to filtered rows only 
			//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
			myDataTable.$('tr').first().trigger('click')
		});		
	}
	else{ //ALLOW only one row to be selected
				
		$('#'+tableid).on('click', 'tbody tr', function() { // tbody tr').click( function( event ) {												
			if ( $(this).hasClass('ui-selected') ) {						
				$(this).removeClass('ui-selected');						
			}                                              
			else {                             
				myDataTable.$('tr').removeClass('ui-selected'); //remove class from ALL rows in table
				$(this).addClass('ui-selected');  												
			}; 								 					
		});
	}

	var timer_id;	
	$(window).resize(function() {				
		console.log('resizing!!')
		clearTimeout(timer_id);
		timer_id = setTimeout(function() {    		        						
			myDataTable.fnAdjustColumnSizing();
		}, 300);
	});
					
	
	//NOW populate the actual table wiht the data variable 
	myDataTable = ModifyHTMLData(myDataTable,header,data,use_qtip=use_row_qtips,qtipname=qTipHeader)		
	
	//Use jquery to change some css styles provided by user 
	var css_styles = {}
	if (minColWidth!='')
		css_styles['min-width'] = String(minColWidth).replace('px','')+'px' //use .replace just incease user passed in a number. so, i.e.,this style will work if user said 200px OR 200 
	if (maxColWidth!='')
		css_styles['max-width'] = String(maxColWidth).replace('px','')+'px' //use .replace just incease user passed in a number. so, i.e.,this style will work if user said 200px OR 200 
	
		//modfiy table header css
	//myDataTable.$('thead td').css(css_styles)
	$('#'+divid+' table.dataTable thead td').css(css_styles)
	
	
	
	if (wrap_text){
		css_styles['white-space'] = 'normal' 
		css_styles['word-wrap'] = 'break-word'		
	}
	else{
		css_styles['white-space'] = 'nowrap' 
		css_styles['text-overflow'] = 'ellipsis'			
	}		
		//modify table body css 
		//use datatable var if, use just table then only modifies first 10 rows 
	myDataTable.$('td').css(css_styles)
				
	return myDataTable;
}


//data =>list, each element of list is a dict where key = colname
//header => list, each element is list. key in each element called colname refers to column name 		
function ModifyHTMLData(myDataTable,header,data,use_qtip,qtipname,row_tooltips){		

	if (typeof use_qtip == 'undefined') {
        use_qtip = false
    }
	
	if (typeof qtipname == 'undefined') {
		qtipname=''
    }
	
	if (typeof row_tooltips == 'undefined') {
        row_tooltips=[]
    }

	//NOW populate the actual table wiht the data variable 
	var rows_of_data = []
	//data =>list, each element of list is a dict where key = colname
	//header => list, each element is list. key in each element called colname refers to column name 	
	for (row_num= 0; row_num<data.length;row_num++){
		data_row = data[row_num]
		row_info = []
		for (i = 0; i<header.length;i++){
			var colname = header[i]['colname']
			if (colname in data_row)
				if(data_row[colname].constructor==Array){					
					row_info.push(data_row[colname].join(','));
				}
				else
					row_info.push(String(data_row[colname]));
				
			else
				row_info.push('')
		}
		rows_of_data.push(row_info );
	}	
	myDataTable = RenderTable(myDataTable,rows_of_data);
	

	//the user wants to use qtips for hovering and showing information at each row 
	if (use_qtip){
		//The user passed DID not pass in list of strings to use as q tips
		//So we will just use json dumps on the current data and display those results as q tips 
		if(row_tooltips=[]){
			for(var row_num = 0; row_num < data.length; row_num+=1){
				//use json dumps and an indentation of 4 spaces to present the data in a prety format 
				//use <pre> to make sure html reads in \t \n and spaces in toolstips
				row_tooltips.push('<pre>'+JSON.stringify(data[row_num],null,4)+'</pre>')
			}
		}
		//Now generate the qtips on the provided datatable 
		AddQTipTitle(myDataTable,row_tooltips,qtipname);
	}
	

	return myDataTable	
	
}

//THIS SORT PART MUST COME AFTER THE ABOVE FUNCTION TO ENSURE THAT ALL ROWS CAN BE SELLECTED....IF IT IS ABOVE THEN ONLY THE FIRST TEN ROWS ARE SELECTD
//pass in the 'div' parent that contains the datatable 	
function getDataTableData( divid ) { 					
	//DO THE FOLLOWING SEARCH: 
		//search for all table elements under divid whose class is datTable, 
		//then for each table only select tables that have a row with calss ui-selected
		//THEN go back a step and get the parent (which is the original table)
		//ONCE we do that then we can retrieve the original datatable 
	var tablefound=$('#'+divid+' table.dataTable tr.ui-selected').closest('table')
		
	if (tablefound.length==0){
		//no tables found 
		replyJson({'datable_results': {'divname':divid,'row_data':[]} } );																		
	}
	else{	
		//retrievent the datatable associated with this div
		dT = tablefound.dataTable()
		
		//get selected rows that ONLY Mmatch the dattable filter 
		var anselected=dT.$('tr.ui-selected', {"filter": "applied"}); // only return selected rows that are filtered		
		var rowData=[];
		//report results here (go through all selected rows and report the results from the table)
		for (i = 0; i<anselected.length;i++){					
			var iRow=dT.fnGetPosition(anselected[i]);     					
			var aData=dT.fnGetData(iRow);   		
			rowData.push([iRow,aData]);
		}
		replyJson({'datatable_results': {'divdname':divid,'row_data':rowData} } );																		
	}	
} 


//When useres sets exp filters using div-btn below, then this function is called 
//for all set filters, this funciton will filter experiments from databae 
function ResetTableFilters(){			
	current_exp_set = []
	var passed_filter
	var filter_type
	for (var i = 0; i<exp_metadata_list.length;i++){ //loop through all metadata
		
		exp_row = exp_metadata_list[i]
		passed_filter = true
		
		$('.exp-filter').each(function(){ //find al lhtml elements taht are filters 			
			
			if($(this).children('input').is(':checked')){ //only consider elements whose filters are chekced				
				filter_type	 = $(this).closest('.filter-parent')
				if(filter_type.hasClass('filter-my-exp')){ //filter by  experiments owned by uesr 
					if (exp_row['OWNERS_OF_EXPERIMENT'].indexOf(exp_user)==-1){
						passed_filter = false //user was not found in this experiment
					}			
				}
				else if(filter_type.hasClass('filter-seq-count')){
					if( !('SEQ_COUNT' in exp_row)||exp_row['SEQ_COUNT'] == 0)
						passed_filter = false
				}
				else if(filter_type.hasClass('filter-paired')){
					if ( !('VH:VL_PAIRED' in exp_row)||exp_row['VH:VL_PAIRED'] == false)
						passed_filter = false
				}
				else if(filter_type.hasClass('filter-annotation')){
					if (!('ANALYSES' in exp_row)) //experiment doesn't have this field which means it hassnt beenannotated
						passed_filter = false 
					else if ($(this).attr('name') != 'any_annotation'){ //if the filter attribute is 'any_annotation' then the previous if statement (analysis in exp_row) satisfies requriment
						if (exp_row['ANALYSES'].indexOf($(this).attr('name'))==-1) //this curernt anotation type was not found in the list in exp_row
							passed_filter = false
					}
					
				}
												
			}
		})
		if (passed_filter) //all filters passed 
			current_exp_set.push(exp_row)
	}
		
	expdT = ModifyHTMLData(expdT,global_table_header,current_exp_set,global_use_qtip,'Experiment Metadata')	
	TableListener()
	FindAnnotatedFields([])
}


//Based on the selected/highlighed experiments, this will got hrough and get all of the unique analyses and recombinations in teh selection
//then based on the selections it will go through the hidden and dsiabled data in the 'results' div-btns. 
//if there are sequences, then .ngsdata becomes enabled; if there are analyses then .annotateddata becomes enabled 
//all analyses and reomcbination typse found will then be made 'unhidden' so user can further select them as filter 
function FindAnnotatedFields(list_of_ids){
	var has_sequences = false
	var annotation_found = []
	var recomb_found = []
	var data
	for(var i = 0 ; i< list_of_ids.length; i++){
		data = exp_metadata_info[list_of_ids[i]]
		if('SEQ_COUNT' in data && data['SEQ_COUNT']>0){
			has_sequences = true			
			if('ANALYSES_COUNT' in data){ //data['ANALYSES_COUNT'] is a dictionary 
				for(each_analysis in data['ANALYSES_COUNT'] ){
					annotation_found.push(each_analysis)
					recomb_found=recomb_found.concat(recomb_found,Object.keys(data['ANALYSES_COUNT'][each_analysis]))
				}
			}			
		}					
	}	
	
		
	if (has_sequences){
		$('.ngsdata').removeClass('disabled')
		$('.alldata').removeClass('disabled')
			//THERE IS SEQUENCE DATA AVAILABLE FOR DOWNLOAD, 
			//BUT CURRENTLY THERE ARE NOT SELECTED OUTPUTS 
			//SO BY DEFAULT SELECT - > REPORT RAW ANNOTATION DATA 
		if($('.report.checked').length==0){
			$('.ngsdata').addClass('checked')
			$('.ngsdata').children('input').prop('checked',true)
			$('#run-seq-dl').prop('disabled',false)
			$('#modify-output-fields').prop('disabled',false)
		}
	}
	else{
		$('.ngsdata').addClass('disabled')
		$('.alldata').addClass('disabled')		
		$('.ngsdata').removeClass('checked')
		$('.alldata').removeClass('checked')	
		$('.ngsdata').children('input').prop('checked',false)
		$('.alldata').children('input').prop('checked',false)
		$('#modify-output-fields').prop('disabled',true)
		$('#run-seq-dl').prop('disabled',true)
		/*$('.annotateddata').removeClass('checked')
		$('.annotateddata').children('input').prop('checked',false)
		$('.annotateddata .subsection .div-btn').removeClass('checked')
		$('.annotateddata .subsection .div-btn').children('input').prop('checked',false) --> probably dont neeed this */
	}
	
	$('.annotateddata .subsection .report').hide()	
	
	if (annotation_found.length>0){
		annotation_found = MakeUnique(annotation_found,true)
		recomb_found = MakeUnique(recomb_found,true)
		$('.annotateddata').removeClass('disabled')		
		
		for (var i = 0; i<annotation_found.length;i++)
			$('.annotateddata .subsection .report.annotation_'+annotation_found[i]).show()
		
		for (var i = 0; i<recomb_found.length;i++)
			$('.annotateddata .subsection .report.recombination_'+recomb_found[i]).show()
		
	}
	else{
		$('.annotateddata').addClass('disabled')
		$('.annotateddata').removeClass('checked')
		$('.annotateddata').children('.subsection').hide()
		$('.annotateddata').children('input').prop('checked',false)
	}	
}

function TableListener(){
	
	//detect any rows from the HTML table that will be selected (remember we used selectable for the data so ui-selectee works)
	expdT.$('tr').click(function(){
		current_selection = []				
		
		//get selected rows that ONLY Mmatch the dattable filter 
		var anselected=expdT.$('tr.ui-selected', {"filter": "applied"}); // only return selected rows that are filtered						
		
		//report results here (go through all selected rows and report the results from the table)
		for (i = 0; i<anselected.length;i++){					
			var iRow=expdT.fnGetPosition(anselected[i]);   
			//get ids from the selected table 			
			current_selection.push(current_exp_set[iRow]['_id'])
		}
								
		FindAnnotatedFields(current_selection)
		FindSchemaAcrossSelection(current_selection)
	})
		
}


//based on the experiment selection , go through each experiments schemma and make them unhidden 
function FindSchemaAcrossSelection(current_exp_list){	

	//the names in the html elements have replace '.' with '___' '.' did not work for making html elements
	var ngs_fields = ['DATA___SEQUENCE','DATA___SEQUENCE_HEADER','DATA___QUALITY_SCORE','EXP_ID','SEQ_ID','DATE_UPDATED','FILENAME']
	var annotated_data_fields = ['ANALYSIS_NAME','RECOMBINATION_TYPE','SETTINGS','SEQ_ID','EXP_ID','DATE_UPDATED']
	
	//first make all schema fields hidden 	
	$('#schema').find('.check_box_parent').hide()
	
	//if ngs data is selected/raw sequence data is checked, then make ngs_Fields visible  
	if($('.ngsdata').hasClass('checked')){
		for (var i = 0; i<ngs_fields.length;i++){
			$(".check_box_parent[name='"+ngs_fields[i]+"']").show()	//show each field 
			$(".check_box_parent[name='"+ngs_fields[i]+"']").parentsUntil('#schema-fields').show()	//ALSO SHOW ALL OF ITS ANCESTORS UNTIL the div schema-fields		
		}
	}
	
	if ($('.annotateddata.checked').length==0) //no annotation is desired
		return 
	
	//NEXT DETERMINE IF USER HAS SELECTED ANY FORM OF ANNOTATION PROGRAM to be returned TO BE RETURNED 
	var annotated_found = $('.annotateddata').children('.subsection').find('.div-btn.checked.analyses_program')
	var check_these_analyses = []
	var show_from_all_annotation = false //show results from all anotation programs in experiment 
	if (annotated_found.length==0){ //no annotation fields were found, so instead check whether any recombination fields were checked
		annotated_found = $('.annotateddata .subsection').find('.div-btn.checked.recombination')
		if (annotated_found.length>0)//a recombination type was selected, so the user does what annotation. by default choose all annotation programs
			show_from_all_annotation = true
		else
			return //just leave function 
	}
	else{	
		//add specific annotation types to check 
		annotated_found.each(function(){
			check_these_analyses.push($(this).attr('name'))
		})
	}	

	for(var i = 0; i<annotated_data_fields.length;i++){
		$(".check_box_parent[name='"+annotated_data_fields[i]+"']").show()	//show each field 
	}
	
	if(show_from_all_annotation){
		//now go through the current list of experiments and for each epxeriment make its individual schema field visible 
		for (var i = 0; i<current_exp_list.length;i++){
			console.log(analysis_schema_annotations_html[current_exp_list[i]])
			for(each_field in analysis_schema_annotations_html[current_exp_list[i]]){//use the variable whose key is ONLY the experiment 								
				//also us ethe _html varaible because we need to replace '.' with '___'
				$(".check_box_parent[name='"+each_field+"']").show()		
				$(".check_box_parent[name='"+each_field+"']").parentsUntil('#schema-fields').show()//ALSO SHOW ALL OF ITS ANCESTORS UNTIL the div schema-fields
			}
		}	
	}
	else{
		
		for (var i = 0; i<current_exp_list.length;i++){		
			console.log(analysis_schema_annotations_html[current_exp_list[i]])
			
			for (var j = 0; j<check_these_analyses.length;j++){
				console.log(analysis_schema_annotations_html_by_analysis[current_exp_list[i]+'_'+check_these_analyses[j]])
				for(each_field in analysis_schema_annotations_html_by_analysis[current_exp_list[i]+'_'+check_these_analyses[j]]){ //instead use the variable whose key is ID AND ANALYSIS NAME 
					//also us ethe _html varaible because we need to replace '.' with '___'
					$(".check_box_parent[name='"+each_field+"']").show()		
					$(".check_box_parent[name='"+each_field+"']").parentsUntil('#schema-fields').show()//ALSO SHOW ALL OF ITS ANCESTORS UNTIL the div schema-fields
				}
			}
		}
			
	}
		
}




$(document).ready(function(){	
	//MakeHTMLTable('table-goes-here','name',[{'colname':'costas'},{'colname':'house'}],[{'costas':'nv1','house':'mate'},{'costas':'nv3','house':'namate'}])		
	//RefreshTableWidths()-->doesnt do anything rigth now 
	
	//anytime you click an div-btn element from the class 'left' then do the following 
	$('.left').on('click','.div-btn',function(){		
		var child=$(this).children('input')
		if ($(this).hasClass('disabled')){ //button is disabled so turn off all settings 
			$(this).removeClass('checked');
			child.prop('checked',false)
			$(this).siblings('.subsection').hide()
			$(this).siblings('.subsection').children().each(function(){
				$(this).removeClass('checked')
				$(this).children('input').prop('checked',false)
			})
		}
		else{				
		
			if($(this).hasClass('alldata')){//its a specific button that we want to treat caefully
					//lets add extra functionality to the 'REPORT' div-btns

				if($(this).children('input').is(':checked')){//we need to UNCHECK it 
					$(this).removeClass('checked')
					$(this).children('input').prop('checked',false)
				}
				else{
					//everything gets highlighted 
					$('.report.div-btn').each(function(){
						$(this).children('input').prop('checked',true)
						$(this).siblings('.subsection').show()
						$(this).addClass('checked')
					})
				}
			}
			else{
			
				child.prop("checked", !child.prop("checked"));				
				if (child.is(':checked')){
					$(this).addClass('checked')
					$(this).siblings('.subsection').show()			
				}
				else{
					$(this).removeClass('checked')
					$(this).siblings('.subsection').hide()
					$(this).siblings('.subsection').children().each(function(){
						$(this).removeClass('checked')
						$(this).children('input').prop('checked',false)
					})
				}
				if($(this).hasClass("exp-filter"))
					ResetTableFilters()
			}
		}	
	})
	
	
	//lets add extra functionality to the 'REPORT' div-btns that are not ALLDATA
	$('.left').on('click','.div-btn.report',function(){				
		if(!($(this).hasClass('disabled'))){					
			if (!($(this).hasClass('alldata')) ){									
				//onlyreporting a few fields, so lets turn off the 'ALLDATA' selected 				
				$('.report.div-btn.alldata').removeClass('checked')
				$('.report.div-btn.alldata').children('input').prop('checked',false)			
			}
		}
		
		var disable_buttons = $('.report.checked').length==0
		
		if(disable_buttons){
			$('#run-seq-dl').prop('disabled',true)
			$('#modify-output-fields').prop('disabled',true)
		}
		else{
			$('#run-seq-dl').prop('disabled',false)
			$('#modify-output-fields').prop('disabled',false)
		}
					
	})
	
	//events for changing the max number sequences per file 
	$('#query-seqs-per-file').on('change',function(){{
		if($(this).attr('value')=='custom'){{
			$('.query-custom-seqs').show();
		}}
		else{{
			$('.query-custom-seqs').hide();
		}}		
	}});
	
	//events for choosing to compress data 
	$('#query-compress').on('click',function(){{	
		if ($(this).is(':checked'))		
			$('#query-com-div').show();
		else
			$('#query-com-div').hide();
		
	}});
	
	//events for choosing to compress data 
	$('#set-limit').on('click',function(){{	
		if ($(this).is(':checked'))		
			$('#seq-limit').show();
		else
			$('#seq-limit').hide();
		
	}});
	
	
	$('#modify-output-fields').on('click',function(){		
		expdT.$('tr').first().trigger('click')
		$('#query-option-pane').hide("slide",{direction:"left"},500)//"slide", {direction: "left" }, 1000);
		$('#schema').show("slide",{ direction: "right" }, 500);
	});
	
	$('.back-btn').on('click',function(){		
		$('#schema').hide("slide",{direction:"right"},500)//"slide", {direction: "left" }, 1000);
		$('#query-option-pane').show("slide",{ direction: "left" }, 500);
	});
	
	//for saving to custom files, show the input box to type a file name 
	$("input[name='save-exp-by']").on('click',function(){

		if($(this).attr('value')=='custom')
			$(this).siblings('span').show()
		else
			$(this).siblings('span').hide()
	});
	
	//for hitting the submit button 
	$('#run-seq-dl').click(function(){
		//the user said they wanted annotationd ata, but there are no details provided (i.e. anlaysis name recombination type)
		var annotated = []
		var recombtype = []
		
		var experiments_to_dl = []
		
		if($('.annotateddata').hasClass('checked')){		
			var annotation_details = $('.annotateddata .subsection .div-btn.checked')
			if(annotation_details.length==0)
				$("span[name='error-dl']").show()
			else{
				$("span[name='error-dl']").hide()
				annotation_details.each(function(){
					if($(this).hasClass('analyses_program'))
						annotated.push($(this).attr('name'))
					else if ($(this).hasClass('recombination'))
						recombtype.push($(this).attr('name'))
				})
			}
		}			
		
		results = {}
		//get all settings for output fil e
		$("[name='output-data']").find('input,select').each(function(){
			if($(this).attr('type') == 'checkbox')
				results[$(this).attr('name')] = $(this).is(':checked')
			else
				results[$(this).attr('name')] = $(this).attr('value')
		});
		
		results['save-exp-by']=$("input[name='save-exp-by']").attr('value')
		
		results['ngs'] = $('.ngsdata').hasClass('checked')
		results['recombination_types'] = recombtype 
		results['analyses_name'] = annotated
		
		//find all selected experiments 
		current_selection = []						
		//get selected rows that ONLY Mmatch the dattable filter 
		var anselected=expdT.$('tr.ui-selected', {"filter": "applied"}); // only return selected rows that are filtered								
		//report results here (go through all selected rows and report the results from the table)
		for (i = 0; i<anselected.length;i++){					
			var iRow=expdT.fnGetPosition(anselected[i]);   
			//get ids from the selected table 			
			current_selection.push(current_exp_set[iRow]['_id'])
		}		
		results['_id'] = current_selection
		results['projection'] = GetMongoProjectionFromCheckboxes()
		console.log(results['projection'])
		//replyJson(results)		 	
	});
	
});

