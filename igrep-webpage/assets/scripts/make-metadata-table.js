/*This is a wrapper script for displaying an experiment table in our UI
/*it will connect to the gglab database

	First attempt at widget making...
	It will be dependent on:
		1) jquery (version 1.11.4)
		2) jquery ui 1.11.4
		3) datatable (version 1.10.9)
		4) qtip (version 2.2.1)
		5) blockui.js(we have a local  version)

	In this first attempt, I was not able to get the 'load javascript functions'
	
	IMPORTANT: BECAUSE I COULD NOT GET THE LOAD-JAVASCRIPT FUCNTIONS TO WORK, I ASSUME THAT THE HTML PAGES ALREADY IMPORT THE REQUIRED DEPENDENCIES!!!

	BUT we will try to load the following scripts in this widget:
		boostrap.min.js and css (local)
		dashboard.css (local)		
		make-igrep-datatable.js (local)
		immunogrep_stylesheet.css (Local)
	
*/
	
(function() {

	var bootstrap_src = "/styles/bootstrap/js/bootstrap.min.js"
	var bootstrap_css = "/styles/bootstrap/css/bootstrap.min.css"
	var dashboard_css = "/styles/dashboard.css"
	//Load javascript functions
	function loadScript(src,on_complete){
		var script_tag = document.createElement('script');
		script_tag.type = "text/javascript";
		if (on_complete!=undefined){
			if (script_tag.readyState) {
				script_tag.onreadystatechange = function () { // For old versions of IE
					if (this.readyState == 'complete' || this.readyState == 'loaded') {
						script.onreadystatechange = null;
						on_complete();
					}
				};
			} else {			
				script_tag.onload = on_complete();
			}
		}
		script_tag.src = src
		// Try to find the head, otherwise default to the documentElement
		document.getElementsByTagName("head")[0].appendChild(script_tag);
	}

	/*LOAD BOOTSTRAP AND BLOCKUI JAVASCRIPT FUNCTIONS*/
	function load_all_scripts(){
		loadScript("/scripts/block_ui.js",loadigrep)		
	}

	function loadigrep(){
		loadScript("/scripts/igrep.script.tools.js")			
		loadScript("/scripts/make-igrep-datatable.js",loadbootstrap)				
	}

	function loadbootstrap(){
		loadScript(bootstrap_src,loadcss)
	}

	function loadcss(){
		var css_link = $("<link>",{ 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: bootstrap_css
		})
		css_link.appendTo('head');

		css_link = $("<link>",{ 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: dashboard_css
		})
		//css_link.appendTo('head');

		css_link = $("<link>",{ 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: "/styles/immunogrep_stylesheet.css"
		})
		css_link.appendTo('head');

		css_link = $("<link>",{ 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: "/styles/metadata.css"
		})
		css_link.appendTo('head');

		createMetadataUI()
	}

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
			if (row1 === sort_tree_in_this_order[i] || row1.startsWith(sort_tree_in_this_order[i]+'.'))
				r1pos = i
			if (row2 === sort_tree_in_this_order[i] || row2.startsWith(sort_tree_in_this_order[i]+'.'))
				r2pos = i
		}		
		//since they are both equally weighted, return the sort order allphabetically by their value
		if (r1pos===r2pos){
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

	//returns a unique array (order will not match starting array...., but can sort afterwards)
	function MakeUnique(input_array,sort){
		if (typeof sort == 'undefined') {
	        sort = false;
	    }
		var uniqueNames = [];

		var temp = {};
    	for (var i = 0; i < input_array.length; i++)
        	temp[input_array[i]] = true;
    	var uniqueNames = [];
    	for (var k in temp)
        	uniqueNames.push(k);    	
		
		if (sort)
			uniqueNames.sort()
		return uniqueNames
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
		if (parent_string === '')
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

	//global variable that will store metadata for each experiment returned from query 
	var exp_metadata_info = {}
	var exp_metadata_list = []
	var query_exp_data = {} //WE USE THIS VARIABLE TO DETERMINE WHICH FIELDS WE CAN QUERY THE DATABASE RESULTS BY
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
	var metadata
	var db_username

	/*THESE ARE SOME HARDCODED TEXT FIELDS THAT WILL APPEAR IN THE LEFT PANE. YOU CAN CHANGE THE TEXT OF THE FIELDS HERE*/
	var my_exp_var_name = 'MY EXPERIMENTS'
	var contains_ngs_data_var_name = 'CONTAINS NGS SEQUENCES'
	var contains_annotated_data_name = 'CONTAINS ANNOTATED DATA'
	var contains_recombination_name = 'RECOMBINATION TYPE'

	//this variable will create a  2-D list of ALL THE FILTERS WE WANT TO ALLOW THE USER TO FILTER EXPEIRMENTS BY.
	//first element of list:
		//When we first load the experiment metadata, it will parse the fields for each experiment and add to this list of filters
		//the first element of this variable will store each of the unique fields we observed
	//the second element of list:
		//for each filter/field, we will store all of the values observed by the metadata
	var list_of_filters = [
		[my_exp_var_name,contains_ngs_data_var_name,contains_annotated_data_name,contains_recombination_name,'VH:VL_PAIRED'],
		[[],[],['ANY'],[],[true,false]] //we will start with some default values for our current fields...
	]
	

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


	/******** Our main function ********/
	function createMetadataUI() {     		    	   
		//create a custom widget for datatables
		$.widget("custom.metadataselect",{
			options:{
				sort_tree:sort_tree_in_this_order,
				field_sort_order:default_field_sort_order,
				last_fields:default_field_sort_order_always_appear_last,
				include_seq_field_selection:false, //if true, then we add extra functionality to the widget
				load_on_start:true
			},
			_create:function(){
				this.element
					.addClass('expmetadata')
				if(this.options.load_on_start)				
					this.load_metadata()
			},

			_initialize_html:function(){

				html_str = "<div class='container-fluid section'>"+
								"<div class='row'>"+
									"<div class='metadata-wrapper sidebar1'>"+
										"<nav name='all-exp-filters' class='sidebar-menu sidebar1'>"+
											"<div>"+
												"<span class='sidebar-header'>"+
													"<strong>Experiment Table Filters</strong>"+
												"</span><hr />"+
												"<div class='filter-parent filter-annotation'>"+
													"<ul class ='nav nav-stacked parent'>"+					
													"</ul>"+
												"</div>"+
											"</div>"+
										"</nav>"		
				if(this.options.include_seq_field_selection){
					html_str+="<nav name='report-seq-data' class='sidebar2'>"+
									"<span class='sidebar-header'>"+
										"<strong>Output these fields</strong>"+
									"</span></hr>"+									
								"</nav>"
				}
				html_str+="<div name='exp-query-section' class='col-xs-12 main col-sm-12 col-md-12 main sidebar-2'>"+
								"<div>"+
									"<div class='space_down'>"+
										"<span class='space_down'>"+
											"<strong>Experiment Table</strong>"+
										"</span>"+
									"</div>"+
									"<div name='exptable-here' class='row'>"+									
									"</div>"+
								"</div>"+
							"</div>"+
					"</div></div></div>"
				return html_str
			},
			//this is a public function that can be called by the webpage. 
			//calling this function will query the dtabase and create the proper interface for working with the experiment table
			//make an ajax call to the database and query for all experiments
			load_metadata:function(proxy_path){				

				$.blockUI({ message: 'Searching metadata...' });
				if(proxy_path!=undefined)
					var args = [proxy_path]
				else
					var args = []
				params = {
					'module':'get_metadata_on_proxy',
					'args':args
				}

				var elem = this
				more_opts = {
					success: function (res) {
						$.unblockUI();
						elem.exp_metadata = res
						elem._initalize_app();
					},					
					error:function (jqXHR, textStatus, errorThrown) {
						$.unblockUI();
						alert(jqXHR.responseText);
						elem._clear_page();
					}					
				}				
				
				run_ajax_call(params,more_opts)		
			},
			//if querying the database doesnt work, then replace all html
			_clear_page:function(){
				this.element.find('.metadata-wrapper').html('<div class="alert alert-danger" role="alert">There was an error when processing the database query. Please report this to an administrator</div>')
			},
			//upon a successful ajax call, this will create a datatable summarizing all experimental metadata
			_initalize_app:function(){
				this.element.html(this._initialize_html())
				var info_from_db = this.exp_metadata
				if(!('metadata' in info_from_db)&&!('user_data' in info_from_db))
					alert('Invalid exeriment metadata object format')
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
				var dont_filter = ['SEQ_COUNT','DESCRIPTION','_id','FILENAMES']

				
				//metadata should be equal to the list of json values corresponding to the database. 
				//we need to parse this data to figure out how to populate the app 				
				//we also need to remove keys/objects we dont want to show
				//also set global variables for the app (schema fields, annotation methods)
					//these global variables will be used in later fields for populating the app 
				for (var i = 0; i<metadata.length;i++){
					analysis_schema_annotations[metadata[i]['_id']] = {}
					analysis_schema_datatypes[metadata[i]['_id']] = {}
					analysis_schema_annotations_html[metadata[i]['_id']] = {}
					exp_metadata_info[metadata[i]['_id']] = {}
					for (each_key in metadata[i]){
						if (each_key === 'ANALYSIS_SCHEMA'){
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
						else if (each_key === 'ANALYSES_COUNT'){
							exp_metadata_info[metadata[i]['_id']]['ANALYSES'] = Object.keys(metadata[i][each_key])
							all_analyses = all_analyses.concat(exp_metadata_info[metadata[i]['_id']]['ANALYSES']) //keep track of all observed analyses
							exp_metadata_info[metadata[i]['_id']]['ANALYSES_COUNT'] = metadata[i][each_key]
							for (each_analysis in metadata[i]['ANALYSES_COUNT']){
								all_recomb=MakeUnique(all_recomb.concat(all_recomb,Object.keys(metadata[i]['ANALYSES_COUNT'][each_analysis])))
							}
							exp_metadata_info[metadata[i]['_id']][contains_annotated_data_name] = ['ANY'].concat(Object.keys(metadata[i][each_key]))
							
						}
						else if(IgnoreFields.indexOf(each_key)===-1){
							exp_metadata_info[metadata[i]['_id']][each_key] = metadata[i][each_key]
														
							//ignore fields from those stored in variable dont_filter
							if (dont_filter.indexOf(each_key)===-1){
								//now check if this field is currently in the list_of_filters
								var pos = list_of_filters[0].indexOf(each_key)
								if(pos===-1){
									pos = list_of_filters[0].length
									list_of_filters[0].push(each_key)
									list_of_filters[1].push([])
								}
								//check if the value is a list:
								if (metadata[i][each_key] instanceof Array){
									//if its a list, then add these values to the current list of filters
									list_of_filters[1][pos] =  MakeUnique(list_of_filters[1][pos].concat(metadata[i][each_key]),true)
								}
								else if(metadata[i][each_key] === true || metadata[i][each_key] === false ){
									list_of_filters[1][pos] = [true,false]	
								}
								else{
									
									list_of_filters[1][pos].push(metadata[i][each_key])
									list_of_filters[1][pos] = MakeUnique(list_of_filters[1][pos],true)
								}
							}
						}
					}
					

					//lets add some extra fields to filter the experiment table by
					if('SEQ_COUNT' in exp_metadata_info[metadata[i]['_id']] && exp_metadata_info[metadata[i]['_id']]['SEQ_COUNT']>0)
						exp_metadata_info[metadata[i]['_id']][contains_ngs_data_var_name]= true
					else
						exp_metadata_info[metadata[i]['_id']][contains_ngs_data_var_name]= false
					
					if(exp_metadata_info[metadata[i]['_id']]['OWNERS_OF_EXPERIMENT'].indexOf(db_username))
						exp_metadata_info[metadata[i]['_id']][my_exp_var_name] = true
					else
						exp_metadata_info[metadata[i]['_id']][my_exp_var_name] = false					

					exp_metadata_list.push(exp_metadata_info[metadata[i]['_id']])
				}
				
				
				//make all analyses unique and then sort
				all_analyses = MakeUnique(all_analyses,true)
				all_recomb = MakeUnique(all_recomb,true)
				current_exp_set = exp_metadata_list
				
				list_of_filters[1][list_of_filters[0].indexOf(contains_annotated_data_name)] = list_of_filters[1][list_of_filters[0].indexOf(contains_annotated_data_name)].concat(all_analyses)
				list_of_filters[1][list_of_filters[0].indexOf(contains_recombination_name)] = list_of_filters[1][list_of_filters[0].indexOf(contains_recombination_name)].concat(all_recomb)

				console.log(list_of_filters)
				//Now draw the experiment datatable
				this.element.find('div[name="exptable-here"]').igreptable({
					header:[{'colname':'PROJECT_NAME'},{'colname':'EXPERIMENT_NAME'},{'colname':'READ_ACCESS'},{'colname':'OWNERS_OF_EXPERIMENT'},{'colname':'DESCRIPTION'},{'colname':'SEQ_COUNT'},{'colname':'ANALYSES'},{'colname':'_id','hidden':true}],
					data:exp_metadata_list,
					multiple_row_selection:true,
					use_cell_tooltips:false,
					use_row_qtips:true,
					qtip_header:'Experiment Metadata',
					max_table_height:'400px',
					/*min_col_width:'200px',*/
					/*wrap_text:false	*/
				})

				//adds annotation section to the html page
				//AddAnnotationFilter()
				//AddAnalysesRecombResults()

				//figure out ALL possible annotation fields based on users experiments. //create checkboxes for selecting specific fields
				//global_annotation= GetGlobalSchema(analysis_schema_annotations)

				//TableListener()
				//FindAnnotatedFields([])

			}

		});
	}	
	load_all_scripts()		
})(); // We call our anonymous function immediately

																																									
				