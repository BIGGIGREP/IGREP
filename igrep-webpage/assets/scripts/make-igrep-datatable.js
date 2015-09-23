/*This is a wrapper script for making tables in our webpages using jquery&datatable. We use these tables in many of our html pages. 
	First attempt at widget making...
	It will be dependent on:
		1) jquery (version 1.11.4)
		2) jquery ui 1.11.4
		3) datatable (version 1.10.9)
		4) qtip (version 2.2.1)

	In this first attempt, I was not able to get the 'load javascript functions'
	
	IMPORTANT: BECAUSE I COULD NOT GET THE LOAD-JAVASCRIPT FUCNTIONS TO WORK, I ASSUME THAT THE HTML PAGES ALREADY IMPORT THE REQUIRED DEPENDENCIES!!!
	
	
	
*/

(function() {

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

	/*
	function loadDataTable(){
		//datatables
		loadScript("//cdn.datatables.net/1.10.9/js/jquery.dataTables.min.js",datatableImported)	
	}
	function datatableImported(){
		var css_link = $("<link>", { 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: "/styles/datatables.css"
		})	
		css_link.appendTo('head'); 
		loadQTip()
	}

	function loadQTip(){
		//qtip
		loadScript("//cdn.jsdelivr.net/qtip2/2.2.1/jquery.qtip.min.js",qtipImported)    	
	}
	function qtipImported(){
		css_link = $("<link>",{ 
	    	rel: "stylesheet", 
			type: "text/css", 
			href: "/styles/jquery.qtip.min.css"
		})
		css_link.appendTo('head');
		//now call the main function
		createTableWidget();
	}
	*/

	/******** Our main function ********/
	function createTableWidget() {     		    	   
	}

	//create a custom widget for datatables
	$.widget("custom.igreptable",{
		//default options
		options:{
			header:[],//header=>list of dictioary containing three keys:
			//colname
			//hidden: true/false
			//width: column width
			data:[],//data => list of json docs. key corresponds to headr/colum name
			multiple_row_selection:true,
			use_cell_tooltips:false,
			use_row_qtips:false,
			qtip_header:'',
			row_tooltips:[],
			class_table_name:'',
			max_table_height:'500px',
			min_col_width:'',
			max_col_width:'',
			wrap_text:false,
			
		},

		_create:function(){
			var datatable_column_visibilty = []//this variable will control visibilty to datatable			

			this.element
				.addClass('igreptable')				
			if(this.options.class_table_name!='')
				this.element.addClass(this.options.class_table_name)
			var header = this.options.header;
			
			//initialize string for making header row
			table_hd = "<tr>"
			//parse through header variable: based on keys, update datacolumn visibility field
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
						datatable_column_visibilty.push(null)
				}
				table_hd+="<td class='"+headername+"'>"+headername+"</td>"
			}
			table_hd+="</tr>"
			var global_table_header = header

			//make html code for the necessary elements we will want
			//we will inject html TABLE code and inject inputs/buttons
			html_code = "<table class='display'>" +
			"<thead>"+table_hd+"</thead>"+
			"<tbody></tbody>"+
			"</table>"+
			"<div align='right'>"+
				"<input type='button' class='all-row-select multirowselect not-multiple' value='Select all filtered rows'>"+
				"<input type='button' class='all-row-deselect multirowselect not-multiple' value='Deselect all rows'>"+
				"<div class='table-errors' style='color:red;'></div>"
			"</div>"
			
			
			this.element.html(html_code);

			//settings used for making datatable (using datatable 1.9...)
			/*datatable_parameters = {
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
			}*/

			//settings uses for making datatable (using datatable 1.10)
			
			datatable_parameters = {
				'ordering':true,
				"scrollX": "100%",
				"scrollCollapse": true,
				"jQueryUI": true,
				"pagingType": "full_numbers",
			//	"autoWidth": false, //does weird things to the scrollbar...it moves the horizontal scrollbar back after sorting..
				"autowidth":false,
			//	"scrollXInner": "100%",
			//	"stateSave": true,
			//	"deferRender": true,
				"retrieve": true,
				"columns": datatable_column_visibilty,
			}

			if (this.options.max_table_height!=''){
				if(typeof this.options.max_table_height!='string')
					this.options.max_table_height = String(this.options.max_table_height)+'px'
				//datatable_parameters["sScrollY"] = this.options.max_table_height
				datatable_parameters["scrollY"] = this.options.max_table_height
			}

			//add this key to the variable above if we want to use cell tool tips
			//basically each of the values in the table will appear as tooltips also
			if (this.options.use_cell_tooltips==true){
				//datatable_parameters['fnCreatedRow']=function(nRow,aData,iDataIndex){
				datatable_parameters['createdRow']=function(nRow,aData,iDataIndex){
					for (var i =0 ;i<aData.length;i++){
						$('td:eq('+i.toString()+')', nRow).attr('title',aData[i]);
					}
				}
			}

			this.myDataTable = this.element.find('table').DataTable(datatable_parameters)
			var elem = this
			
			if (this.options.multiple_row_selection){
				this.element.find('.multirowselect').removeClass('not-multiple')

				//Allow the HTML table to be selectable using the jquery selectable	
				this.element.find('table tbody').selectable({			
					items:'tr',
					stop: function( event, ui ) {					
						var current_rows
						if(event.ctrlKey){//the user is pressing down the ctrl key for selection (so perhaps we want to save all hihglihged rows in all pages 
							//find all selected rows that are ONLY in the filtered table 
							current_rows = elem.myDataTable.$('tr.ui-selected', {"filter": "applied"}); 
						}
						else{
							//only find selected rows that IN THE CURRENT PAGE!							
							current_rows = elem.myDataTable.$('tr.ui-selected', {"page":"current"});
							
						}				
						//remove all rows in the entire datatable (including rows not filtered)
						elem.myDataTable.$('tr').removeClass('ui-selected');
						//now re-select the rows from first statement (selected rows only in filtered table) 				
						current_rows.addClass('ui-selected');		
						//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
						elem.myDataTable.$('tr').first().trigger('click')						
					}
				}).disableSelection()
				
				//add javascript/jquery to allow multiple selection of fields
				this.element.find('.all-row-deselect').on('click',function(){					
					//FIND ALL ROWS in the datatable var (must use datatable() function to include non -filtered rows
					elem.myDataTable.$('tr').removeClass('ui-selected');	
					//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
					elem.myDataTable.$('tr').first().trigger('click')
				});				


				this.element.find('.all-row-select').on('click',function(){
					//FIND ROWS THAT ARE FILTERED
					var filteredrows = elem.myDataTable.$('tr', {"filter": "applied"});					
					elem.myDataTable.$('tr').removeClass('ui-selected'); //remove class from ALL rows in table
					filteredrows.addClass('ui-selected');//add class to filtered rows only 
					//trigger a click event so that other HTML or other functions can detect SELECTING of datatable rows (not just clicking of rows)
					elem.myDataTable.$('tr').first().trigger('click')
				});		
			}
			else{ //ALLOW only one row to be selected
						
				this.element.find('tbody tr', function() { // tbody tr').click( function( event ) {												
					if ( $(this).hasClass('ui-selected') ) {						
						$(this).removeClass('ui-selected');						
					}                                              
					else {                             
						elem.myDataTable.$('tr').removeClass('ui-selected'); //remove class from ALL rows in table
						$(this).addClass('ui-selected');  												
					}; 								 					
				});
			}
				
			//NOW populate the actual table wiht the data variable 			
			this.modify_data()

			//Use jquery to change some css styles provided by user 
			var css_styles = {}
			if (this.options.min_col_width!='')
				css_styles['min-width'] = String(this.options.min_col_width).replace('px','')+'px' //use .replace just incease user passed in a number. so, i.e.,this style will work if user said 200px OR 200 
			if (this.options.max_col_width!='')
				css_styles['max-width'] = String(this.options.max_col_width).replace('px','')+'px' //use .replace just incease user passed in a number. so, i.e.,this style will work if user said 200px OR 200 
			
			//modfiy table header css
			//myDataTable.$('thead td').css(css_styles)
			theadstyles = css_styles
			theadstyles['padding-right']='25px'
			this.element.find('table thead td').css(theadstyles)
												
			if (this.options.wrap_text){
				css_styles['white-space'] = 'normal' 
				css_styles['word-wrap'] = 'break-word'						
			}
			else{
				css_styles['white-space'] = 'nowrap' 
				css_styles['text-overflow'] = 'ellipsis'			
			}		
				
			this.element.find('td').css(css_styles)
			

			//create a callback trigger so that users can respond to interaction with database
			//listen for any time a row is pressed, when this happens then send a callback
			$(this.element).find('table').on('click','tr',function(){
				selected_data = elem.get_table_data()				
				//create a callback called "select"
				elem._trigger("change",null, {rows:selected_data[0],data:selected_data[1]})				
			});


		},

		//public function for getting the table element
		get_table:function(){
			return this.element.find('table')
		},

		//public function for returning the selected fields
		get_table_data:function(){
			//DO THE FOLLOWING SEARCH: 
			//search for all table elements
			//then for each table only select tables that have a row with calss ui-selected
			//THEN go back a step and get the parent (which is the original table)
			//ONCE we do that then we can retrieve the original datatable 
			var dT = this.myDataTable

			//get selected rows that ONLY Mmatch the dattable filter 
			var anselected=dT.$('tr.ui-selected', {"filter": "applied"}); // only return selected rows that are filtered					
			var rows = [];
			var data = [];
			//report results here (go through all selected rows and report the results from the table)
			for (i = 0; i<anselected.length;i++){					
				var iRow=dT.row(anselected[i]).index();     					
				var aData=dT.row(anselected[i]).data();// fnGetData(iRow);   		
				rows.push(iRow)
				data.push(aData)				
			}
			var rowData= [rows,data]			
			return rowData
		},
			
							
		//function for adding rows to datatable. can be called externally by user
		modify_data:function(data,qtips,qtipheader){			
			
			if (data!=undefined)
				this.options.data = data
			
			if (qtips!=undefined){
				this.options.row_tooltips = qtips			
				this.options.use_row_qtips = true
			}

			if (qtipheader!=undefined)
				this.options.qtip_header = qtipheader

			//NOW populate the actual table wiht the data variable 
			var rows_of_data = []			
			//data =>list, each element of list is a dict where key = colname
			//header => list, each element is list. key in each element called colname refers to column name 			
			var header = this.options.header
			var colname
			
			for (row_num= 0; row_num<this.options.data.length;row_num++){
				data_row = this.options.data[row_num]
				row_info = []
				for (i = 0; i<header.length;i++){
					colname = header[i]['colname']
					if (colname in data_row)
						if(data_row[colname].constructor==Array){					
							row_info.push(data_row[colname].join(','));
						}
						else{
							row_info.push(String(data_row[colname]));
						}						
					else
						row_info.push('')
				}
				rows_of_data.push(row_info );
			}	
			
			this._render_table(rows_of_data)

			//the user wants to use qtips for hovering and showing information at each row 
			if (this.options.use_row_qtips){
				//The user passed DID not pass in list of strings to use as q tips
				//So we will just use json dumps on the current data and display those results as q tips 
				if( this.options.row_tooltips.length==0){
					for(var row_num = 0; row_num < this.options.data.length; row_num+=1){
						//use json dumps and an indentation of 4 spaces to present the data in a prety format 
						//use <pre> to make sure html reads in \t \n and spaces in toolstips
						this.options.row_tooltips.push('<pre>'+JSON.stringify(this.options.data[row_num],null,4)+'</pre>')
					}	
				}
				//Now generate the qtips on the provided datatable 				
				this._add_qtip_title()
			}
		},

		//private function for rendering datatable
		_render_table:function(rows_of_data){			
			var elem = this
			var timer_id;
			elem.myDataTable.clear()// var dT = datatable_var;			
			elem.myDataTable.rows.add(rows_of_data).draw();
			//add a listening event to each time window is resized. adjust table accordingly
			
			$(window).resize(function() {								
				clearTimeout(timer_id);
				timer_id = setTimeout(function() {    		        						
					elem.myDataTable.columns.adjust()//fnAdjustColumnSizing();
				}, 300);
			});
			
			$(window).trigger( "resize" );//ADJUST TABLE HEIGHT AND COLUMNS SEE window resize function below 				
		},

		_refresh_table_widths:function(){
			$(window).trigger( "resize" );//ADJUST TABLE HEIGHT AND COLUMNS SEE window resize function below 			
		},

		//tool_tip_text_array => a list of text to associate as qtips to each of the rows 
		_add_qtip_title:function(){
			
			var elem = this		
			//$('#'+tableid+' tbody tr').each(function(){{	--> DO NOT USE THIS, USE DATATABLE ROWS. TABLEID WILL ONLY LIST THE FIRST TEN ROWS
			
			//first add a title attribute to each row			
			elem.myDataTable.$('tr').each(function(){								
				$(this).attr('title',elem.options.row_tooltips[elem.myDataTable.row(this).index()])
			});	

			//finally associate a qtip jquery to all of the rows in the datatable
			elem.myDataTable.$('tr').qtip({
				content:{	
					title: elem.options.qtip_header,
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
					classes:"expqtipcontent"

				}
			});	
								
		}
		

	});
		
		
})(); // We call our anonymous function immediately

