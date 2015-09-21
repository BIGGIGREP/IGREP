/*This is a wrapper script for making tables in our webpages using jquery&datatable. We use these tables in many of our html pages. 
	First attempt at widget making...
	It will be dependent on:
		1) jquery (version 1.11.3)
		2) datatable (version 1.10.9)
		3) qtip (version 2.2.1)



	<script src="//code.jquery.com/jquery-1.11.3.js"></script>
	<script src="//cdn.datatables.net/1.10.9/js/jquery.dataTables.min.js"></script>	
	<script src="//cdn.jsdelivr.net/qtip2/2.2.1/jquery.qtip.min.js"></script>

	<script type="text/javascript" src="/scripts/jquery-1.11.3.min.js"></script> 
    <script type="text/javascript" src="/scripts/datatables.min.js"></script>
    <script type="text/javascript" src="/scripts/jquery-ui-1.11.4.js"></script>
    <script type="text/javascript" src="/scripts/qtip/jquery.qtip.js"></script> 

    <!-- <script type="text/javascript" src="/scripts/exptable.js"></script>-->
    <link rel="stylesheet" type="text/css" href="/styles/datatables.css" >
	<link rel="stylesheet" type="text/css" href="/styles/jquery.qtip.min.css">
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

/******** Our main function ********/
function createTableWidget() { 
    css_link = $("<link>",{ 
    	rel: "stylesheet", 
		type: "text/css", 
		href: "/styles/immunogrep_stylesheet.css"
	})
	css_link.appendTo('head');
	    
    $(document).ready(function() { 
    	/****Load the following scripts***/
    	setTimeout(function() {alert('hello');},1250);
    	$('#table_id').DataTable()
		alert('loaded!')		
    });
}

loadDataTable();
})(); // We call our anonymous function immediately