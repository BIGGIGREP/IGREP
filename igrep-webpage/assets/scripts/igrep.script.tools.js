var host_location = String(window.location.host)
var pos_location = location.protocol+"//"+host_location
//this will return some default settings for running ajax calls to the web-server
//params should be an object with the following keys: module (value = string), args value = list, and kwargs value= object
var mongoproxy_path


//params:object with following keys (only module is mandatory)
	//module = name of method to run 
	//args = args for method
	//kwargs = keyword args for method
//more_paraams => object to pass into remainder of ajax call
function run_ajax_call(params,more_params){
	var opts = {
			url: encodeURI(pos_location),
			type: "POST",
			contentType: "application/json",
			dataType: "json",
			data: JSON.stringify(params)		
	}
	
	if(more_params != undefined)
		$.extend( opts, more_params);	
	console.log(more_params)
	console.log(opts)	
	// DO ajax here
	xhr = $.ajax(opts);
}
function get_mongoproxy_address(on_complete_function) {	
	if (mongoproxy_path != undefined) {
		alert('itsalreadydefined!')
		if (on_complete_function)
			on_complete_function(mongoproxy_path)
	}
	else {
		var params = {
			module: 'get_default_mongoproxy_address', //run this function
			args: [], //use this proxy path to database
			kwargs: {}
		}
		var opts = get_ajax_params(params)
		//opts.async=false	
		opts.success = function (res) {
			mongoproxy_path = res
		}
		opts.error = function (jqXHR, textStatus, errorThrown) {
			alert(jqXHR.responseText);
			alert('Could not find the default monogproxy address so will not be able to access the database')
			mongoproxy_path = ''
		}
		opts.complete = function () {
			if (on_complete_function)
				on_complete_function(mongoproxy_path)
		}
		xhr = $.ajax(opts);
	}
}
//window.onload = get_mongoproxy_address()

function createCookie(name,value,days) {
	if (days) {
		var date = new Date();
		date.setTime(date.getTime()+(days*24*60*60*1000));
		var expires = "; expires="+date.toGMTString();
	}
	else var expires = "";
	document.cookie = name+"="+value+expires+"; path=/";
}

function readCookie(name) {
	var nameEQ = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++) {
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
	}
	return null;
}

function eraseCookie(name) {
	createCookie(name,"",-1);
}