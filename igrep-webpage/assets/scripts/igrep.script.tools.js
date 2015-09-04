var host_location = String(window.location.host)
var pos_location = location.protocol+"//"+host_location
//this will return some default settings for running ajax calls to the web-server
//params should be an object with the following keys: module (value = string), args value = list, and kwargs value= object
var mongoproxy_path
function get_ajax_params(params){
	var opts = {
		url: encodeURI(pos_location),
		type: "POST",
		contentType: "application/json",
		dataType: "json",
		data: JSON.stringify(params)
	}
	return opts
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
window.onload = get_mongoproxy_address()
