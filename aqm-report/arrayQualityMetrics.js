// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSM631755", "Non-diabetic islets, rep1", "human islets, non-diabetic", "disease state: non-diabetic", "age: 47 yrs", "gender: male", "bmi (kg/m2): 27.7", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "2", "GSM631756", "Non-diabetic islets, rep2", "human islets, non-diabetic", "disease state: non-diabetic", "age: 33 yrs", "gender: male", "bmi (kg/m2): 22.9", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "3", "GSM631757", "Non-diabetic islets, rep3", "human islets, non-diabetic", "disease state: non-diabetic", "age: 47 yrs", "gender: male", "bmi (kg/m2): 28.4", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "4", "GSM631758", "Non-diabetic islets, rep4", "human islets, non-diabetic", "disease state: non-diabetic", "age: 54 yrs", "gender: male", "bmi (kg/m2): 23.1", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "5", "GSM631759", "Non-diabetic islets, rep5", "human islets, non-diabetic", "disease state: non-diabetic", "age: 76 yrs", "gender: female", "bmi (kg/m2): 25.9", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "6", "GSM631760", "Non-diabetic islets, rep6", "human islets, non-diabetic", "disease state: non-diabetic", "age: 77 yrs", "gender: female", "bmi (kg/m2): 23.8", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "7", "GSM631761", "Non-diabetic islets, rep7", "human islets, non-diabetic", "disease state: non-diabetic", "age: 73 yrs", "gender: female", "bmi (kg/m2): 22", "ctrl", "Gene expression data from non-diabetic human islets." ], [ "8", "GSM631762", "Type 2 diabetic islets, rep1", "human islets, diabetic", "disease state: type 2 diabetes", "age: 79 yrs", "gender: male", "bmi (kg/m2): 27.5", "diab", "Gene expression data from type 2 diabetic human islets." ], [ "9", "GSM631763", "Type 2 diabetic islets, rep2", "human islets, diabetic", "disease state: type 2 diabetes", "age: 76 yrs", "gender: male", "bmi (kg/m2): 26", "diab", "Gene expression data from type 2 diabetic human islets." ], [ "10", "GSM631764", "Type 2 diabetic islets, rep3", "human islets, diabetic", "disease state: type 2 diabetes", "age: 73 yrs", "gender: female", "bmi (kg/m2): 29", "diab", "Gene expression data from type 2 diabetic human islets." ], [ "11", "GSM631765", "Type 2 diabetic islets, rep4", "human islets, diabetic", "disease state: type 2 diabetes", "age: 75 yrs", "gender: female", "bmi (kg/m2): 26.5", "diab", "Gene expression data from type 2 diabetic human islets." ], [ "12", "GSM631766", "Type 2 diabetic islets, rep5", "human islets, diabetic", "disease state: type 2 diabetes", "age: 54 yrs", "gender: female", "bmi (kg/m2): 23.9", "diab", "Gene expression data from type 2 diabetic human islets." ], [ "13", "GSM631767", "Type 2 diabetic islets, rep6", "human islets, diabetic", "disease state: type 2 diabetes", "age: 66 yrs", "gender: male", "bmi (kg/m2): 23.1", "diab", "Gene expression data from type 2 diabetic human islets." ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
	var success = false;
	i = 0; 
	/* Some of this looping could already be cached in reportInit() */
	while( (!success) & (i < ssrules.length) ) {
	    selector = ssrules[i].selectorText;  // The selector 
            if (!selector) 
		continue; // Skip @import and other nonstyle rules
            if (selector == (".aqm" + reportObjId)) {
		success = true; 
		ssrules[i].style.cssText = cssText[0+status];
	    } else {
		i++;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
