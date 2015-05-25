//The following are d3 scaling functions.
//We are using a linear scale (scale.linear())
//There are 4 values in both the domain and range representing
//the four different colors represented in the color gradient
var homzScale = d3.scale.linear()
    .domain( [0,         .33,        .66,        1])
    .range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])

var freqScale = d3.scale.linear()
.domain( [0,         .33,        .66,        1])
.range(['#66bd63', '#fee08b', '#f46d43', '#d73027'])

//This is a quick function to take a two value array of max and min values
//and spit out a four element array of evenly spaced values between those two points.
//If this was in R it would be easier (and may be easier in JS, I just dont know how.)
var domainGen = function(ar){
//Javascript indexes starting at 0, so this is the first value, or min.
var start = ar[0]
var end = ar[1]
var range = end - start
//return the finished array.
return [start, start + (0.33*range), start + (0.66*range), end]
}

//Waits to execute the javascript code 200 miliseconds. This is needed because shiny loads the javascript before it loads the table
//and thus has nothing to color.
window.setInterval(function() {

//Reach in and grab the max and min values from the hidden output renderText()
//See lines ~235-236 for the location of these values.
var maxFreq = parseFloat(d3.select('#maxVal_f').text())
var minFreq = parseFloat(d3.select('#minVal_f').text())
freqColorRange = [minFreq, maxFreq]
//modify the domain of the frequency color scale from earlier with these newly
//obtained max and min values
freqScale.domain(domainGen(freqColorRange))

//Now we repeat this with the homz values
var maxHomz = parseFloat(d3.select('#maxVal_h').text())
var minHomz = parseFloat(d3.select('#minVal_h').text())
homzColorRange = [minHomz, maxHomz]
homzScale.domain(domainGen(homzColorRange))

//Now we use d3 to select all ('selectAll') of the frequency values by
//specifically targeting the 4th column of the html table.
//Then it sets the css style (.style) background-color by first grabbing the
//particular cell's value (d3.select(this).text()) and then putting that value
//into the color scaling function we set up in the above code.
//Set the colors for the allele freq
d3.selectAll('#asf_table tbody tr td:nth-child(4)')
.style('background-color', function() {
//grab the cells value
var cellValue = d3.select(this).text();
//return for the background-color value the result of the value run through color scale.
return (freqScale(cellValue))
})

//we then repeat this for the homz values
d3.selectAll('#asf_table tbody tr td:nth-child(5)')
.style('background-color', function() {
var cellValue = d3.select(this).text();
return (homzScale(cellValue))
})
//this is the amount of miliseconds between runs of the function. 
}, 200);
