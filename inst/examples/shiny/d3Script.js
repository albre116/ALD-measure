<script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>

<script>

console.log('I ran!')

var homzScale = d3.scale.linear()
            .domain([0,1])
            .range(['#4daf4a', '#e41a1c'])

var freqScale = d3.scale.linear()
            .domain([0,.3])
            .range(['#4daf4a', '#e41a1c'])

//This is a crazy hacky way to do this, but it works!
//It waits to execute the javascript code 200 miliseconds. This is needed because shiny loads the javascript before it loads the table
//and thus has nothing to color.
window.setInterval( function(){

console.log($('.tab-pane').filter("[data-value='Testing!']").hasClass('active'))
if(false){
     $('tbody tr td:nth-child(5)').each(function() {
      var cellValue = $(this).text();

      $(this).css('background-color', homzScale(cellValue));
     })

    $('tbody tr td:nth-child(4)').each(function() {
      var cellValue = $(this).text();

      $(this).css('background-color', freqScale(cellValue));
     })
  }


}
, 200);

</script>
