import * as d3 from 'd3'
// Top color bar

const BottomColor = ({
  heatmapsvg,
  bottomColorStartY,
  bottomColorEndY,
  colLabels,
  x,
  colorHeight2,
  color,
} = {}) => {

  // colors for each factor
  const botColorY = d3.scaleBand()
    .domain(["Factor Colors"])
    .range([bottomColorStartY, bottomColorEndY])
    .padding(0.1)
  heatmapsvg.append('g')

  // bottom factor color row
  heatmapsvg.selectAll()
    .data(colLabels)
    .join("rect")
      .attr("x", function(d) {return x(d)})
      .attr("y", botColorY("Factor Colors"))
      .attr("width", x.bandwidth())
      .attr("height", colorHeight2)
      .style("fill", function(d) {return color(d)})
      .style("stroke-width", 4)
      .style("stroke", "none")
  
}
export default BottomColor