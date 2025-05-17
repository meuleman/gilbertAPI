import * as d3 from 'd3'
// Top color bar

const TopColor = ({
  heatmapsvg,
  topColorStartY,
  topColorEndY,
  plotStartX,
  colLabels,
  x,
  colorHeight1,
  color

} = {}) => {
  // colors for each factor
  const topColorY = d3.scaleBand()
    .domain(["Factor Colors"])
    .range([topColorStartY, topColorEndY])
    .padding(0.1)
  heatmapsvg.append('g')
    .attr("transform", "translate(" + plotStartX + ", 0)")
    .call(d3.axisLeft(topColorY).tickValues([]).tickSize(0))

  // top factor color row
  heatmapsvg.selectAll()
    .data(colLabels)
    .join("rect")
      .attr("x", function(d) {return x(d)})
      .attr("y", topColorY("Factor Colors"))
      .attr("width", x.bandwidth())
      .attr("height", colorHeight1)
      .style("fill", function(d) {return color(d)})
      .style("stroke-width", 4)
      .style("stroke", "none")
}
export default TopColor