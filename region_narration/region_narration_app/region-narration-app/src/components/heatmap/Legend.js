import * as d3 from 'd3'
// legend for heatmap

const Legend = ({
  heatmapsvg,
  heatmapColor,
  heatmapStartY,
  heatmapEndY,
  legendEndX,
  legendStartX,
  legendWidth
} = {}) => {
  // colorbar for heatmap legend
  const heatmapLegend = heatmapsvg.append("linearGradient")
    .attr("id", "heatmap-legend")
    .attr("x1", "0%")
    .attr("y1", "100%")
    .attr("x2", "0%")
    .attr("y2", "0%")
  heatmapLegend
    .append("stop")
      .attr("offset", "0%")
      .attr("stop-color", heatmapColor(0))
  heatmapLegend
    .append("stop")
      .attr("offset", "100%")
      .attr("stop-color", heatmapColor(1))

  // heatmap legend y axis
  const legendY = d3.scaleLinear()
    .domain([1, 0])
    .range([heatmapStartY, heatmapEndY])
  heatmapsvg.append('g')
    .style('font', '8px sans-serif')
    .attr("transform", "translate(" + (legendEndX) + ", 0)")
    .call(d3.axisRight(legendY).tickValues([0, 0.25, 0.50, 0.75, 1.0]).tickFormat(d3.format(".0%")))

  // add heatmap colorbar legend
  heatmapsvg.append("rect")
    .attr("x", legendStartX)
    .attr("y", legendY(1))
    .attr("height", legendY(0) - legendY(1))
    .attr("width", legendWidth)
    .style("fill", "url(#heatmap-legend)")
  
  return {
    legendY: legendY
  }
}
export default Legend