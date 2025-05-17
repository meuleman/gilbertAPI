// index buttons for query region

const QueryIndexButton = ({
  heatmapsvg,
  data,
  queryButtonX,
  queryY,
  linkClick,
  y,
  browserButtonWidth
} = {}) => {
  const queryData = data[0]

  // add query index button
  heatmapsvg.selectAll()
    .data([queryData])
    .join("rect")
      .attr("x", queryButtonX)
      .attr("y", queryY(0.5) - 3 * y.bandwidth() / 2 - 2.5)
      .attr("width", browserButtonWidth)
      .attr("height", y.bandwidth())
      .style("fill", "#E0E0E0")
      .style("stroke-width", 1)
      .style("stroke", "black")
      .style("rx", 1)
      .style("ry", 1)
      .style("cursor", "pointer")
      .datum(function(d) {return {
        "coords": d.coordinates,
        "browser": "index"
      }})
      .on("click", (d) => {linkClick(d)})
  heatmapsvg.selectAll()
    .data([queryData]) 
    .join("text")
      .attr("x", queryButtonX + browserButtonWidth / 2)
      .attr("y", queryY(0.5) - y.bandwidth() / 2 - 7.5)
      .attr("text-anchor", "middle")
      .style('font', '15px sans-serif')
      .style("cursor", "pointer")
      .text("i")
      .datum(function(d) {return {
        "coords": d.coordinates,
        "browser": "index"
      }})
      .on("click", (d) => {linkClick(d)})
}
export default QueryIndexButton