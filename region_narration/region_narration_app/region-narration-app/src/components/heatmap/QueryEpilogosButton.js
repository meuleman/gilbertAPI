// epilogos buttons for query region

const QueryEpilogosButton = ({
  heatmapsvg,
  data,
  queryButtonX,
  queryY,
  linkClick,
  y,
  browserButtonWidth
} = {}) => {
  const queryData = data[0]

  // add query epilogos button
  heatmapsvg.selectAll()
    .data([queryData])
    .join("rect")
      .attr("x", queryButtonX)
      .attr("y", queryY(0.5) - y.bandwidth() / 2)
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
        "browser": "epilogos"
      }})
      .on("click", (d) => {linkClick(d)})
  heatmapsvg.selectAll()
    .data([queryData])
    .join("text")
      .attr("x", queryButtonX + browserButtonWidth / 2)
      .attr("y", queryY(0.5) + (y.bandwidth() / 2) - 5)
      .attr("text-anchor", "middle")
      .style('font', '15px sans-serif')
      .style("cursor", "pointer")
      .text("e")
      .datum(function(d) {return {
        "coords": d.coordinates,
        "browser": "epilogos"
      }})
      .on("click", (d) => {linkClick(d)})
}
export default QueryEpilogosButton