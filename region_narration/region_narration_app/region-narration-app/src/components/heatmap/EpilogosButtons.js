// epilogos buttons for similar regions

const EpilogosButtons = ({
  heatmapsvg,
  data,
  epilogosButtonX,
  y,
  browserButtonWidth,
  linkClick,
  searchedByFactor,
} = {}) => {
  let simData = null
  if (searchedByFactor) {
    simData = data
  } else {
    simData = data.slice(1)
  }

  // add epilogos buttons
  heatmapsvg.selectAll()
    .data(simData)
    .join("rect")
      .attr("x", epilogosButtonX)
      .attr("y", function(d) {return y(d.rank)})
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

  // console.log(simData)
  heatmapsvg.selectAll()
    .data(simData)
    .join("text")
      .attr("x", epilogosButtonX + browserButtonWidth / 2)
      .attr("y", function(d) {return y(d.rank) + y.bandwidth() - 5})
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
export default EpilogosButtons