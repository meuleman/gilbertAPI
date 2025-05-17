// epilogos buttons for similar regions

const SearchFromHereButtons = ({
  heatmapsvg,
  data,
  searchButtonX,
  searchButtonWidth,
  y,
  search,
  searchedByFactor,
} = {}) => {
  let simData = null
  if (searchedByFactor) {
    simData = data
  } else {
    simData = data.slice(1)
  }

  // search from here button
  heatmapsvg.selectAll()
    .data(simData)
    .join("rect")
      .attr("x", searchButtonX)
      .attr("y", function(d) {return y(d.rank)})
      .attr("width", searchButtonWidth)
      .attr("height", y.bandwidth())
      .style("fill", "#E0E0E0")
      .style("stroke-width", 1)
      .style("stroke", "black")
      .style("rx", 1)
      .style("ry", 1)
      .style("cursor", "pointer")
      .datum(function(d) {return {
        "regionCoords": d.coordinates
      }})
      .on("click", (d) => {
        search(d.srcElement.__data__.regionCoords)
      })
  heatmapsvg.selectAll()
    .data(simData)
    .join("text")
      .attr("x", searchButtonX + searchButtonWidth / 2)
      .attr("y", function(d) {return y(d.rank) + y.bandwidth() - 5})
      .attr("text-anchor", "middle")
      .style('font', '10px sans-serif')
      .style("cursor", "pointer")
      .text("Go!")
      .datum(function(d) {return {
        "regionCoords": d.coordinates
      }})
      .on("click", (d) => {
        search(d.srcElement.__data__.regionCoords)
      })
}
export default SearchFromHereButtons