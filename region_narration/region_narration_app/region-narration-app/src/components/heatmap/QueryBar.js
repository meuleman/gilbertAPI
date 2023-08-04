import * as d3 from 'd3'
// Bar chart to describe the query region

const QueryBarChart = ({
  heatmapsvg,
  tooltip,
  data,
  colLabels,
  x,
  y,
  queryStartY,
  queryEndY,
  mousemove,
  mouseleave,
  mouseoverLabel,
  mousemoveLabel,
  mouseleaveLabel,
  queryLabelX,
  queryLabelY,
  plotStartX,
  plotEndX,
  linkClick,
  queryButtonX,
  browserButtonWidth,
  detailTargetFactors,
} = {}) => {
  const queryData = data[0]
  const queryDataFlat = (data.slice(0,1).flatMap((region) => {
    const rank = region.rank
    const percentiles = region.percentiles
    const coords = region.coordinates
    const score = region.score
    const regionInd = region.region_ind

    return (
      percentiles.map((p, index) => {
        return ({
          percentile: p,
          factor: colLabels[index],
          rank: rank,
          coords: coords,
          score: score,
          regionInd: regionInd,
          factorInd: index,
          originalFactorInd: queryData.factor_order[index]
        })
      })
    )
  }))

  // mouseover for query region bar chart
  const mouseoverQuery = function (d, cellData) {
    const factor = cellData.factor
    const percentile = Math.round((cellData.percentile * 100) * 10) / 10
    tooltip
      .html(factor + ": " + percentile + "%")
      .style('opacity', 1)
      .style('z-index', 1)
    d3.select(this)
      .style("stroke", "black")
      .style("stroke-width", "1")
  }

  // query region bar chart y axis
  const queryY = d3.scaleLinear()
    .domain([1, 0])
    .range([queryStartY, queryEndY])
  heatmapsvg.append('g')
    .style('font', '8px sans-serif')
    .attr("transform", "translate(" + plotStartX + ", 0)")
    .call(d3.axisLeft(queryY).tickValues([1.0, 0.5, 0]).tickFormat(d3.format(".0%")))

  // query region bar chart label
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "end")
    .attr("y", queryLabelY)
    .attr("x", queryLabelX)
    // .attr("transform", "rotate(-90)")
    .style('font', '13px sans-serif')
    .text("Query Region")
    .on('mouseover', () => mouseoverLabel("Ranks for factors in the query region"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)

  // add roiURL button
  heatmapsvg.selectAll()
    .data([queryData])
    .join("rect")
      .attr("x", queryButtonX)
      .attr("y", queryY(0.5) + y.bandwidth() / 2 + 2.5)
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
        "browser": "roiURL"
      }})
      .on("click", (d) => {linkClick(d)})


  heatmapsvg.selectAll()
    .data([queryData])
    .join("text")
      .attr("x", queryButtonX + browserButtonWidth / 2)
      .attr("y", queryY(0.5) + 3 * y.bandwidth() / 2 - 2.5)
      .attr("text-anchor", "middle")
      .style('font', '15px sans-serif')
      .style("cursor", "pointer")
      .text("r")
      .datum(function(d) {return {
        "coords": d.coordinates,
        "browser": "roiURL"
      }})
      .on("click", (d) => {linkClick(d)})


  // add query region bar chart
  heatmapsvg.selectAll()
  .data(queryDataFlat)
  .join("rect")
    .attr("x", function(d) {return x(d.factor)})
    .attr("y", function(d) {if (d.percentile > 0.5) {return queryY(d.percentile)} else {return queryY(0.5)}})
    .attr("width", x.bandwidth())
    .attr("height", function(d) {if (d.percentile > 0.5) {return queryY(0.5) - queryY(d.percentile)} else {return queryY(d.percentile) - queryY(0.5)}})
    .style("fill", function(d) {if (detailTargetFactors.includes(d.originalFactorInd)) {return "#555555"} else {return "#999999"}})
    .style("stroke-width", 4)
    .style("stroke", "none")
    .on('mouseover', mouseoverQuery)
    .on('mousemove', mousemove)
    // .on('mouseleave', mouseleave)
    .on('mouseleave', function() {mouseleave.bind(this)("none")})
    // .on('click', mouseclick)

  // line at 50% in query region bar chart
  heatmapsvg
    .append("line")
      .style("stroke", "black")
      .style("stroke-width", 0.5)
      .attr("x1", plotStartX)
      .attr("y1", queryY(0.5))
      .attr("x2", plotEndX)
      .attr("y2", queryY(0.5))
  return {
    queryY: queryY
  }
}
export default QueryBarChart