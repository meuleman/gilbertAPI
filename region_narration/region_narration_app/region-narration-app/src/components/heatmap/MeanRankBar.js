import * as d3 from 'd3'
// Bar chart to describe the mean ranks of similar regions

const MeanRankBarChart = ({
  heatmapsvg,
  tooltip,
  meanStartY,
  meanEndY,
  plotStartX,
  meanLabelX,
  meanLabelY,
  mouseoverLabel,
  mousemoveLabel,
  mouseleaveLabel,
  data,
  x,
  detailTargetFactors,
  mousemove,
  mouseleave,
  searchedByFactor,
  colLabels
} = {}) => {
  const dataFlat = (data.flatMap((region) => {
    const rank = region.rank
    const percentiles = region.percentiles

    return (
      percentiles.map((p, index) => {
        return ({
          percentile: p,
          rank: rank,
          factorInd: index,
        })
      })
    )
  }))

  let heatmapData = ''
  if (searchedByFactor) {
    heatmapData = dataFlat
  } else {
    heatmapData = dataFlat.slice(34)
  }

  // const factorMean
  let factorSums = new Array(colLabels.length).fill(0)
  heatmapData.forEach((d) => {
    factorSums[d.factorInd] += d.percentile
  })
  const factorMeans = factorSums.map((value, index) => {
    return {
      mean: value / (heatmapData[heatmapData.length - 1].rank + 1),
      factor: colLabels[index],
      originalFactorInd: data[0].factor_order[index]
    }
})

  // mouseover for factor mean bar plot
  const mouseoverMean = function (d, cellData) {
    const mean = Math.round((cellData.mean * 100) * 10) / 10
    const factor = cellData.factor
    tooltip
      .html(factor + " Mean: " + mean + "%")
      .style('opacity', 1)
      .style('z-index', 1)
    d3.select(this)
      .style("stroke", "black")
      .style("stroke-width", "1")
  }

  // mean rank region bar chart y axis
  const factorMeanY = d3.scaleLinear()
    .domain([1, 0])
    .range([meanStartY, meanEndY])
  heatmapsvg.append('g')
    .style('font', '8px sans-serif')
    .attr("transform", "translate(" + plotStartX + ", 0)")
    .call(d3.axisLeft(factorMeanY).tickValues([1.0, 0.5, 0]).tickFormat(d3.format(".0%")))
  
  // add significance bar chart label
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "end")
    .attr("y", meanLabelY)
    .attr("x", meanLabelX)
    // .attr("transform", "rotate(-90)")
    .style('font', '13px sans-serif')
    .text("Rank Enrichment")
    .on('mouseover', () => mouseoverLabel("Average of factor ranks in similar regions"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)

  // factor mean bar chart
  heatmapsvg.selectAll()
    .data(factorMeans)
    .join("rect")
      .attr("x", function(d) {return x(d.factor)})
      .attr("y", function(d) {if (d.mean > 0.5) {return factorMeanY(d.mean)} else {return factorMeanY(0.5)}})
      .attr("width", x.bandwidth())
      .attr("height", function(d) {if (d.mean > 0.5) {return factorMeanY(0.5) - factorMeanY(d.mean)} else {return factorMeanY(d.mean) - factorMeanY(0.5)}})
      .style("fill", function(d) {if (detailTargetFactors.includes(d.originalFactorInd)) {return "#555555"} else {return "#999999"}})
      .style("stroke-width", 4)
      .style("stroke", "none")
      .on('mouseover', mouseoverMean)
      .on('mousemove', mousemove)
      // .on('mouseleave', mouseleave)
      .on('mouseleave', function() {mouseleave.bind(this)("none")})
      // .on('click', mouseclick)
}
export default MeanRankBarChart