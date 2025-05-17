import * as d3 from 'd3'
import factorsArr from '../../factors.json'

// Heatmap
const Heatmap = ({
  heatmapsvg,
  tooltip,
  colLabels,
  plotStartX,
  plotEndX,
  figureStartY,
  searchedByFactor,
  data,
  heatmapStartY,
  heatmapEndY,
  mouseoverLabel,
  mousemoveLabel,
  mouseleaveLabel,
  heatmapLabelY,
  heatmapLabelX,
  mousemove,
  mouseleave,
  detailTargetFactors,
} = {}) => {
  const originalLabelOrder = factorsArr.map((factor) => {
    return factor.name
  })

  const dataFlat = (data.flatMap((region) => {
    const rank = region.rank
    const percentiles = region.percentiles
    const score = region.score
    return (
      percentiles.map((p, index) => {
        return ({
          percentile: p,
          factor: colLabels[index],
          rank: rank,
          score: score,
        })
      })
    )
  }))

  let heatmapData = ""
  let simData = ""
  let heatmapLabel = ""
  if (searchedByFactor) {
    heatmapData = dataFlat
    simData = data
    heatmapLabel = "Top Regions"
  } else {
    heatmapData = dataFlat.slice(34)
    simData = data.slice(1)
    heatmapLabel = "Similar Regions"
  }

  // mouseover for heatmap 
  const mouseover = function (d, cellData) {
    const factor = cellData.factor
    const percentile = Math.round((cellData.percentile * 100) * 10) / 10
    const score = Math.round(cellData.score * 10000) / 10000
    tooltip
      .html("Similarity Score: " + score +
      "<br>" + factor + ": " + percentile + "%")
      .style('opacity', 1)
      .style('z-index', 1)
    d3.select(this)
      .style("stroke", "black")
      .style("stroke-width", "1")
  }
  
  // X axis
  const x = d3.scaleBand()
    .domain(colLabels)
    .range([plotStartX, plotEndX])
    .padding(0.05)
  heatmapsvg.append('g')
    .style('font', '10px sans-serif')
    .attr('class', 'x axis')
    .attr('transform', `translate(0, ` + figureStartY + `)`)
    .call(d3.axisTop(x))
    .selectAll('text')
    .attr('transform', `rotate(-45)`)
    .style("text-anchor", "start")
    .style("font-weight", (factor) => {
      const factorInd = originalLabelOrder.indexOf(factor)
      if (detailTargetFactors.includes(factorInd)) {
        return "900"
      } else {
        return "normal"
      }
    })
    .attr('dx', '5')

  // heatmap y axis
  const y = d3.scaleBand()
    .domain(heatmapData.map(d => d.rank))
    .range([heatmapStartY, heatmapEndY])
    .padding(0.1)
  heatmapsvg.append('g')
    .style('font', '10px sans-serif')
    .attr("transform", "translate(" + plotStartX + ", 0)")
    .call(d3.axisLeft(y).tickFormat((regionInd) => {
      if (!searchedByFactor) {
        regionInd -= 1
      }
      const coords = simData[regionInd].coordinates
      return coords //+ ' - ' + regionInd
    }))

  heatmapsvg.append('g')
    .style('font', '10px sans-serif')
    .attr("transform", "translate(-60, 0)")
    .call(d3.axisLeft(y).tickSize(0))
    .call((g) => {return g.select(".domain").remove()})

  // heatmap color scale
  const heatmapColor = d3.scaleLinear()
    .range(["#FFC20A", "#0C7BDC"])
    .domain([1, 0])

  // add similar regions label
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "end")
    .attr("y", heatmapLabelY)
    .attr("x", heatmapLabelX)
    // .attr("transform", "rotate(-90)")
    .style('font', '13px sans-serif')
    .text(heatmapLabel)
    .on('mouseover', () => mouseoverLabel("Ranks for factors in similar regions"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)

  // add factor heatmap
  heatmapsvg.selectAll()
    .data(heatmapData)
    .join("rect")
      .attr("x", function(d) {return x(d.factor)})
      .attr("y", function(d) {return y(d.rank)})
      .attr("width", x.bandwidth())
      .attr("height", y.bandwidth())
      .style("fill", function(d) {return heatmapColor(d.percentile)})
      .style("stroke-width", 4)
      .style("stroke", "none")
      .on('mouseover', mouseover)
      .on('mousemove', mousemove)
      // .on('mouseleave', mouseleave)
      .on('mouseleave', function() {mouseleave.bind(this)("none")})
  
  return {
    x: x,
    y: y,
    heatmapColor: heatmapColor,
  }

}
export default Heatmap