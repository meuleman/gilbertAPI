import * as d3 from 'd3'
// Bars showing the high scoring convolution bars

const ConvolutionBars = ({
  heatmapsvg,
  tooltip,
  colLabels,
  simThresh,
  data,
  simFactorXStart,
  simFactorXWidth,
  queryY,
  y,
  color,
  mousemove,
  mouseleave
} = {}) => {
  const dataFlat = (data.flatMap((region) => {
    const rank = region.rank
    const percentiles = region.percentiles

    return (
      percentiles.map((p, index) => {
        return ({
          percentile: p,
          factor: colLabels[index],
          rank: rank,
          factorInd: index,
        })
      })
    )
  }))
  const queryDataFlat = dataFlat.slice(0,34)

  const simFactorData = dataFlat.filter((d) => {
    return d.percentile * queryDataFlat[d.factorInd].percentile > simThresh
  }) // includes query region data

  // Factors similar to query
  const simFactorX = d3.scaleBand()
    .domain([...Array(colLabels.length).keys()])
    .range([simFactorXStart, simFactorXStart + simFactorXWidth])

  // mouseover for convolution bars
  const mouseoverConv = function (d, cellData) {
    const factor = cellData.factor
    tooltip
      .html(factor)
      .style('opacity', 1)
      .style('z-index', 1)
    d3.select(this)
      .style("stroke", "black")
      .style("stroke-width", "1")
  }

  let factorNum = 0
  let simRank = null
  heatmapsvg.selectAll()
    .data(simFactorData)
    .join("rect")
      .attr("x", function(d) {
        if (d.rank === simRank) {
          factorNum += 1
        } else {
          factorNum = 0
          simRank = d.rank
        }
        return simFactorX(factorNum)
      })
      .attr("y", function(d) {
        if (d.rank === 0) {
          return queryY(0.5) - (y.bandwidth() / 2)
        } else{
          return y(d.rank)
        }
      })
      .attr("width", simFactorX.bandwidth())
      .attr("height", y.bandwidth())
      .style("fill", function(d) {return color(d.factor)})
      .style("fill", function(d) {return color(d.factor)})
      .style("stroke-width", 0.1)
      .style("stroke", "black")
      .on('mouseover', mouseoverConv)
      .on("mousemove", mousemove)
      .on('mouseleave', function() {mouseleave.bind(this)(0.1)})
  
}
export default ConvolutionBars