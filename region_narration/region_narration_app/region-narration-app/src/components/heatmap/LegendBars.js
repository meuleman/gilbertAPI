// bars for heatmap legend

const LegendBars = ({
  heatmapsvg,
  legendStartX,
  legendEndX,
  data,
  detailTargetFactors,
  legendY,
  mouseoverLabel,
  mousemoveLabel,
  mouseleaveLabel,
  searchedByFactor,
} = {}) => {
  const dataFlat = (data.flatMap((region) => {
    const percentiles = region.percentiles
    return (
      percentiles.map((p, index) => {
        return ({
          percentile: p,
          originalFactorInd: data[0].factor_order[index]
        })
      })
    )
  }))

  let heatmapData = null
  if (searchedByFactor) {
    heatmapData = dataFlat
  } else {
    heatmapData = dataFlat.slice(34)
  }

  let allSimRanks = heatmapData.map((d) => {
    if (detailTargetFactors.includes(d.originalFactorInd)) {
      return d.percentile
    } else {
      return null
    }
  })
  allSimRanks = allSimRanks.filter((value) => {
    return value != null
  })
  const maxRank = Math.max(...allSimRanks)
  const minRank = Math.min(...allSimRanks)
  const meanRank = allSimRanks.reduce((rankSum, val) => rankSum + val, 0) / allSimRanks.length


  // add heatmap colorbar max rank
  heatmapsvg.append("line")
    .style("stroke", "green")
    .style("stroke-width", 5)
    .attr("x1", legendStartX)
    .attr("y1", legendY(maxRank))
    .attr("x2", legendEndX)
    .attr("y2", legendY(maxRank))
    .on('mouseover', () => mouseoverLabel("Max Rank in Similar Regions"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)

  // add heatmap colorbar mean rank
  heatmapsvg.append("line")
    .style("stroke", "orange")
    .style("stroke-width", 5)
    .attr("x1", legendStartX)
    .attr("y1", legendY(meanRank))
    .attr("x2", legendEndX)
    .attr("y2", legendY(meanRank))
    .on('mouseover', () => mouseoverLabel("Mean Rank in Similar Regions"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)

  // add heatmap colorbar min rank
  heatmapsvg.append("line")
    .style("stroke", "red")
    .style("stroke-width", 5)
    .attr("x1", legendStartX)
    .attr("y1", legendY(minRank))
    .attr("x2", legendEndX)
    .attr("y2", legendY(minRank))
    .on('mouseover', () => mouseoverLabel("Min Rank in Similar Regions"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)
}
export default LegendBars