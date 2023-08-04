// Mean Squared Error report

const MeanSquaredError = ({
  heatmapsvg,
  colLabels,
  data,
  detailTargetFactors,
  MSELabelX,
  MSELabelY,
} = {}) => {
  const queryData = data[0]

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


  let heatmapData = dataFlat.slice(34)

  // const factorMean
  let factorSums = new Array(colLabels.length).fill(0)
  heatmapData.forEach((d) => {
    factorSums[d.factorInd] += d.percentile
  })
  const factorMeans = factorSums.map((value, index) => {
    return {
      mean: value / heatmapData[heatmapData.length - 1].rank,
      factor: colLabels[index],
      originalFactorInd: data[0].factor_order[index]
    }
  })

  let errSum = 0
  factorMeans.forEach((d, index) => {
    const queryOFI = queryData.factor_order[index]
    const simOFI = d.originalFactorInd
    if (detailTargetFactors.includes(queryOFI) && detailTargetFactors.includes(simOFI)) {
      errSum += Math.pow(queryData.percentiles[index] - d.mean, 2)
    }
  })
  const MSE = errSum / detailTargetFactors.length

  // add MSE
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "start")
    .attr("y", MSELabelY)
    .attr("x", MSELabelX)
    .style('font', '13px sans-serif')
    .text("MSE=" + String(MSE.toExponential(3)))
  
}
export default MeanSquaredError