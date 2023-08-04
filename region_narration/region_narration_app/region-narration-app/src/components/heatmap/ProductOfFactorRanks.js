// Product of Factor Ranks

const ProductOfFactorRanks = ({
  heatmapsvg,
  data,
  prodFactorRankX,
  queryY,
  y,
  mouseoverLabel,
  mousemoveLabel,
  mouseleaveLabel,
} = {}) => {
  const prodFactorRanks = data.map((d) => {
    return {
      rank: d.rank,
      prodFactorRank: d.prod_rank
    }
  })

  // Product of factor ranks
  heatmapsvg.selectAll()
    .data(prodFactorRanks)
      .join("text")
        .attr("x", prodFactorRankX)
        .attr("y", function(d) {
          if (d.rank === 0) {
            return queryY(0.5) + 4
          } else {
            return y(d.rank) + y.bandwidth() - 5
          }
            
        })
        .attr("text-anchor", "start")
        .style('font', '12px sans-serif')
        .text(function(d) {
          return String(
            Math.round(d.prodFactorRank * 100) / 100
          ) + "%"
        })

  // Product of Factor Ranks Label
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("y", queryY(1))
    .attr("x", prodFactorRankX + 18)
    .style('font', '13px sans-serif')
    .text("Factor Combination")
    .on('mouseover', () => mouseoverLabel("Regional ranks for factors used in similarity search"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)
  heatmapsvg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("y", queryY(0.75))
    .attr("x", prodFactorRankX + 18)
    .style('font', '13px sans-serif')
    .text("Ranks")
    .on('mouseover', () => mouseoverLabel("Regional ranks for factors used in similarity search"))
    .on('mousemove', mousemoveLabel)
    .on('mouseleave', mouseleaveLabel)
  
}
export default ProductOfFactorRanks