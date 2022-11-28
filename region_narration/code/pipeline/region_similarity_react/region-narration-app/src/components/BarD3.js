import React from 'react'
import * as d3 from 'd3'

class BarD3 extends React.Component {
  render () {
    const data = this.props.barProps.regionalVals
    // console.log(data)

    if (data.length !== 0) {
      const margin = { top: 0, right: 50, left: 50, bottom: 10 }
      const height = 1000
      const width = 800

      d3.select('svg').remove()

      if (document.getElementById('barChart').getElementsByTagName('*').length !== 0) {
        document.getElementById('barChart').removeChild(document.getElementById('barChart').lastElementChild)
      }
      //
      // d3.select("div").remove()
      const barContainer = d3.select('#barChart')

      const barsvg = barContainer
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .append('g')
        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')

      const metrics = Object.keys(data[0]).slice(0, 3)
      // var regions = d3.map(data, function(d){return(d.region)}).keys()
      const regionSums = data.map(x => {
        return Object.values(x).slice(0, 3).reduce(
          (regionTotal, metricVal) => regionTotal + metricVal, 0
        )
      })

      const regionScores = data.map(x => {
        return 1 - Object.values(x).slice(0, 3).reduce(
          (scoreIterator, metricVal) => scoreIterator * (1 - metricVal), 1
        )
      })

      const scaledData = data.map((x, r) => {
        const regionScaled = {}
        metrics.map((m) => {
          const s = regionScaled[m] = x[m] / regionSums[r] * regionScores[r]
          return s
        })
        regionScaled.region = r
        return regionScaled
      })

      // X axis
      const x = d3.scaleLinear()
        .domain([0, Math.max(...regionSums)])
        .range([0, width - margin.left - margin.right])
      barsvg.append('g')
        .style('font', '8px times')
        .attr('transform', `translate(0, ${height})`)
        .call(d3.axisBottom(x))

      // Y axis
      const y = d3.scaleBand()
        .domain(data.map(d => d.region))
        .range([0, height - margin.top - margin.bottom])
        .padding(0.1)
      barsvg.append('g')
        .style('font', '8px times')
      // .attr("transform", "translate(0," + width + ")")
        .call(d3.axisLeft(y))

      const color = d3.scaleOrdinal()
        .domain(metrics)
        .range(['#a6cee3', '#1f78b4', '#b2df8a'])

      const stackedData = d3.stack().keys(metrics)(scaledData)

      const tooltip = barContainer
        .append('div')
        .style('opacity', 0)
        .attr('class', 'tooltip')
        .style('background-color', 'white')
        .style('border', 'solid')
        .style('border-width', '1px')
        .style('border-radius', '5px')
        .style('padding', '10px')
        .style('position', 'absolute')
        .style('display', 'inline')
        .style('z-index', 1)

      const mouseover = function (d, i) {
        const metricName = d3.select(this.parentNode).datum().key
        const regionNum = i.data.region
        tooltip
          .html('Region: ' + regionNum + '<br>Metric: ' + metricName)
          .style('opacity', 1)
      }

      const mousemove = function (e) {
        tooltip.style('left', (e.pageX + 20) + 'px')
        tooltip.style('top', (e.pageY) + 'px')
        // console.log(tooltip.style.top, tooltip.style.left)
      }
      const mouseleave = function (d) {
        tooltip
          .style('opacity', 0)
      }

      barsvg.append('g')
        .selectAll('g')
        .data(stackedData)
        .enter().append('g')
        .attr('fill', function (d) { return color(d.key) })
        .selectAll('rect')
        .data(function (d) {
          return d
        })
        .enter().append('rect')
        .attr('x', function (d) {
          return x(d[0])
        })
        .attr('y', function (d) {
          // console.log(y(d.data.region))
          return y(d.data.region)
        })
        .attr('width', function (d) {
          // console.log(x(d[1]) - x(d[0]))
          return x(d[1]) - x(d[0])
        })
        .attr('height', y.bandwidth)
        .on('mouseover', mouseover)
        .on('mousemove', mousemove)
        .on('mouseleave', mouseleave)
    }

    return (

      <div id='barChart'>

        {/* <div>
                    {document.getElementById('barChart')}
                </div> */}
      </div>

    )
  }
};

export default BarD3
