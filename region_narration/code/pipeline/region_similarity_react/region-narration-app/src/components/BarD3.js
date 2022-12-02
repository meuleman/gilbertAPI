import React, { useEffect } from 'react'
import * as d3 from 'd3'

function BarD3 (props) {
  useEffect(() => {
    if (props.barProps.regionalVals !== '') {
      const dataAll = props.barProps.regionalVals
      const data = dataAll.slice(1)

      const margin = { top: 0, right: 50, left: 50, bottom: 10 }
      const height = 1000
      const width = 800

      if (document.getElementById('barChart').getElementsByTagName('*').length !== 0) {
        d3.select('svg').remove()
        document.getElementById('barChart').removeChild(document.getElementById('barChart').lastElementChild)
      }

      const barContainer = d3.select('#barChart')
      const barsvg = barContainer
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .append('g')
        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')

      const metrics = Object.keys(data[0]).slice(0, 3)
      // var regions = d3.map(data, function(d){return(d.rank)}).keys()
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
        // console.log(r)
        regionScaled.rank = (r + 1)
        return regionScaled
      })

      // X axis
      const x = d3.scaleLinear()
        .domain([0, 1])
        .range([0, width - margin.left - margin.right])

      barsvg.append('g')
        .style('font', '8px times')
        .attr('transform', `translate(0, ${height})`)
        .call(d3.axisBottom(x))

      // Y axis
      const y = d3.scaleBand()
        .domain(data.map(d => d.rank))
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

      // tool tip and functions
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

      let clicked = false

      const mouseover = function (d, i) {
        if (!(clicked)) {
          const metricStacked = d3.select(this.parentNode).datum()
          const metricName = metricStacked.key
          const rankNum = i.data.rank
          const coords = data[(rankNum - 1)].coordinates
          const regionChr = coords.split(':')[0]
          const regionStart = coords.split(':')[1].split('-')[0]
          const regionEnd = coords.split(':')[1].split('-')[1]
          const regionSentence = data[(rankNum - 1)].sentence
          const regionScore = regionScores[(rankNum - 1)]
          const indexPath = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&dtt=200000&dtc=av3&gatt=cv&gac=off'
          const epilogosPath = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&gatt=cv'
          tooltip
            .html('Coordinates: ' + coords + '<br>Region Rank: ' + rankNum + '<br>Fusion Score: ' + regionScore + '<br>Metric: ' +
              metricName + '<br>Summary: ' + regionSentence + '<br><a href=' + indexPath + '>' + 'Index Browser' + '</a>' + 
              '<br><a href=' + epilogosPath + '>' + 'Epilogos Browser' + '</a>')
            .style('opacity', 1)
        }
      }

      const mousemove = function (e) {
        if (!(clicked)) {
          tooltip.style('left', (e.pageX + 20) + 'px')
          tooltip.style('top', (e.pageY) + 'px')
        };
      }

      const mouseleave = function (d) {
        if (!(clicked)) {
          tooltip
            .style('opacity', 0)
        }
      }

      const mouseclick = function (d) {
        if (clicked) {
          clicked = false
        } else {
          clicked = true
        };
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
          // console.log((d.data.rank))
          return y(d.data.rank)
        })
        .attr('width', function (d) {
          // console.log(x(d[1]) - x(d[0]))
          return x(d[1]) - x(d[0])
        })
        .attr('height', y.bandwidth)
        .on('mouseover', mouseover)
        .on('mousemove', mousemove)
        .on('mouseleave', mouseleave)
        .on('click', mouseclick)
    }
  })
  return (
    <div id='barChart' />
  )
};

export default BarD3
