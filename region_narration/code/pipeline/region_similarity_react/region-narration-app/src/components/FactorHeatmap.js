import React, { useEffect } from 'react'
import * as d3 from 'd3'
import { scaleLinear } from 'd3'

function FactorHeatmap (props) {
    useEffect(() => {
        if (props.heatmapProps.data !== "") {
            const data = props.heatmapProps.data
            
            const factorColors = data[0].factor_colors
            factorColors[33] = "#e0e0e0" // changing Quies color to something visable
        
            const dataFlat = (data.flatMap((region) => {
                const factors = region.factors
                const rank = region.rank
                const percentiles = region.percentiles
                const coords = region.coordinates
                const sentence = region.sentence
                const score = region.score
                return (
                    percentiles.map((p, index) => {
                        return ({
                            percentile: p,
                            factor: factors[index],
                            rank: rank,
                            coords: coords,
                            sentence: sentence,
                            score: score
                        })
                    })
                )
            }))

            const margin = { top: 75, right: 50, left: 50, bottom: 10 }
            const height = 600
            const width = 800

            if (document.getElementById('factorHeatmap').getElementsByTagName('*').length !== 0) {
                d3.select('svg').remove()
                document.getElementById('factorHeatmap').removeChild(document.getElementById('factorHeatmap').lastElementChild)
            }
        
            const heatmapContainer = d3.select('#factorHeatmap')
            const heatmapsvg = heatmapContainer
                .append('svg')
                .attr('width', width)
                .attr('height', height)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')

            // const rowLabels = [data.map((d) => {
            //     return d.rank
            // })]

            const colLabels = data[0].factors

            // X axis
            const x = d3.scaleBand()
                .domain(colLabels)
                .range([0, width - margin.left - margin.right])
                .padding(0.05)

            heatmapsvg.append('g')
                .style('font', '10px times')
                .attr('class', 'x axis')
                // .attr('transform', `translate(0, ${height})`)
                .call(d3.axisTop(x))
                .selectAll('text')
                .attr('transform', `rotate(-45)`)
                .style("text-anchor", "start")
                .attr('dx', '5')

            // Y axis
            const y = d3.scaleBand()
                .domain(data.map(d => d.rank))
                .range([0, height - margin.top - margin.bottom])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(y))
            

            // Build color scale
            const color = d3.scaleOrdinal()
                .range(factorColors)
                .domain(colLabels)

            const opacity = d3.scaleLinear()
                .range([0, 1])
                .domain([0.5, 1])

            // tool tip and functions
            const tooltip = heatmapContainer
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
                .style('max-width', '500px')

            let clicked = false

            

            const mouseover = function (d, cellData) {
                if (!(clicked)) {
                    const factor = cellData.factor
                    const rank = cellData.rank
                    const percentile = Math.round(cellData.percentile * 10000) / 10000
                    const coords = cellData.coords
                    const regionChr = coords.split(':')[0]
                    const regionStart = coords.split(':')[1].split('-')[0]
                    const regionEnd = coords.split(':')[1].split('-')[1]
                    const queryCoords = data[0].coordinates
                    const queryChr = queryCoords.split(':')[0]
                    const queryStart = queryCoords.split(':')[1].split('-')[0]
                    const queryEnd = queryCoords.split(':')[1].split('-')[1]
                    const sentence = cellData.sentence
                    const score = Math.round(cellData.score * 10000) / 10000
                    var searchFromSim = () => props.heatmapProps.search(coords)
        
                    const indexPath = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&dtt=200000&dtc=av3&gatt=cv&gac=off'
                    const epilogosPath = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&gatt=cv'
                    const queryIndexPath = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=' + queryChr + '&chrRight=' + queryChr + '&start=' + queryStart + '&stop=' + queryEnd + '&dtt=200000&dtc=av3&gatt=cv&gac=off'
                    const queryEpilogosPath = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=' + queryChr + '&chrRight=' + queryChr + '&start=' + queryStart + '&stop=' + queryEnd + '&gatt=cv'
                    tooltip
                        .html("Coordinates: " + coords + "<br>Region Rank: " + rank + "<br>Similarity Score: " + score +
                        "<br>Factor: " + factor + "<br>Factor Percentile: " + percentile + "<br>" + sentence +
                        "<br><a href=" + indexPath + ">Index Browser</a>" +
                        "<br><a href=" + epilogosPath + ">Epilogos Browser</a>" +
                        "<br><a href=" + queryIndexPath + ">Query: Index Browser</a>" +
                        "<br><a href=" + queryEpilogosPath + ">Query: Epilogos Browser</a>" +
                        // "<br><button onclick='s(coords)'>Search This Region</button>")
                        "<br><button id='tooltipButton'>Search From This Region</button>")
                        .style('opacity', 1)
                }
                document.getElementById("tooltipButton").addEventListener("click", searchFromSim)
                // const element = document.getElementById("tooltipButton")
                // element.addEventListener("click", console.log(coords));
                
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
        
            const mouseclick = function () {
                clicked = !clicked
            }


            heatmapsvg.selectAll()
                .data(dataFlat)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", function(d) {return y(d.rank)})
                    .attr("width", x.bandwidth())
                    .attr("height", y.bandwidth())
                    .style("fill", function(d) {return color(d.factor)})
                    .style("stroke-width", 4)
                    .style("stroke", "none")
                    .style("opacity", function(d) {return opacity(d.percentile)})

                    .on('mouseover', mouseover)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)

                
        }
    })
    return (
        <div id='factorHeatmap' />
    )
};

export default FactorHeatmap
