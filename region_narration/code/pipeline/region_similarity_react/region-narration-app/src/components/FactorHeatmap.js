import React, { useEffect } from 'react'
import * as d3 from 'd3'
import { scaleLinear, timeFormatDefaultLocale } from 'd3'

function FactorHeatmap (props) {
    useEffect(() => {
        if (props.heatmapProps.data !== "") {

            const margin = { top: 75, right: 50, left: 50, bottom: 20 }
            const height = 700
            const width = 800
            const sigHeight = 60
            const barMax = 50
            const queryHeight = 60
            const colorHeight = 5
            const gap = 50

            const data = props.heatmapProps.data
            
            const queryData = data[0]

            const factorColors = queryData.factor_colors            
            factorColors[queryData.factors.indexOf('Quies')] = "#e0e0e0" // changing Quies color to something visable

            const dataFlat = (data.flatMap((region) => {
                const factors = region.factors
                const rank = region.rank
                const percentiles = region.percentiles
                const coords = region.coordinates
                const sentence = region.sentence
                const score = region.score
                const regionInd = region.region_ind
                return (
                    percentiles.map((p, index) => {
                        return ({
                            percentile: p,
                            factor: factors[index],
                            rank: rank,
                            coords: coords,
                            sentence: sentence,
                            score: score,
                            regionInd: regionInd,
                            factorInd: index
                        })
                    })
                )
            }))

            const queryDataFlat = dataFlat.slice(0,34)
            const heatmapData = dataFlat.slice(34)

            const factorSig = queryData.top_region_factor_Z_score
            const factorSigData = factorSig.map((value, index) => {
                if (value < 0) {
                    var barValue = 0
                } else if (value > barMax) {
                    var barValue = barMax
                } else {
                    var barValue = value
                }
                return ({
                    zscore: value,
                    sig: barValue,
                    factor: queryData.factors[index]
                })
            })

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

            const colLabels = queryData.factors

            // X axis
            const x = d3.scaleBand()
                .domain(colLabels)
                .range([0, width - margin.left - margin.right])
                .padding(0.05)

            heatmapsvg.append('g')
                .style('font', '10px times')
                .attr('class', 'x axis')
                .attr('transform', `translate(0, ${0})`)
                .call(d3.axisTop(x))
                .selectAll('text')
                .attr('transform', `rotate(-45)`)
                .style("text-anchor", "start")
                .attr('dx', '5')

            const topColorY = d3.scaleBand()
                .domain(["Factor Colors"])
                .range([0, colorHeight])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(topColorY).tickValues([]).tickSize(0))

            const queryY = d3.scaleLinear()
                .domain([1, 0])
                .range([colorHeight, queryHeight])
                // .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(queryY).tickValues([0, 0.5, 1]))

            

            const botColorY = d3.scaleBand()
                .domain(["Factor Colors"])
                .range([queryHeight + gap - colorHeight, queryHeight + gap])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(botColorY).tickValues([]).tickSize(0))

            // Y axis
            const y = d3.scaleBand()
                .domain(heatmapData.map(d => d.rank))
                .range([queryHeight + gap, height - margin.top - margin.bottom - sigHeight])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(y).tickValues([10, 20, 30, 40, 50, 60, 70, 80, 90]))

            

            const sigY = d3.scaleLinear()
                .domain([barMax, 0])
                .range([height - margin.top - margin.bottom - sigHeight, height - margin.top - margin.bottom])
                // .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '6px times')
                .call(d3.axisLeft(sigY).tickValues([0, 10, 20, 30, 40, 50]))
            
            

            // console.log(queryY(queryData.percentiles))
            
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -(height) / 2 + 15)
                .attr("y", -30)
                .attr("transform", "rotate(-90)")
                .style('font', '12px times')
                .text("Similar Regions")
            
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -35)
                .attr("y", -30)
                .attr("transform", "rotate(-90)")
                .style('font', '12px times')
                .text("Query Region")

            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -(height - margin.bottom) + 100)
                .attr("y", -30)
                .attr("transform", "rotate(-90)")
                .style('font', '12px times')
                .text("Factor Z-scores")
            
            // Build color scale
            const color = d3.scaleOrdinal()
                .range(factorColors)
                .domain(colLabels)

            const heatmapColor = d3.scaleLinear()
                .range(["cyan", "white", "yellow"])
                .domain([0, 0.5, 1])

            

            const opacity = d3.scaleLinear()
                .range([0, 1])
                .domain([0.5, 1])

            // const sigOpacity = d3.scaleLinear()
            //     .range([0, 1])
            //     .domain([0, Math.max(...factorSig)])



            // console.log(factorSigData.map((d) => {
            //     return(
            //         sigOpacity(d.sig)
            //     )
            // }
            // ))

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
                    const regionInd = cellData.regionInd
                    const regionChr = coords.split(':')[0]
                    const regionStart = coords.split(':')[1].split('-')[0]
                    const regionEnd = coords.split(':')[1].split('-')[1]
                    const queryCoords = queryData.coordinates
                    const queryChr = queryCoords.split(':')[0]
                    const queryStart = queryCoords.split(':')[1].split('-')[0]
                    const queryEnd = queryCoords.split(':')[1].split('-')[1]
                    const sentence = cellData.sentence
                    const score = Math.round(cellData.score * 10000) / 10000
                    var searchFromSim = () => props.heatmapProps.search(regionInd)
        
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
                        "<br><button id='tooltipButton'>Search From This Region</button>")
                        .style('opacity', 1)
                }
                document.getElementById("tooltipButton").addEventListener("click", searchFromSim)
                
            }

            const mouseoverSig = function (d, cellData) {
                if (!(clicked)) {
                    const Zscore = Math.round(cellData.zscore * 10000) / 10000
                    const factor = cellData.factor
                    tooltip
                        .html("Factor: " + factor + "<br>Factor Z-score: " + Zscore)
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
        
            const mouseclick = function () {
                clicked = !clicked
            }


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
                    // .style("opacity", function(d) {return opacity(d.percentile)})

                    .on('mouseover', mouseover)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)

            
            // console.log(queryDataFlat)
            heatmapsvg.selectAll()
                .data(factorSigData)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", function(d) {return sigY(d.sig)})
                    .attr("width", x.bandwidth())
                    .attr("height", function(d) {return height - margin.top - margin.bottom - sigY(d.sig)})
                    .style("fill", "#555555")
                    .style("stroke-width", 4)
                    .style("stroke", "none")
                    // .style("opacity", function(d) {return sigOpacity(d.sig)})

                    .on('mouseover', mouseoverSig)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)


            
            // console.log(props.heatmapProps.detailLevel)
            // console.log(queryDataFlat.map(
            //     function(d) {return color(d.factor)}
            //     // return d.factorInd
            // ))
            heatmapsvg.selectAll()
                .data(queryDataFlat)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", topColorY("Factor Colors"))
                    .attr("width", x.bandwidth())
                    .attr("height", colorHeight)
                    .style("fill", function(d) {return color(d.factor)})
                    .style("stroke-width", 4)
                    .style("stroke", "none")

            heatmapsvg.selectAll()
                .data(queryDataFlat)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", botColorY("Factor Colors"))
                    .attr("width", x.bandwidth())
                    .attr("height", colorHeight)
                    .style("fill", function(d) {return color(d.factor)})
                    .style("stroke-width", 4)
                    .style("stroke", "none")

            heatmapsvg.selectAll()
                .data(queryDataFlat)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", function(d) {return queryY(d.percentile)})
                    .attr("width", x.bandwidth())
                    .attr("height", function(d) {return queryHeight - queryY(d.percentile)})
                    .style("fill", function(d) {if (d.factorInd < props.heatmapProps.detailLevel) {return "#555555"} else {return "#999999"}})
                    .style("stroke-width", 4)
                    .style("stroke", "none")
                    // .style("opacity", function(d) {return sigOpacity(d.sig)})

                    // .on('mouseover', mouseoverSig)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)
            
            heatmapsvg
                .append("line")
                    .style("stroke", "black")
                    .style("stroke-width", 0.5)
                    .attr("x1", 0)
                    .attr("y1", queryY(0.5))
                    .attr("x2", width - margin.left - margin.right)
                    .attr("y2", queryY(0.5))
        }
    })
    return (
        <div id='factorHeatmap' />
    )
};

export default FactorHeatmap
