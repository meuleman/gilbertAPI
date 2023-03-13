import React, { useEffect, useState } from 'react'
import * as d3 from 'd3'
import { scaleLinear } from 'd3'
// import { scaleLinear, timeFormatDefaultLocale } from 'd3'

function FactorHeatmap (props) {
    useEffect(() => {
        if (props.heatmapProps.data !== "") {

            const margin = { top: 75, right: 120, left: 25, bottom: 220 }
            const height = 900
            const width = 800
            const sigHeight = 60
            const zMax = 50
            const zMin = -10
            const queryHeight = 60
            const colorHeight = 5
            const gap = 50


            // data
            const data = props.heatmapProps.data
            
            const queryData = data[0]
            const colLabels = queryData.factors
            const roiURL = queryData.roiURL
            const detailTargetFactors = queryData.target_factors.slice(0, props.heatmapProps.detailLevel)

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
                            factorInd: index,
                            originalFactorInd: queryData.query_factor_order[index]
                        })
                    })
                )
            }))

            const queryDataFlat = dataFlat.slice(0,34)
            const heatmapData = dataFlat.slice(34)

            const allSimRanks = heatmapData.map((d) => {
                return d.percentile
            })
            const maxRank = Math.max(...allSimRanks)
            const minRank = Math.min(...allSimRanks)
            const meanRank = allSimRanks.reduce((rankSum, val) => rankSum + val, 0) / allSimRanks.length

            const factorSig = queryData.top_region_factor_Z_score
            
            const factorSigData = factorSig.map((value, index) => {
                var barValue = 0
                if (value < zMin) {
                    barValue = zMin
                } else if (value > zMax) {
                    barValue = zMax
                } else {
                    barValue = value
                }
                return ({
                    zscore: value,
                    sig: barValue,
                    factor: queryData.factors[index],
                    originalFactorInd: queryData.query_factor_order[index]
                })
            })

            // clear document
            if (document.getElementById('factorHeatmap').getElementsByTagName('*').length !== 0) {
                d3.select('svg').remove()
                document.getElementById('factorHeatmap').removeChild(document.getElementById('factorHeatmap').lastElementChild)
            }
        
            // container
            const heatmapContainer = d3.select('#factorHeatmap')
            const heatmapsvg = heatmapContainer
                .append('svg')
                .attr('width', width)
                .attr('height', height)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')

            // X axis
            const x = d3.scaleBand()
                .domain(colLabels)
                .range([margin.left, width - margin.right])
                .padding(0.05)
            heatmapsvg.append('g')
                .style('font', '10px sans-serif')
                .attr('class', 'x axis')
                .attr('transform', `translate(0, ${0})`)
                .call(d3.axisTop(x))
                .selectAll('text')
                .attr('transform', `rotate(-45)`)
                .style("text-anchor", "start")
                .attr('dx', '5')

            // colors for each factor
            const topColorY = d3.scaleBand()
                .domain(["Factor Colors"])
                .range([0, colorHeight])
                .padding(0.1)
            heatmapsvg.append('g')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(topColorY).tickValues([]).tickSize(0))

            // query region bar chart y axis
            const queryY = d3.scaleLinear()
                .domain([1, 0])
                .range([colorHeight, queryHeight])
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(queryY).tickValues([1.0, 0.5, 0]).tickFormat(d3.format(".0%")))

            // colors for each factor
            const botColorY = d3.scaleBand()
                .domain(["Factor Colors"])
                .range([queryHeight + gap - colorHeight, queryHeight + gap])
                .padding(0.1)
            heatmapsvg.append('g')

            // heatmap y axis
            const y = d3.scaleBand()
                .domain(heatmapData.map(d => d.rank))
                .range([queryHeight + gap, height - margin.bottom - sigHeight])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(y).tickValues([10, 20, 30, 40, 50, 60, 70, 80, 90]))

            // z score y axis
            const sigY = d3.scaleLinear()
                .domain([zMax, zMin])
                .range([height - margin.bottom - sigHeight, height - margin.bottom])
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(sigY).tickValues([-10, 0, 10, 20, 30, 40, 50]))

            // colors for each factor
            const color = d3.scaleOrdinal()
                .range(factorColors)
                .domain(colLabels)

            // heatmap color scale
            const heatmapColor = d3.scaleLinear()
                .range(["#FFC20A", "#0C7BDC"])
                .domain([1, 0])
            
            // colorbar for heatmap legend
            const heatmapLegend = heatmapsvg.append("linearGradient")
                .attr("id", "heatmap-legend")
                .attr("x1", "0%")
                .attr("y1", "100%")
                .attr("x2", "0%")
                .attr("y2", "0%")
            heatmapLegend
                .append("stop")
                    .attr("offset", "0%")
                    .attr("stop-color", heatmapColor(0))
            heatmapLegend
                .append("stop")
                    .attr("offset", "100%")
                    .attr("stop-color", heatmapColor(1))

            // heatmap legend y axis
            const legendY = d3.scaleLinear()
                .domain([1, 0])
                .range([queryHeight + gap, height - margin.bottom - sigHeight])
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + (width - margin.right + 30) + ", 0)")
                .call(d3.axisRight(legendY).tickValues([0, 0.25, 0.50, 0.75, 1.0]).tickFormat(d3.format(".0%")))

            // tool tip
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

            // click tracker
            let clicked = false
            const mouseclick = function () {
                clicked = !clicked
                if (!clicked) {
                    heatmapsvg.selectAll("rect")
                        .style("stroke", "none")
                }
            }
            // const [clicked, setClicked] = useState(
            //     false
            // );
            // const mouseclick = function () {
            //     setClicked(!clicked)
            // }

            // mouseover for heatmap 
            const mouseover = function (d, cellData) {
                if (!(clicked)) {
                    const factor = cellData.factor
                    const rank = cellData.rank
                    const percentile = Math.round((cellData.percentile * 100) * 10) / 10
                    const coords = cellData.coords
                    const regionInd = cellData.regionInd
                    const regionChr = coords.split(':')[0]
                    const regionStart = coords.split(':')[1].split('-')[0]
                    const regionEnd = coords.split(':')[1].split('-')[1]
                    const sentence = cellData.sentence
                    const score = Math.round(cellData.score * 10000) / 10000
                    var searchFromSim = () => props.heatmapProps.search(regionInd)
                    const indexPath = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&dtt=200000&dtc=av3&gatt=cv&gac=off'
                    const epilogosPath = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&gatt=cv'
                    if (roiURL != "roiURL not available") {
                        var roiURLTag = "<br><a href=" + roiURL + " target=_blank rel=noopener noreferrer>roiURL</a>"
                    } else {
                        var roiURLTag = ""
                    }
                    tooltip
                        .html("Region: " + rank + "<br>Coordinates: " + coords + "<br>Similarity Score: " + score +
                        "<br>Factor: " + factor + "<br>Factor Enrichment: " + percentile + "%<br>" + sentence +
                        "<br><a href=" + indexPath + " target=_blank rel=noopener noreferrer>Index Browser</a>" +
                        "<br><a href=" + epilogosPath + " target=_blank rel=noopener noreferrer>Epilogos Browser</a>" +
                        roiURLTag +
                        "<br><button id='tooltipButton'>Search From This Region</button>")
                        .style('opacity', 1)
                    document.getElementById("tooltipButton").addEventListener("click", searchFromSim)
                    d3.select(this)
                        .style("stroke", "black")
                        .style("stroke-width", "1")
                }
            }

            // mouseover for z score bar plot
            const mouseoverSig = function (d, cellData) {
                if (!(clicked)) {
                    const Zscore = Math.round(cellData.zscore * 10000) / 10000
                    const factor = cellData.factor
                    tooltip
                        .html("Factor: " + factor + "<br>Factor Z-score: " + Zscore)
                        .style('opacity', 1)
                    d3.select(this)
                        .style("stroke", "black")
                        .style("stroke-width", "1")
                }
            }

            // mouseover for query region bar chart
            const mouseoverQuery = function (d, cellData) {
                if (!(clicked)) {
                    const factor = cellData.factor
                    // const rank = cellData.rank
                    const percentile = Math.round((cellData.percentile * 100) * 10) / 10
                    const coords = cellData.coords
                    const regionChr = coords.split(':')[0]
                    const regionStart = coords.split(':')[1].split('-')[0]
                    const regionEnd = coords.split(':')[1].split('-')[1]
                    const sentence = cellData.sentence
        
                    const indexPath = 'https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&dtt=200000&dtc=av3&gatt=cv&gac=off'
                    const epilogosPath = 'https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=' + regionChr + '&chrRight=' + regionChr + '&start=' + regionStart + '&stop=' + regionEnd + '&gatt=cv'
                    tooltip
                        .html("Coordinates: " + coords +
                        "<br>Factor: " + factor + "<br>Factor Enrichment: " + percentile + "%<br>" + //sentence +
                        "<br><a href=" + indexPath + " target=_blank rel=noopener noreferrer>Index Browser</a>" +
                        "<br><a href=" + epilogosPath + " target=_blank rel=noopener noreferrer>Epilogos Browser</a>")
                        .style('opacity', 1)
                    d3.select(this)
                        .style("stroke", "black")
                        .style("stroke-width", "1")
                }
            }

            // mousemove
            const mousemove = function (e) {
                if (!(clicked)) {
                    tooltip.style('left', (e.pageX + 20) + 'px')
                    tooltip.style('top', (e.pageY) + 'px')
                };
            }

            // mouseleave
            const mouseleave = function (d) {
                if (!(clicked)) {
                    tooltip
                        .style('opacity', 0)
                    d3.select(this)
                        .style("stroke", "none")
                }
            }

            const labelTooltip = heatmapContainer
                .append('div')
                .style('opacity', 0)
                .attr('class', 'tooltip')
                .style('background-color', 'white')
                .style('border', 'solid')
                .style('border-width', '1px')
                .style('padding', '2px')
                .style('position', 'absolute')
                .style('display', 'inline')
                .style('z-index', 1)
                .style('max-width', '500px')


            // mouseover for query region bar chart label
            const mouseoverQueryLabel = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Ranks for factors in the query region</div>")
                    .style('opacity', 1)
            }

            // mouseover for heatmap label
            const mouseoverHeatmapLabel = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Ranks for factors in similar regions</div>")
                    .style('opacity', 1)
            }

            // mouseover for z score label
            const mouseoverSigLabel = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Significance of factor ranks in similar regions</div>")
                    .style('opacity', 1)
            }

            // mouseover for scale max line
            const mouseoverMax = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Max Rank in Similar Regions</div>")
                    .style('opacity', 1)
            }

            // mouseover for scale min line
            const mouseoverMin = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Min Rank in Similar Regions</div>")
                    .style('opacity', 1)
            }

            // mouseover for scale mean line
            const mouseoverMean = function () {
                labelTooltip
                    .html("<div style= font-size:10px; >Mean Rank in Similar Regions</div>")
                    .style('opacity', 1)
            }

            // mousemove
            const mousemoveLabel = function (e) {
                labelTooltip.style('left', (e.pageX + 20) + 'px')
                labelTooltip.style('top', (e.pageY) + 'px')
            }

            // mouseleave
            const mouseleaveLabel = function () {
                labelTooltip
                    .style('opacity', 0)
            }

            // query region bar chart label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -30)
                .attr("y", -10)
                .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Query Region")
                .on('mouseover', mouseoverQueryLabel)
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)

            // add similar regions label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -(height) / 3 - 65)
                .attr("y", -10)
                .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Similar Regions")
                .on('mouseover', mouseoverHeatmapLabel)
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)
            
            // add colorbar label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", (height) / 3 + 65)
                .attr("y", -width + margin.right - 65)
                .attr("transform", "rotate(90)")
                .style('font', '13px sans-serif')
                .text("Rank")

            // add z score bar chart label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "middle")
                .attr("x", -(height - margin.bottom) + 30)
                .attr("y", -10)
                .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Factor Z-scores")
                .on('mouseover', mouseoverSigLabel)
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
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)

            // add heatmap colorbar legend
            heatmapsvg.append("rect")
                .attr("x", width - margin.right + 10)
                .attr("y", legendY(1))
                .attr("height", legendY(0) - legendY(1))
                .attr("width", 20)
                .style("fill", "url(#heatmap-legend)")

            // add heatmap colorbar max rank
            heatmapsvg.append("line")
                .style("stroke", "green")
                .style("stroke-width", 5)
                .attr("x1", width - margin.right + 10)
                .attr("y1", legendY(maxRank))
                .attr("x2", width - margin.right + 30)
                .attr("y2", legendY(maxRank))
                .on('mouseover', mouseoverMax)
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)

            // add heatmap colorbar mean rank
            heatmapsvg.append("line")
                .style("stroke", "orange")
                .style("stroke-width", 5)
                .attr("x1", width - margin.right + 10)
                .attr("y1", legendY(meanRank))
                .attr("x2", width - margin.right + 30)
                .attr("y2", legendY(meanRank))
                .on('mouseover', mouseoverMean)
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)
            
            // add heatmap colorbar min rank
            heatmapsvg.append("line")
                .style("stroke", "red")
                .style("stroke-width", 5)
                .attr("x1", width - margin.right + 10)
                .attr("y1", legendY(minRank))
                .attr("x2", width - margin.right + 30)
                .attr("y2", legendY(minRank))
                .on('mouseover', mouseoverMin)
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)

            // add z score bar chart
            heatmapsvg.selectAll()
                .data(factorSigData)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", function(d) {if (d.sig > 0) {return sigY(d.sig)} else {return sigY(0)}})
                    .attr("width", x.bandwidth())
                    .attr("height", function(d) {if (d.sig > 0) {return sigY(0) - sigY(d.sig)} else {return sigY(d.sig) - sigY(0)}})
                    .style("fill", function(d) {if (detailTargetFactors.includes(d.originalFactorInd)) {return "#555555"} else {return "#999999"}})
                    .style("stroke-width", 4)
                    .style("stroke", "none")
                    .on('mouseover', mouseoverSig)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)

            // add line at 0 for z score bar chart
            heatmapsvg
                .append("line")
                    .style("stroke", "black")
                    .style("stroke-width", 0.5)
                    .attr("x1", margin.left)
                    .attr("y1", sigY(0))
                    .attr("x2", width - margin.right)
                    .attr("y2", sigY(0))

            // top factor color row
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

            // bottom factor color row
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

            // add query region bar chart
            heatmapsvg.selectAll()
                .data(queryDataFlat)
                .join("rect")
                    .attr("x", function(d) {return x(d.factor)})
                    .attr("y", function(d) {if (d.percentile > 0.5) {return queryY(d.percentile)} else {return queryY(0.5)}})
                    .attr("width", x.bandwidth())
                    .attr("height", function(d) {if (d.percentile > 0.5) {return queryY(0.5) - queryY(d.percentile)} else {return queryY(d.percentile) - queryY(0.5)}})
                    .style("fill", function(d) {if (detailTargetFactors.includes(d.originalFactorInd)) {return "#555555"} else {return "#999999"}})
                    .style("stroke-width", 4)
                    .style("stroke", "none")
                    .on('mouseover', mouseoverQuery)
                    .on('mousemove', mousemove)
                    .on('mouseleave', mouseleave)
                    .on('click', mouseclick)

            // line at 50% in query region bar chart
            heatmapsvg
                .append("line")
                    .style("stroke", "black")
                    .style("stroke-width", 0.5)
                    .attr("x1", margin.left)
                    .attr("y1", queryY(0.5))
                    .attr("x2", width - margin.right)
                    .attr("y2", queryY(0.5))
        }
    })
    return (
        <div id='factorHeatmap' />
    )
};

export default FactorHeatmap