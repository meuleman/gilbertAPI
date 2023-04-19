import React, { useEffect } from 'react'
import * as d3 from 'd3'
// import { scaleLinear } from 'd3'
import factorsArr from '../factors.json'
// import { scaleLinear, timeFormatDefaultLocale } from 'd3'

function FactorHeatmap (props) {
    useEffect(() => {
        if (props.heatmapProps.data !== "") {

            const margin = { top: 100, right: 275, left: 150, bottom: 100 }
            const height = 700
            const width = 1200
            const sigHeight = 60
            const zMax = 20
            const zMin = -10
            const queryHeight = 60
            const colorHeight = 5
            const gap1 = 50
            const gap2 = 25
            const browserButtonWidth = 25
            const indexButtonX = -50
            const epilogosButtonX = indexButtonX + browserButtonWidth + 5
            const queryButtonX = margin.left - 60
            
            const legendX = width - margin.right + 10
            const legendWidth = 20
            const prodFactorRankX = legendX + legendWidth + 50
            const searchButtonX = prodFactorRankX + 50
            const searchButtonWidth = 25
            const simFactorXStart = searchButtonX + searchButtonWidth + 10
            const simFactorXWidth = colorHeight * 34
            const topQueryFactorXStart = legendX
            const topQueryFactorXWidth = colorHeight * 34
            const simThresh = 0.9

            // data
            const data = props.heatmapProps.data
            
            const queryData = data[0]
            
            const colLabels = queryData.query_factor_order.map((ind) => {
                return factorsArr[ind].name
            })
            const roiURL = queryData.roiURL
            const detailTargetFactors = queryData.target_factors.slice(0, props.heatmapProps.detailLevel)

            const factorColors = queryData.query_factor_order.map((ind) => {
                return factorsArr[ind].color
            })
            factorColors[colLabels.indexOf('Quies')] = "#e0e0e0" // changing Quies color to something visable

            const dataFlat = (data.flatMap((region) => {
                const rank = region.rank
                const percentiles = region.percentiles
                const coords = region.coordinates
                const score = region.score
                const regionInd = region.region_ind

                return (
                    percentiles.map((p, index) => {
                        return ({
                            percentile: p,
                            factor: colLabels[index],
                            rank: rank,
                            coords: coords,
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

            let allSimRanks = heatmapData.map((d) => {
                if (detailTargetFactors.includes(d.originalFactorInd)) {
                    return d.percentile
                } else {
                    return null
                }
            })
            // console.log(allSimRanks)
            allSimRanks = allSimRanks.filter((value) => {
                return value != null
            })
            // console.log(allSimRanks)
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
                    factor: colLabels[index],
                    originalFactorInd: queryData.query_factor_order[index]
                })
            })

            const prodFactorRanks = data.map((d) => {
                return {
                    rank: d.rank,
                    prodFactorRank: d.prod_rank
                }
            })



            // clear document
            if (document.getElementById('factorHeatmap').getElementsByTagName('*').length !== 0) {
                d3.select('svg').remove()
                d3.selectAll("div.tooltip").remove()
            }

            // container
            const heatmapContainer = d3.select('#factorHeatmap')
            const heatmapsvg = heatmapContainer
                .append('svg')
                .attr('width', width + margin.left)
                .attr('height', height + margin.top)
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

            // Factors similar to query
            const simFactorX = d3.scaleBand()
                .domain([...Array(colLabels.length).keys()])
                .range([simFactorXStart, simFactorXStart + simFactorXWidth])


            const topQueryFactorX = d3.scaleBand()
                .domain([...Array(colLabels.length).keys()])
                .range([topQueryFactorXStart, topQueryFactorXStart + topQueryFactorXWidth])

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
                .range([queryHeight + gap1 - colorHeight, queryHeight + gap1])
                .padding(0.1)
            heatmapsvg.append('g')

            // console.log(data)
            // heatmap y axis
            const y = d3.scaleBand()
                .domain(heatmapData.map(d => d.rank))
                .range([queryHeight + gap1, height - margin.bottom - sigHeight])
                .padding(0.1)
            heatmapsvg.append('g')
                .style('font', '10px sans-serif')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(y).tickFormat((regionInd) => {
                    const coords = data[regionInd].coordinates
                    return coords //+ ' - ' + regionInd
                }))

            heatmapsvg.append('g')
                .style('font', '10px sans-serif')
                .attr("transform", "translate(-60, 0)")
                .call(d3.axisLeft(y).tickSize(0))
                .call((g) => {return g.select(".domain").remove()})

            // z score y axis
            const sigTicks = [...Array(((zMax - zMin) / 10) + 1).keys()].map((val) => {
                return val * 10 - 10
            })
            const sigY = d3.scaleLinear()
                .domain([zMax, zMin])
                .range([height - margin.bottom - sigHeight + gap2, height - margin.bottom + gap2])
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + margin.left + ", 0)")
                .call(d3.axisLeft(sigY).tickValues(sigTicks))
            
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
                .range([queryHeight + gap1, height - margin.bottom - sigHeight])
            heatmapsvg.append('g')
                .style('font', '8px sans-serif')
                .attr("transform", "translate(" + (legendX + legendWidth) + ", 0)")
                .call(d3.axisRight(legendY).tickValues([0, 0.25, 0.50, 0.75, 1.0]).tickFormat(d3.format(".0%")))

            // tooltip
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
                .style('z-index', -1)
                .style('max-width', '500px')

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

            // mouseover for z score bar plot
            const mouseoverSig = function (d, cellData) {
                const Zscore = Math.round(cellData.zscore * 10000) / 10000
                const factor = cellData.factor
                tooltip
                    .html(factor + " Z-score: " + Zscore)
                    .style('opacity', 1)
                    .style('z-index', 1)
                d3.select(this)
                    .style("stroke", "black")
                    .style("stroke-width", "1")
            }

            // mouseover for query region bar chart
            const mouseoverQuery = function (d, cellData) {
                const factor = cellData.factor
                const percentile = Math.round((cellData.percentile * 100) * 10) / 10
                tooltip
                    .html(factor + ": " + percentile + "%")
                    .style('opacity', 1)
                    .style('z-index', 1)
                d3.select(this)
                    .style("stroke", "black")
                    .style("stroke-width", "1")
            }

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

            // mousemove
            const mousemove = function (e) {
                tooltip.style('left', (e.pageX + 20) + 'px')
                tooltip.style('top', (e.pageY) + 'px')
            }

            // mouseleave
            const mouseleave = function () {
                tooltip
                    .style('opacity', 0)
                    .style('z-index', -1)
                d3.select(this)
                    .style("stroke", "none")
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
                .style('z-index', -1)
                .style('max-width', '500px')

            // mouseover for labels
            const mouseoverLabel = function (label) {
                labelTooltip
                    .html("<div style= font-size:10px; >" + label + "</div>")
                    .style('opacity', 1)
                    .style('z-index', 1)
            }

            // mousemove for labels
            const mousemoveLabel = function (e) {
                labelTooltip.style('left', (e.pageX + 20) + 'px')
                labelTooltip.style('top', (e.pageY) + 'px')
            }

            // mouseleave for labels
            const mouseleaveLabel = function () {
                labelTooltip
                    .style('opacity', 0)
                    .style('z-index', -1)
            }

            // query region bar chart label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "end")
                .attr("y", -10)
                .attr("x", margin.left - 20)
                // .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Query Region")
                .on('mouseover', () => mouseoverLabel("Ranks for factors in the query region"))
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)

            // add similar regions label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "end")
                .attr("y", 100)
                .attr("x", margin.left - 20)
                // .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Similar Regions")
                .on('mouseover', () => mouseoverLabel("Ranks for factors in similar regions"))
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)
            
            // // add colorbar label
            // heatmapsvg.append("text")
            //     .attr("class", "y label")
            //     .attr("text-anchor", "middle")
            //     .attr("x", (height) / 2)
            //     .attr("y", -legendX - 55)
            //     .attr("transform", "rotate(90)")
            //     .style('font', '13px sans-serif')
            //     .text("Rank")

            // add z score bar chart label
            heatmapsvg.append("text")
                .attr("class", "y label")
                .attr("text-anchor", "end")
                .attr("y", 560)
                .attr("x", margin.left - 20)
                // .attr("transform", "rotate(-90)")
                .style('font', '13px sans-serif')
                .text("Rank Enrichment")
                .on('mouseover', () => mouseoverLabel("Z-score of factor ranks in similar regions"))
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

            
            // const simFactorData = heatmapData.filter((d) => {
            //     return d.percentile * queryDataFlat[d.factorInd].percentile > simThresh
            // })
            const simFactorData = dataFlat.filter((d) => {
                return d.percentile * queryDataFlat[d.factorInd].percentile > simThresh
            }) // includes query region data

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
                        if (d.rank === 0) {
                            return topQueryFactorX(factorNum)
                        } else {
                            return simFactorX(factorNum)
                        }
                        
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
                    .on("mouseleave", mouseleave)

            
            const linkClick = (d) => {
                const buttonData = d.srcElement.__data__
                const coords = buttonData.coords
                const seqname = coords.split(":")[0]
                const start  = coords.split(":")[1].split("-")[0]
                const stop = coords.split(":")[1].split("-")[1]
                const browser = buttonData.browser
                let url = null
                if (browser === "index") {
                    url = "https://index.altius.org:1234/?application=viewer&mode=pd&chrLeft=" + seqname + "&chrRight=" + seqname + "&start=" + start + "&stop=" + stop + "&dtt=200000&dtc=av3&gatt=cv&gac=off"
                } else if (browser === "epilogos") {
                    url = "https://epilogos.altius.org:1234/?application=viewer&sampleSet=vC&mode=single&genome=hg38&model=18&complexity=KL&group=all&chrLeft=" + seqname + "&chrRight=" + seqname + "&start=" + start + "&stop=" + stop + "&gatt=cv"
                } else if (browser === "roiURL") {
                    url = roiURL
                }
                if ((url !== "roiURL not available") & (url !== null)) {
                    window.open(url, "_blank")
                }
            }
                
            // add index buttons
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("rect")
                    .attr("x", indexButtonX)
                    .attr("y", function(d) {return y(d.rank)})
                    .attr("width", browserButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "index"
                    }})
                    .on("click", (d) => {linkClick(d)})
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("text")
                    .attr("x", indexButtonX + browserButtonWidth / 2)
                    .attr("y", function(d) {return y(d.rank) + y.bandwidth() - 5})
                    .attr("text-anchor", "middle")
                    .style('font', '15px sans-serif')
                    .style("cursor", "pointer")
                    .text("i")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "index"
                    }})
                    .on("click", (d) => {linkClick(d)})
            
            // add epilogos buttons
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("rect")
                    .attr("x", epilogosButtonX)
                    .attr("y", function(d) {return y(d.rank)})
                    .attr("width", browserButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "epilogos"
                    }})
                    .on("click", (d) => {linkClick(d)})
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("text")
                    .attr("x", epilogosButtonX + browserButtonWidth / 2)
                    .attr("y", function(d) {return y(d.rank) + y.bandwidth() - 5})
                    .attr("text-anchor", "middle")
                    .style('font', '15px sans-serif')
                    .style("cursor", "pointer")
                    .text("e")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "epilogos"
                    }})
                    .on("click", (d) => {linkClick(d)})

            // add query index button
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("rect")
                    .attr("x", queryButtonX)
                    .attr("y", queryY(0.5) - 3 * y.bandwidth() / 2 - 2.5)
                    .attr("width", browserButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "index"
                    }})
                    .on("click", (d) => {linkClick(d)})
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("text")
                    .attr("x", queryButtonX + browserButtonWidth / 2)
                    .attr("y", queryY(0.5) - y.bandwidth() / 2 - 7.5)
                    .attr("text-anchor", "middle")
                    .style('font', '15px sans-serif')
                    .style("cursor", "pointer")
                    .text("i")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "index"
                    }})
                    .on("click", (d) => {linkClick(d)})
            
            // add query epilogos button
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("rect")
                    .attr("x", queryButtonX)
                    .attr("y", queryY(0.5) - y.bandwidth() / 2)
                    .attr("width", browserButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "epilogos"
                    }})
                    .on("click", (d) => {linkClick(d)})
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("text")
                    .attr("x", queryButtonX + browserButtonWidth / 2)
                    .attr("y", queryY(0.5) + (y.bandwidth() / 2) - 5)
                    .attr("text-anchor", "middle")
                    .style('font', '15px sans-serif')
                    .style("cursor", "pointer")
                    .text("e")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "epilogos"
                    }})
                    .on("click", (d) => {linkClick(d)})

            // add roiURL button
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("rect")
                    .attr("x", queryButtonX)
                    .attr("y", queryY(0.5) + y.bandwidth() / 2 + 2.5)
                    .attr("width", browserButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "roiURL"
                    }})
                    .on("click", (d) => {linkClick(d)})
            heatmapsvg.selectAll()
                .data(data.slice(0, 1))
                .join("text")
                    .attr("x", queryButtonX + browserButtonWidth / 2)
                    .attr("y", queryY(0.5) + 3 * y.bandwidth() / 2 - 2.5)
                    .attr("text-anchor", "middle")
                    .style('font', '15px sans-serif')
                    .style("cursor", "pointer")
                    .text("r")
                    .datum(function(d) {return {
                        "coords": d.coordinates,
                        "browser": "roiURL"
                    }})
                    .on("click", (d) => {linkClick(d)})


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

            // search from here button
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("rect")
                    .attr("x", searchButtonX)
                    .attr("y", function(d) {return y(d.rank)})
                    .attr("width", searchButtonWidth)
                    .attr("height", y.bandwidth())
                    .style("fill", "#E0E0E0")
                    .style("stroke-width", 1)
                    .style("stroke", "black")
                    .style("rx", 1)
                    .style("ry", 1)
                    .style("cursor", "pointer")
                    .datum(function(d) {return {
                        "regionCoords": d.coordinates
                    }})
                    .on("click", (d) => {
                        props.heatmapProps.search(d.srcElement.__data__.regionCoords)
                    })
            heatmapsvg.selectAll()
                .data(data.slice(1))
                .join("text")
                    .attr("x", searchButtonX + searchButtonWidth / 2)
                    .attr("y", function(d) {return y(d.rank) + y.bandwidth() - 5})
                    .attr("text-anchor", "middle")
                    .style('font', '10px sans-serif')
                    .style("cursor", "pointer")
                    .text("Go!")
                    .datum(function(d) {return {
                        "regionCoords": d.coordinates
                    }})
                    .on("click", (d) => {
                        props.heatmapProps.search(d.srcElement.__data__.regionCoords)
                    })

            // add heatmap colorbar legend
            heatmapsvg.append("rect")
                .attr("x", legendX)
                .attr("y", legendY(1))
                .attr("height", legendY(0) - legendY(1))
                .attr("width", legendWidth)
                .style("fill", "url(#heatmap-legend)")

            // add heatmap colorbar max rank
            heatmapsvg.append("line")
                .style("stroke", "green")
                .style("stroke-width", 5)
                .attr("x1", legendX)
                .attr("y1", legendY(maxRank))
                .attr("x2", legendX + legendWidth)
                .attr("y2", legendY(maxRank))
                .on('mouseover', () => mouseoverLabel("Max Rank in Similar Regions"))
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)

            // add heatmap colorbar mean rank
            heatmapsvg.append("line")
                .style("stroke", "orange")
                .style("stroke-width", 5)
                .attr("x1", legendX)
                .attr("y1", legendY(meanRank))
                .attr("x2", legendX + legendWidth)
                .attr("y2", legendY(meanRank))
                .on('mouseover', () => mouseoverLabel("Mean Rank in Similar Regions"))
                .on('mousemove', mousemoveLabel)
                .on('mouseleave', mouseleaveLabel)
            
            // add heatmap colorbar min rank
            heatmapsvg.append("line")
                .style("stroke", "red")
                .style("stroke-width", 5)
                .attr("x1", legendX)
                .attr("y1", legendY(minRank))
                .attr("x2", legendX + legendWidth)
                .attr("y2", legendY(minRank))
                .on('mouseover', () => mouseoverLabel("Min Rank in Similar Regions"))
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
                    // .on('click', mouseclick)

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
                    // .on('click', mouseclick)

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
        <div className='factor-heatmap' id='factorHeatmap' />
    )
};

export default FactorHeatmap