import React, { useEffect } from 'react'
import * as d3 from 'd3'
// import { scaleLinear } from 'd3'
import factorsArr from '../../factors.json'
import QueryBarChart from './QueryBar'
import MeanRankBarChart from './MeanRankBar'
import ConvolutionBars from './ConvolutionBars'
import Legend from './Legend'
import LegendBars from './LegendBars'
import IndexButtons from './IndexButtons'
import EpilogosButtons from './EpilogosButtons'
import QueryIndexButton from './QueryIndexButton'
import QueryEpilogosButton from './QueryEpilogosButton'
import SearchFromHereButtons from './SearchFromHereButtons'
import ProductOfFactorRanks from './ProductOfFactorRanks'
import Heatmap from './Heatmap'
import TopColor from './TopColor'
import BottomColor from './BottomColor'
import MeanSquaredError from './MeanSquaredError'
// import topRegionPerFactor from '../hilbert_orders_4t11.top_regions_per_factor.json'
// import topRegionPerFactor from '../hilbert_orders_4t11.sfc.top_regions_per_factor.closest_order.json'
// import { scaleLinear, timeFormatDefaultLocale } from 'd3'

function FactorHeatmapContainer (props) {
  useEffect(() => {
    if (props.heatmapProps.data !== "") {

      let margin = { top: 100, right: 250, left: 200, bottom: 20 }
      let height = 800 - margin.top + margin.bottom
      let width = 1350 - margin.left

      let figureStartY = 90
      let colorHeight1 = 5
      let colorHeight2 = 5
      let queryLength = 60
      let heatmapLength = 400
      let meanBarLength = 60
      let sectionLabelYOffset = -10
      let sectionLabelXOffset = 180
      let gap1 = 25
      let gap2 = 25
      let browserButtonWidth = 25
      let indexButtonX = 1
      let goXAdjustment = 50

      // data
      let data = props.heatmapProps.data
      const colLabels =  data[0].factor_order.map((ind) => {
        return factorsArr[ind].name
      })

      const roiURL =  data[0].roiURL
      // console.log(data)
      let detailTargetFactors = null
      if (data[0].target_factors) {
        detailTargetFactors =  data[0].target_factors.slice(0, props.heatmapProps.detailLevel)
      } else if (data[0].selected_factors) {
        detailTargetFactors =  data[0].selected_factors
      }
      
      // console.log(detailTargetFactors)

      const factorColors =  data[0].factor_order.map((ind) => {
        return factorsArr[ind].color
      })
      factorColors[colLabels.indexOf('Quies')] = "#e0e0e0" // changing Quies color to something visable

      // colors for each factor
      const color = d3.scaleOrdinal()
        .range(factorColors)
        .domain(colLabels)

      const searchedByFactor = props.heatmapProps.searchedByFactor
      if (searchedByFactor) {
        gap1 = 0
        queryLength = 0
        colorHeight2 = 0
        goXAdjustment = 0
      }

      // set sizes
      const topColorStartY = figureStartY
      const topColorEndY = topColorStartY + colorHeight1

      const queryStartY = figureStartY + colorHeight1
      const queryEndY = queryStartY + queryLength
      const queryLabelY = queryStartY + sectionLabelYOffset
      const queryLabelX = sectionLabelXOffset

      const bottomColorStartY = queryEndY + gap1
      const bottomColorEndY = bottomColorStartY + colorHeight2

      const heatmapStartY = bottomColorEndY
      const heatmapEndY = heatmapStartY + heatmapLength
      const heatmapLabelY = heatmapStartY + sectionLabelYOffset
      const heatmapLabelX = sectionLabelXOffset

      const meanStartY = heatmapEndY + gap2
      const meanEndY = meanStartY + meanBarLength
      const meanLabelY = meanStartY + sectionLabelYOffset
      const meanLabelX = sectionLabelXOffset
      const MSELabelY = meanLabelY + 20
      const MSELabelX = sectionLabelXOffset - 100

      const plotStartX = margin.left
      const plotEndX = width - margin.right

      const epilogosButtonX = indexButtonX + browserButtonWidth + 5
      const queryButtonX = margin.left - 60
      
      const legendWidth = 20
      const legendStartX = plotEndX + 10
      const legendEndX = legendStartX + legendWidth
      const prodFactorRankX = legendEndX + 50
      const searchButtonX = prodFactorRankX + goXAdjustment
      const searchButtonWidth = 25
      const simFactorXStart = searchButtonX + searchButtonWidth + 10
      const simFactorXWidth = colorHeight1 * 34
      const simThresh = 0.9
      
      // clear document
      if (document.getElementById('factorHeatmap').getElementsByTagName('*').length !== 0) {
        // console.log(d3.selectAll('svg#heatmapSvg'))
        d3.selectAll('svg#heatmapSvg').remove()
        // d3.selectAll('#factorHeatmap').remove()
        d3.selectAll("div.tooltip").remove()
      }

      // SVG container
      const heatmapContainer = d3.select('#factorHeatmap')
      const heatmapsvg = heatmapContainer
        .append('svg')
        .attr('id', 'heatmapSvg')
        .attr('width', width + margin.left*0)
        .attr('height', height + margin.top*0)
        .append('g')
        .attr('transform', 'translate(' + margin.left*0 + ',' + margin.top*0 + ')')

      // tooltips
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

      const findOffset = function () {
        // offset of parent element
        const parentElement = document.getElementById('results');
        const parentOffset = parentElement.getBoundingClientRect()
        return parentOffset
      }

      // mousemove
      const mousemove = function (e) {
        const parentOffset = findOffset()
        tooltip.style('left', (e.clientX - parentOffset.x + 20) + 'px')
        tooltip.style('top', (e.clientY - parentOffset.y) + 'px')
      }

      // mouseleave
      const mouseleave = function (strokeWidth="none") {
        tooltip
          .style('opacity', 0)
          .style('z-index', -1)
        if (strokeWidth !== "none") {
          d3.select(this)
            .style("stroke", 'black')
            .style("stroke-width",strokeWidth)
        } else {
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
        .style('z-index', -1)
        .style('max-width', '500px')
        .style('min-width', '100px')

      // mouseover for labels
      const mouseoverLabel = function (label) {
        labelTooltip
          .html("<div style= font-size:10px; >" + label + "</div>")
          .style('opacity', 1)
          .style('z-index', 1)
      }

      // mousemove for labels
      const mousemoveLabel = function (e) {
        const parentOffset = findOffset()
        labelTooltip.style('left', (e.clientX - parentOffset.x + 20) + 'px')
        labelTooltip.style('top', (e.clientY - parentOffset.y) + 'px')
      }

      // mouseleave for labels
      const mouseleaveLabel = function () {
        labelTooltip
          .style('opacity', 0)
          .style('z-index', -1)
      }
      
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

      const { x, y, heatmapColor } = Heatmap({
        heatmapsvg: heatmapsvg,
        tooltip: tooltip,
        colLabels: colLabels,
        plotStartX: plotStartX,
        plotEndX: plotEndX,
        figureStartY: figureStartY,
        searchedByFactor: searchedByFactor,
        data: data,
        heatmapStartY: heatmapStartY,
        heatmapEndY: heatmapEndY,
        mouseoverLabel: mouseoverLabel,
        mousemoveLabel: mousemoveLabel,
        mouseleaveLabel: mouseleaveLabel,
        heatmapLabelY: heatmapLabelY,
        heatmapLabelX: heatmapLabelX,
        mousemove: mousemove,
        mouseleave: mouseleave,
        detailTargetFactors: detailTargetFactors,
      })
      MeanRankBarChart({
        heatmapsvg: heatmapsvg,
        tooltip: tooltip,
        meanStartY: meanStartY,
        meanEndY: meanEndY,
        plotStartX: plotStartX,
        meanLabelX: meanLabelX,
        meanLabelY: meanLabelY,
        mouseoverLabel: mouseoverLabel,
        mousemoveLabel: mousemoveLabel,
        mouseleaveLabel: mouseleaveLabel,
        data: data,
        x: x,
        detailTargetFactors: detailTargetFactors,
        mousemove: mousemove,
        mouseleave: mouseleave,
        searchedByFactor: searchedByFactor,
        colLabels: colLabels
      })
      const { legendY } = Legend({
        heatmapsvg: heatmapsvg,
        heatmapColor: heatmapColor,
        heatmapStartY: heatmapStartY,
        heatmapEndY: heatmapEndY,
        legendEndX: legendEndX,
        legendStartX: legendStartX,
        legendWidth: legendWidth
      })
      LegendBars({
        heatmapsvg: heatmapsvg,
        legendStartX: legendStartX,
        legendEndX: legendEndX,
        data: data,
        detailTargetFactors: detailTargetFactors,
        legendY: legendY,
        mouseoverLabel: mouseoverLabel,
        mousemoveLabel: mousemoveLabel,
        mouseleaveLabel: mouseleaveLabel,
        searchedByFactor: searchedByFactor,
      })
      IndexButtons({
        heatmapsvg: heatmapsvg,
        data: data,
        indexButtonX: indexButtonX,
        y: y,
        browserButtonWidth: browserButtonWidth,
        linkClick: linkClick,
        searchedByFactor: searchedByFactor,
      })
      EpilogosButtons({
        heatmapsvg: heatmapsvg,
        data: data,
        epilogosButtonX: epilogosButtonX,
        y: y,
        browserButtonWidth: browserButtonWidth,
        linkClick: linkClick,
        searchedByFactor: searchedByFactor,
      })
      SearchFromHereButtons({
        heatmapsvg: heatmapsvg,
        data: data,
        searchButtonX: searchButtonX,
        searchButtonWidth: searchButtonWidth,
        y: y,
        search: props.heatmapProps.search,
        searchedByFactor: searchedByFactor,
      })
      TopColor({
        heatmapsvg: heatmapsvg,
        topColorStartY: topColorStartY,
        topColorEndY: topColorEndY,
        plotStartX: plotStartX,
        colLabels: colLabels,
        x: x,
        colorHeight1: colorHeight1,
        color: color
      })
      
      if (!searchedByFactor) {
        const {queryY} = QueryBarChart({
          heatmapsvg: heatmapsvg,
          tooltip: tooltip,
          data: data,
          colLabels: colLabels,
          x: x,
          y: y,
          queryStartY: queryStartY,
          queryEndY: queryEndY,
          mousemove: mousemove,
          mouseleave: mouseleave,
          mouseoverLabel: mouseoverLabel,
          mousemoveLabel: mousemoveLabel,
          mouseleaveLabel: mouseleaveLabel,
          queryLabelX: queryLabelX,
          queryLabelY: queryLabelY,
          plotStartX: plotStartX,
          plotEndX: plotEndX,
          linkClick: linkClick,
          queryButtonX: queryButtonX,
          browserButtonWidth: browserButtonWidth,
          detailTargetFactors: detailTargetFactors
        })
        ConvolutionBars({
          heatmapsvg: heatmapsvg,
          tooltip: tooltip,
          data: data,
          colLabels: colLabels,
          simThresh: simThresh,
          simFactorXStart: simFactorXStart,
          simFactorXWidth: simFactorXWidth,
          queryY: queryY,
          y: y,
          color: color,
          mousemove: mousemove,
          mouseleave: mouseleave
        })
        QueryIndexButton({
          heatmapsvg: heatmapsvg,
          data: data,
          queryButtonX: queryButtonX,
          queryY: queryY,
          linkClick: linkClick,
          y: y,
          browserButtonWidth: browserButtonWidth,
        })
        QueryEpilogosButton({
          heatmapsvg: heatmapsvg,
          data: data,
          queryButtonX: queryButtonX,
          queryY: queryY,
          linkClick: linkClick,
          y: y,
          browserButtonWidth: browserButtonWidth,
        })
        ProductOfFactorRanks({
          heatmapsvg: heatmapsvg,
          data: data,
          prodFactorRankX: prodFactorRankX,
          queryY: queryY,
          y: y,
          mouseoverLabel: mouseoverLabel,
          mousemoveLabel: mousemoveLabel,
          mouseleaveLabel: mouseleaveLabel,
        })
        BottomColor({
          heatmapsvg: heatmapsvg,
          bottomColorStartY: bottomColorStartY,
          bottomColorEndY: bottomColorEndY,
          colLabels: colLabels,
          x: x,
          colorHeight2: colorHeight2,
          color: color,
        })
        MeanSquaredError({
          heatmapsvg: heatmapsvg,
          colLabels: colLabels,
          data: data,
          detailTargetFactors: detailTargetFactors,
          MSELabelX: MSELabelX,
          MSELabelY: MSELabelY,
        })
      }
    }
  })
  return (
    <div className='factor-heatmap' id='factorHeatmap' />
  )
};

export default FactorHeatmapContainer