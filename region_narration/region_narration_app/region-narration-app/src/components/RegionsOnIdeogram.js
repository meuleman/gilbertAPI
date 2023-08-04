import React, { useEffect } from 'react'
// import * as d3 from 'd3'
import Ideogram from 'ideogram';


function RegionsOnIdeogram (props) {
  useEffect(() => {
    if (props.ideogramProps.data !== '') {
      const data = props.ideogramProps.data
      // const ideogramElement = document.getElementById("regionsOnIdeogram");
      // ideogramElement.style.paddingTop = "90px"
      // if (props.ideogramProps.searchedByFactor) {
      //   ideogramElement.style.paddingTop = "0px"
      // }
      // clear document
      // if (document.getElementById('factorHeatmap').getElementsByTagName('*').length !== 0) {
      //   d3.select('svg').remove()
      //   d3.selectAll("div.tooltip").remove()
      // }

      const regionAnnotations = data.map((d) => {
        let c = ''
        let n = ''
        let start = parseInt(d.coordinates.split(':')[1].split('-')[0])
        if (start === 0) {
          start = 1
        }
        let stop = parseInt(d.coordinates.split(':')[1].split('-')[1])
        if (!props.ideogramProps.searchedByFactor) {
          if (d.rank === 0) {
            c = '#F00'
            n = 'Query Region'
          } else {
            c = '#8D4'
            n = 'Similar Region ' + d.rank
          }
        } else {
          c = '#8D4'
          n = 'Top Region ' + (parseInt(d.rank) + 1)
        }
        
        return {
          name: n,
          chr: d.coordinates.split(':')[0].split('chr')[1],
          start: start,
          stop: stop,
          color: c
        }
      })
      let chrWidth = 7
      // if (props.ideogramProps.searchedByFactor) {
      //   chrWidth = 7
      // }

      var ideogram = new Ideogram({
        organism: 'human',
        annotations: regionAnnotations,
        container: '#regionsOnIdeogram',
        assembly: 'GRCh38',
        orientation: 'horizontal',
        // rows: 2,
        chrHeight: 100,
        chrWidth: chrWidth,
        annotationHeight: 3,
        // otatable: false,
      });
      ideogram.onClickAnnot = (annot) => {
        console.log(annot)
      } 
        
    }
  })
  // d3 = Ideogram.d3;
  

  return (
    <div className='regions-on-ideogram' id='regionsOnIdeogram' />
  )
};

export default RegionsOnIdeogram