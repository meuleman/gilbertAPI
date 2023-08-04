import { useState } from "react"

function QueryRegion (props) {

  if (props.queryProps.region !== '') {
    let regionStr = ''
    let regionApproxLengthStr = ''
    if (!props.queryProps.searchedByFactor) {
      const regionObj = props.queryProps.region
      regionStr = regionObj.chrom + ':' + regionObj.start + '-' + regionObj.stop
      const regionLength = regionObj.stop - regionObj.start


      if (regionLength < 1000) {
        regionApproxLengthStr = '(~' + (Math.round(regionLength / 100) * 100) + 'nt)'
      } else if (regionLength < 1000000) {
        regionApproxLengthStr = '(~' + (Math.round(regionLength / 1000)) + 'kb)'
      } else {
        regionApproxLengthStr = '(~' + (Math.round(regionLength / 1000000)) + 'Mb)'
      }

    }

    return (
      <div className='query-region'>
        <div>{regionStr} {regionApproxLengthStr}</div>
      </div>
    )
  } else {
    return (
      <div className='query-region'>
      </div>
    )
  }
}
export default QueryRegion