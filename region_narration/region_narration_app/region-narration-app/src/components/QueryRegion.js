function QueryRegion (props) {
  if (props.queryProps.region !== '') {
    const regionObj = props.queryProps.region
    const regionStr = regionObj.chrom + ':' + regionObj.start + '-' + regionObj.stop
    const regionLength = regionObj.stop - regionObj.start
    let regionApproxLengthStr = ''
    if (regionLength < 1000) {
      regionApproxLengthStr = (Math.round(regionLength / 100) * 100) + 'nt'
    } else if (regionLength < 1000000) {
      regionApproxLengthStr = (Math.round(regionLength / 1000)) + 'kb'
    } else {
      regionApproxLengthStr = (Math.round(regionLength / 1000000)) + 'Mb'
    }
    // console.log(regionLength, regionApproxLengthStr)
    return (
      <div className='query-region'>
        <div>{regionStr} (~{regionApproxLengthStr})</div>
        <div>Order: {props.queryProps.scale}</div>
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default QueryRegion