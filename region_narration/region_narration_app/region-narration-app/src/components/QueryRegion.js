import { useState } from "react"

function QueryRegion (props) {
  const [previousOrder, setPreviousOrder] = useState(null)

  if (props.queryProps.region !== '') {
    const regionObj = props.queryProps.region
    const regionStr = regionObj.chrom + ':' + regionObj.start + '-' + regionObj.stop
    const regionLength = regionObj.stop - regionObj.start
    const order = props.queryProps.scale

    let regionApproxLengthStr = ''
    if (regionLength < 1000) {
      regionApproxLengthStr = (Math.round(regionLength / 100) * 100) + 'nt'
    } else if (regionLength < 1000000) {
      regionApproxLengthStr = (Math.round(regionLength / 1000)) + 'kb'
    } else {
      regionApproxLengthStr = (Math.round(regionLength / 1000000)) + 'Mb'
    }
    
    const orderText = document.getElementById('order-text');
    if (order !== previousOrder) {
      setPreviousOrder(order)
      orderText.classList.add('animation');
      setTimeout(() => {
        orderText.classList.remove('animation');
      }, 2000); 
    }

    return (
      <div className='query-region'>
        <div>{regionStr} (~{regionApproxLengthStr})</div>
        <div className='order-text' id='order-text'>Order: {order}</div>
      </div>
    )
  } else {
    return (
      <div className='query-region'>
        <div></div>
        <div className='order-text' id='order-text'></div>
      </div>
    )
  }
}
export default QueryRegion