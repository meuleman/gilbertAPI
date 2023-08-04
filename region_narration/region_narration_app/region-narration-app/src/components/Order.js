import { useState } from "react"

function QueryRegion (props) {
  const [previousOrder, setPreviousOrder] = useState(null)

  if ((props.orderProps.region !== '') || (props.orderProps.searchedByFactor)) {
    const scales = props.orderProps.possibleScales
    const queriedLocation = props.orderProps.queriedLocation
    const setGenomicLocation = props.orderProps.setGenomicLocation
    const queryBuilding = props.orderProps.queryBuilding
    const searchedByFactor = props.orderProps.searchedByFactor
    let order = props.orderProps.currentScale
    
    // annimation for order text
    const orderText = document.getElementById('order-text');
    if (order !== previousOrder) {
      setPreviousOrder(order)
      orderText.classList.add('animation');
      setTimeout(() => {
        orderText.classList.remove('animation');
      }, 500); 
    }

    // change order
    const currentScaleInd = scales.indexOf(order)
    const handleClick = (direction) => {
      let newScaleInd = ''
      if (direction === 'in') {
        newScaleInd = Math.min(currentScaleInd + 1, scales.length - 1)
      } else if (direction === 'out') {
        newScaleInd = Math.max(currentScaleInd - 1, 0)
      }
      const newScale = scales[newScaleInd]
      if (newScale !== order) {
        if (searchedByFactor) {
          queryBuilding(newScale)
        } else {
          setGenomicLocation(queriedLocation, false, newScale, false)
        }
      }
    }

    return (
      <div className='order-container'>
        <a className='order-text' id='order-text'>Order: {order}</a>
        <button className='zoom-button' onClick={() => handleClick('out')}>-</button>
        <button className='zoom-button' onClick={() => handleClick('in')}>+</button>
      </div>

    )
  } else {
    return (
      <div className='order-container'>
        <div className='order-text' id='order-text'></div>
      </div>
    )
  }
}
export default QueryRegion