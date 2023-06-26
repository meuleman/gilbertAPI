function Zoom (props) {
  const scales = props.zoomProps.possibleScales
  const queriedLocation = props.zoomProps.queriedLocation
  const returnedLocation = props.zoomProps.returnedLocation
  const currentScale = props.zoomProps.currentScale
  const setGenomicLocation = props.zoomProps.setGenomicLocation

  if (returnedLocation !== '') {
    const currentScaleInd = scales.indexOf(currentScale)

    const handleClick = (direction) => {
      let newScaleInd = ''
      if (direction === 'in') {
        newScaleInd = Math.min(currentScaleInd + 1, scales.length - 1)
      } else if (direction === 'out') {
        newScaleInd = Math.max(currentScaleInd - 1, 0)
      }
      const newScale = scales[newScaleInd]
      
      if (newScale !== currentScale) {
        setGenomicLocation(queriedLocation, false, newScale, false)
      }
    }

    return (
      <div className='zoom'>
        <span className='zoom-text'>Zoom:</span>
        <button className='zoom-button' onClick={() => handleClick('in')}>+</button>
        <button className='zoom-button' onClick={() => handleClick('out')}>-</button>
      </div>
    )

  } else {
    return (
      <div />
    )
  }
}
export default Zoom