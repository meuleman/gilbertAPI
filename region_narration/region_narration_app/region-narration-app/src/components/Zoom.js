import React from 'react'

function Zoom (props) {
    const scales = props.zoomProps.possibleScales
    const currentLocation = props.zoomProps.currentLocation
    const currentScale = props.zoomProps.currentScale
    const setGenomicLocation = props.zoomProps.setGenomicLocation

    const currentScaleInd = scales.indexOf(currentScale)

    const handleClick = (direction) => {
        // console.log(props.zoomProps)
        let newScaleInd = ''
        if (direction === 'in') {
            newScaleInd = Math.min(currentScaleInd + 1, scales.length - 1)
        } else if (direction === 'out') {
            newScaleInd = Math.max(currentScaleInd - 1, 0)
        }
        const newScale = scales[newScaleInd]
        
        if (newScale !== currentScale) {
            setGenomicLocation(currentLocation, false, newScale)
        }
    }

    return (

        <div className='zoom'>
            <span className='zoom-text'>Zoom:</span>
            <button className='zoom-button' onClick={() => handleClick('in')}>+</button>
            <button className='zoom-button' onClick={() => handleClick('out')}>-</button>
        </div>
    )
}
export default Zoom