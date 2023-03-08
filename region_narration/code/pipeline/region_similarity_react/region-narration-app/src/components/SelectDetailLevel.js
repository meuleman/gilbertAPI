import React, { useState } from 'react'
// import factors from '../factors.json'

function SelectDetailLevel (props) {

    const [detailLevel, setDetailLevel] = useState(props.detailProps.detailLevel);

    // const [mouseUp, setMouseUp] = useState(true);

    const handleChange = (value) => {
        setDetailLevel(value)
    }

    // const handleMouseUp = (value) => {
    //     // setMouseUp(value)
    //     if (value) {
    //         props.detailProps.setDetailLevel(detailLevel)
    //     }
    // }
    

    return (

        <div>
            <div className='deatil-level'>Detail Level: {detailLevel}</div>
            <input 
            className='detail-slider'
            type="range" 
            id="detailRange"
            min="1" 
            max="34" 
            step="1"
            value={detailLevel}
            // onChange={() => this.handleChange(document.getElementById("detailRange").value)}
            onChange={() => handleChange(document.getElementById("detailRange").value)}
            // onMouseUp={() => handleMouseUp(true)}
            // onMouseDown={() => handleMouseUp(false)}
            onMouseUp={() => props.detailProps.setDetailLevel(detailLevel)}
            />
        </div>
    )
}
export default SelectDetailLevel