import React, { useState } from 'react'
// import factors from '../factors.json'

function SelectDetailLevel (props) {

    // const [detailLevel, setDetailLevel] = useState(props.detailProps.detailLevel);

    const handleChange = (value) => {
        // setDetailLevel(value)
        props.detailProps.setDetailLevel(value)
    }

    return (

        <div>
            <div className='deatil-level'>Detail Level: {props.detailProps.detailLevel}</div>
            <input 
            className='detail-slider'
            type="range" 
            id="detailRange"
            min="1" 
            max="34" 
            step="1"
            value={props.detailProps.detailLevel}
            onChange={() => handleChange(document.getElementById("detailRange").value)}
            />
        </div>
    )
}
export default SelectDetailLevel