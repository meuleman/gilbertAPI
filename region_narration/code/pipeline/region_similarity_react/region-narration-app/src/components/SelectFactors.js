import React, { useState } from 'react'
import factors from '../factors.json'

function SelectFactors (props) {
    
    var possibleFactorInds = props.factorProps[0]

    const [allDHS, setAllDHS] = useState(true)
    const [allEpi, setAllEpi] = useState(true)

    const handleAllChange = (metric) => {
        if (metric === 0) {
            setAllDHS(!allDHS)

        } else if (metric === 1) {
            setAllEpi(!allEpi)
        }
        
    }

    const [included, setIncluded] = useState(
        new Array(factors.length).fill(false)
    );
    const handleChange = (position) => {
        const updatedIncluded = included.map((value, index) =>
            {if (index === position){
                return !value;
            } else {
                return value;
            }}
        );
        
        setIncluded(updatedIncluded);
        var selectedFactorInds = possibleFactorInds.filter((factorInd, i) => updatedIncluded[i])

        props.factorProps[1](selectedFactorInds)
    }; 

    const [open, setOpen] = useState(
        new Array(2).fill(false)
    );
    const handleOpen = (metric) => {
        const updatedOpen = open.map((value, index) => {
            if (index === metric) {
                return !value;
            } else {
                return value;
            }
        });

        setOpen(updatedOpen)
    };

    return (
        <div>
            <div>
                <button onClick={() => handleOpen(0)}>DHS Factors for Simsearch</button>
                {open[0] ? (
                    factors.slice(0, 16).map(({name, ind}, index) => {
                        return(
                            <div key={ind}>
                                <input 
                                type="checkbox" 
                                id={ind} 
                                name={name} 
                                value={name}
                                checked={included[ind]}
                                onChange={() => handleChange(ind)}
                                /> 
                                <label className='checkbox-label'>
                                    {name}
                                </label>
                            </div>
                        )
                    })
                ) : null}
                <br></br>
                <button onClick={() => handleOpen(1)}>Chromatin State Factors for Simsearch</button>
                {open[1] ? (
                    factors.slice(16).map(({name, ind}, index) => {
                        return(
                            <div key={ind}>
                                <input 
                                type="checkbox" 
                                id={ind} 
                                name={name} 
                                value={name}
                                checked={included[ind]}
                                onChange={() => handleChange(ind)}
                                /> 
                                <label className='checkbox-label'>
                                    {name}
                                </label>
                            </div>
                        )
                    })) : null}
            </div>
        </div>
    )
}
export default SelectFactors