import React, { useState } from 'react'
import factors from '../factors.json'

function SelectFactors (props) {
    
    var possibleFactorInds = props.factorProps[0]

    // const [allDHS, setAllDHS] = useState(true)
    // const [allEpi, setAllEpi] = useState(true)

    // const handleAllChange = (metric) => {
    //     if (metric === 0) {
    //         setAllDHS(!allDHS)

    //     } else if (metric === 1) {
    //         setAllEpi(!allEpi)
    //     }
        
    // }

    const factors_ext = factors.concat([{name: <b>All DHS Factors</b>, ind: 34}, {name: <b>All Chromatin State Factors</b>, ind: 35}]);

    const [included, setIncluded] = useState(
        new Array(factors_ext.length).fill(true)
    );
    const handleChange = (position) => {

        var updatedIncluded = new Array(factors_ext.length)

        if (position === 34) {
            updatedIncluded = included.slice(0, 16).fill(!included[34]).concat(included.slice(16))
            updatedIncluded[34] = !included[34]
        } else if (position === 35) {
            updatedIncluded = included.slice(0, 16).concat(included.slice(16, 34).fill(!included[35])).concat(included.slice(34))
            updatedIncluded[35] = !included[35]
        } else {
            updatedIncluded = included.map((value, index) =>
            { if (index === position){
                return !value;
            } else {
                return value;
            }});

            if ((updatedIncluded.slice(0, 16)).reduce((pSum, val) => pSum + val, 0) !== 16) {
                updatedIncluded[34] = false
            } else {
                updatedIncluded[34] = true
            }

            if ((updatedIncluded.slice(16, 34)).reduce((pSum, val) => pSum + val, 0) !== 18) {
                updatedIncluded[35] = false
            } else {
                updatedIncluded[35] = true
            }
        }

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
                    factors_ext.slice(0, 16).concat(factors_ext[34]).map(({name, ind}, index) => {
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
                    factors_ext.slice(16, 34).concat(factors_ext[35]).map(({name, ind}, index) => {
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

                ) : null
                }
                <div></div>
                
            </div>
        </div>
    )
}
export default SelectFactors