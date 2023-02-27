import React, { useState } from 'react'
import factors from '../factors.json'

function SelectFactors (props) {
    var factorOrder = [...Array(factors.length).keys()]
    const data = props.factorProps.data
    // var autoFactorSelection = true


    if (data !== '') {
        factorOrder = data[0].query_factor_order
        // console.log(factorOrder)
    }
    const possibleFactorInds = props.factorProps.possibleFactors

    // const [allDHS, setAllDHS] = useState(true)
    // const [allEpi, setAllEpi] = useState(true)

    // const handleAllChange = (metric) => {
    //     if (metric === 0) {
    //         setAllDHS(!allDHS)

    //     } else if (metric === 1) {
    //         setAllEpi(!allEpi)
    //     }
        
    // }

    // const factors_ext = factors.concat([{name: <b>All DHS Factors</b>, ind: 34}, {name: <b>All Chromatin State Factors</b>, ind: 35}]);
    // const factors_ext = factors.concat([{name: <b>All Factors</b>, ind: 34}]);
    
    
    const factors_ext = factorOrder.map((i) => {
        return (factors[i])
    }).concat([{name: <b>All Factors</b>, ind: 34}]);

    // const [included, setIncluded] = useState(
    //     new Array(factors_ext.length).fill(true)
    // );
    

    // console.log(props)
    // props.factorProps.preselectFactors(setIncluded)
    
    // var preIncluded = new Array(factors_ext.length).fill(false)
    // preIncluded = preIncluded.map((value, index) => 
    // { if (props.factorProps.selectedFactors.includes(index)) {
    //     return !value
    // } else {
    //     return value;
    // }})
    // if (preIncluded.slice(0, 34).reduce((pSum, val) => pSum + val, 0) === 34) {
    //     preIncluded[34] = true
    // }
    // console.log(preIncluded)
    // if (props.factorProps.preselectUpdate) {
    //     setIncluded(preIncluded)
    //     props.factorProps.setPreselectUpdate()
    // }
    


    

    const handleChange = (position) => {
        props.factorProps.setAutoFactorSelection(false)
        var updatedIncluded = new Array(factors_ext.length)

        // if (position === 34) {
        //     updatedIncluded = included.slice(0, 16).fill(!included[34]).concat(included.slice(16))
        //     updatedIncluded[34] = !included[34]
        // } else if (position === 35) {
        //     updatedIncluded = included.slice(0, 16).concat(included.slice(16, 34).fill(!included[35])).concat(included.slice(34))
        //     updatedIncluded[35] = !included[35]

        if (position === 34) {
            updatedIncluded = props.factorProps.included.fill(!props.factorProps.included[34])
        } else {
            updatedIncluded = props.factorProps.included.map((value, index) =>
            { if (index === position){
                return !value;
            } else {
                return value;
            }});

            // if ((updatedIncluded.slice(0, 16)).reduce((pSum, val) => pSum + val, 0) !== 16) {
            //     updatedIncluded[34] = false
            // } else {
            //     updatedIncluded[34] = true
            // }

            // if ((updatedIncluded.slice(16, 34)).reduce((pSum, val) => pSum + val, 0) !== 18) {
            //     updatedIncluded[35] = false
            // } else {
            //     updatedIncluded[35] = true
            // }

            if ((updatedIncluded.slice(0, 34)).reduce((pSum, val) => pSum + val, 0) !== 34) {
                updatedIncluded[34] = false
            } else {
                updatedIncluded[34] = true
            }
        }

        props.factorProps.setIncluded(updatedIncluded);
        // included = updatedIncluded

        var selectedFactorInds = possibleFactorInds.filter((factorInd, i) => updatedIncluded[i])

        props.factorProps.setSelectedFactors(selectedFactorInds)
        
    }; 

    // const [open, setOpen] = useState(
    //     new Array(2).fill(false)
    // );

    const [open, setOpen] = useState(
        new Array(1).fill(false)
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
        // <div className='factor-select-grid'>
        //     <div className='dropdown'>
        //         <button className='dropdown-button' onClick={() => handleOpen(0)}>DHS Factors for Simsearch</button>
        //         {open[0] ? (
        //             <div className='dropdown-option-container'>
        //                 {factors_ext.slice(0, 16).concat(factors_ext[34]).map(({name, ind}, index) => {
        //                     return(
        //                         <div className='dropdown-check' key={ind}>
        //                             <input 
        //                             type='checkbox' 
        //                             id={ind} 
        //                             name={name} 
        //                             value={name}
        //                             checked={included[ind]}
        //                             onChange={() => handleChange(ind)}
        //                             /> 
        //                             <label className='dropdown-label'>
        //                                 {name}
        //                             </label>
        //                         </div>
        //                     )
        //                 })}
        //             </div>
        //         ) : null}
        //     </div>
        //     <div className='dropdown'>
        //         <button className='dropdown-button' onClick={() => handleOpen(1)}>Chromatin State Factors for Simsearch</button>
        //         {open[1] ? (
        //             <div className='dropdown-option-container'>
        //             {factors_ext.slice(16, 34).concat(factors_ext[35]).map(({name, ind}, index) => {
        //                 return(
        //                     <div className='dropdown-check' key={ind}>
        //                         <input 
        //                         type='checkbox'
        //                         id={ind} 
        //                         name={name} 
        //                         value={name}
        //                         checked={included[ind]}
        //                         onChange={() => handleChange(ind)}
        //                         /> 
        //                         <label className='dropdown-label'>
        //                             {name}
        //                         </label>
        //                     </div>
        //                 )
        //             })}
        //             </div>

        //         ) : null}
                
        //     </div>
        // </div>
        <div className='dropdown'>
                <button className='dropdown-button' onClick={() => handleOpen(0)}>Factors for Simsearch</button>
                {open[0] ? (
                    <div className='dropdown-option-container'>
                        {factors_ext.map(({name, ind}, index) => {
                            return(
                                <div className='dropdown-check' key={ind}>
                                    <input 
                                    type='checkbox' 
                                    id={ind} 
                                    name={name} 
                                    value={name}
                                    checked={props.factorProps.included[ind]}
                                    onChange={() => handleChange(ind)}
                                    /> 
                                    <label className='dropdown-label'>
                                        {name}
                                    </label>
                                </div>
                            )
                        })}
                    </div>
                ) : null}
            </div>
    )
}
export default SelectFactors