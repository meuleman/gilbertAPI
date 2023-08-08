import React, { useState } from 'react'
import factors from '../factors.json'
// import topRegionPerFactor from '../hilbert_orders_4t11.sfc.top_regions_per_factor.closest_order.json'

function TopFactorDropdown (props) {
  // let order = props.topFactorDropdownProps.scale
  // const search = props.topFactorDropdownProps.search
  // const searchedFactorInd = props.topFactorDropdownProps.searchedFactorInd
  const setQueryBuildingSelectedFactors = props.topFactorDropdownProps.setQueryBuildingSelectedFactors
  const queryBuildingSelectedFactors = props.topFactorDropdownProps.queryBuildingSelectedFactors
  const factorOrder = props.topFactorDropdownProps.factorOrder

  const [open, setOpen] = useState(false);
  const factorsOrdered = factorOrder.map((i) => {
    return factors[i]
  })

  let checked = Array(factorsOrdered.length).fill(false).map((value, i) => {
    if (queryBuildingSelectedFactors.includes(i)) {
      return true
    } else {
      return false
    }
  })
  // const [checked, setChecked] = useState(initArray)

  const handleButtonClick = function () {
    setOpen(!open)
  }

  const handleCheck = function (ind) {
    // factors for query building
    let newQueryBuildingSelectedFactors = null
    if (queryBuildingSelectedFactors.includes(ind)) {
      newQueryBuildingSelectedFactors = queryBuildingSelectedFactors.filter((value) => value !== ind)
    } else {
      newQueryBuildingSelectedFactors = queryBuildingSelectedFactors
      newQueryBuildingSelectedFactors.push(ind)
    }
    setQueryBuildingSelectedFactors(newQueryBuildingSelectedFactors)

    // close dropdown
    setOpen(false)
  }

  let buttonText = "Search by Factors"

  return (
    <div className='dropdown'>
      <button className='dropdown-button' onClick={() => handleButtonClick()} >{buttonText}</button>
      {open ? (
        <div className='dropdown-container'>
          {factorsOrdered.map(({name, ind}) => {
            return(
              // <div className='dropdown-item' key={ind} onClick={() => handleFactorClick(ind)}>
              <div className='dropdown-item' key={ind}>
                <input 
                  className='dropdown-checkbox'
                  type='checkbox' 
                  id={ind} 
                  name={name} 
                  value={name}
                  checked={checked[ind]}
                  onChange={() => handleCheck(ind)}
                /> 
                <label className='dropdown-factor-label'>
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
export default TopFactorDropdown