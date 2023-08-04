import React, { useState } from 'react'
import factors from '../factors.json'
// import topRegionPerFactor from '../hilbert_orders_4t11.sfc.top_regions_per_factor.closest_order.json'

function TopFactorDropdown (props) {
  // let order = props.topFactorDropdownProps.scale
  // const search = props.topFactorDropdownProps.search
  // const searchedFactorInd = props.topFactorDropdownProps.searchedFactorInd
  const setQueryBuildingSelectedFactors = props.topFactorDropdownProps.setQueryBuildingSelectedFactors
  const queryBuildingSelectedFactors = props.topFactorDropdownProps.queryBuildingSelectedFactors

  const [open, setOpen] = useState(false);

  let checked = Array(factors.length).fill(false).map((value, i) => {
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
    // checkboxes
    // let newChecked = checked.map((value, i) => {
    //   if (i === ind) {
    //     return !value
    //   } else {
    //     return value
    //   }
    // })
    // setChecked(newChecked)

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
          {factors.map(({name, ind}) => {
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