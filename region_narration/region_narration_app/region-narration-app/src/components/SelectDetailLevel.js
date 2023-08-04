function SelectDetailLevel (props) {
  const detailLevel = props.detailProps.detailLevel
  const returnedLocation = props.detailProps.returnedLocation
  const maxDetailLevel = props.detailProps.maxDetailLevel
  const setDetailLevel = props.detailProps.setDetailLevel
  const searchedByFactor = props.detailProps.searchedByFactor

  if ((returnedLocation !== '') && (!searchedByFactor)) {
    const handleChange = (value) => {
      setDetailLevel(value)
    }

    return (
      <div className='detail-parent'>
        <div className='detail-level'>Detail Level: {detailLevel}</div>
        <input 
        className='detail-slider'
        type="range" 
        id="detailRange"
        min="1" 
        max={maxDetailLevel}
        step="1"
        value={detailLevel}
        onChange={() => handleChange(document.getElementById("detailRange").value)}
        />
      </div>
    )

  } else {
    return (
      <div />
    )
  }
}
export default SelectDetailLevel