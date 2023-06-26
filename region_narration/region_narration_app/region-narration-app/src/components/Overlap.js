function Overlap (props) {
  const queriedLocation = props.overlapProps.queriedLocation
  const returnedLocation = props.overlapProps.returnedLocation

  if (returnedLocation !== '') {
    const bpMissingLeft = Math.max(0, returnedLocation.start - queriedLocation.start)
    const bpMissingRight = Math.max(0, queriedLocation.stop - returnedLocation.stop)
    const queriedLength = queriedLocation.stop - queriedLocation.start

    const percentContained = Math.max(
      0, Math.round(
        (
          (
            queriedLength - (
              bpMissingLeft + bpMissingRight
              )
            ) / queriedLength * 100
          ) * 100
        ) / 100
      )

    return (
      <div className='overlap'>
        {percentContained}% of Queried Region
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default Overlap