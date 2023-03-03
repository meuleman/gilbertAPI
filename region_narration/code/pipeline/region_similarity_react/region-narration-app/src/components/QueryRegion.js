import React from 'react'

function QueryRegion (props) {
  if (props.queryProps !== '') {
    return (
      <div>
        <div>{props.queryProps}</div>
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default QueryRegion