import React from 'react'

function QueryRegion (props) {
  if (props.queryProps.region !== '') {
    return (
      <div className='query-region'>
        <div>{props.queryProps.region}</div>
        <div>Scale: {props.queryProps.scale / 1000}kb</div>
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default QueryRegion