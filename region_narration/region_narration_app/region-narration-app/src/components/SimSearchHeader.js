import React from 'react'

function SimSearchHeader (props) {
  if (props.simSearchHeaderProps !== '') {
    return (
      <div className='boxed'>
        <div className='sentence-header'><b>Similarity Search</b></div>
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default SimSearchHeader