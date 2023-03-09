import React from 'react'

function DisplaySentence (props) {
  if (props.sentenceProps !== '') {
    return (
      <div className='sentence-box'>
        <div className='sentence-header'><b>Region Summary</b></div>
        {/* <div className='sentence-body'>{props.sentenceProps}</div> */}
        <div className='sentence-body' dangerouslySetInnerHTML={{__html: props.sentenceProps}} />
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default DisplaySentence
