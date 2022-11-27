import React from 'react'

function DisplayBar (props) {
  if (props.barProps !== '') {
    /* let blob = new Blob([props.barProps], {type: 'image/svg+xml'});
        let url = URL.createObjectURL(blob); */
    const url = props.barProps
    return (
    /* <div className="bar">
                <svg
                //height="100%"
                //width="100%"
                //className="scaling-svg"
                viewBox="0 0 1300 1600"
                //preserveAspectRatio="xMinYMin meet"
                style={{
                    position: "absolute",
                    height: "100%",
                    width: "100%",
                    left: "0",
                    top: "0"
                }}
                dangerouslySetInnerHTML={{__html: props.barProps}} />
            </div> */
      <div className='bar'>
        <img
          src={url}
                    // viewBox="0 0 1200 1600"
          width='100%'
          height='auto'
          alt=''
        />
      </div>
    )
  } else {
    return (
      <div className='bar' />
    )
  }
}
export default DisplayBar
