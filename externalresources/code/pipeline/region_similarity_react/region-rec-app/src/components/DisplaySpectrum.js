import React from "react";

function DisplaySpectrum(props) {
    if (props.spectrumProps !== "") {
        return(
            <div>
                <div className="figure" dangerouslySetInnerHTML={{__html: props.spectrumProps}} />
            </div>
        )
    } else {
        return(
            <div></div>
        )
    }
}
export default DisplaySpectrum