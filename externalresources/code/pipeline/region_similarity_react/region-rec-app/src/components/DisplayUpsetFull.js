import React from "react";

function DisplayUpsetFull(props) {
    if (props.upsetFullProps !== "") {
        return(
            <div>
                <div className="figure" dangerouslySetInnerHTML={{__html: props.upsetFullProps}} />
            </div>
        )
    } else {
        return(
            <div></div>
        )
    }
}
export default DisplayUpsetFull