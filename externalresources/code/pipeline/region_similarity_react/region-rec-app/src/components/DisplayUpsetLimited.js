import React from "react";

function DisplayUpsetLimited(props) {
    if (props.upsetLimitedProps !== "") {
        return(
            <div>
                <div className="figure" dangerouslySetInnerHTML={{__html: props.upsetLimitedProps}} />
            </div>
        )
    } else {
        return(
            <div></div>
        )
    }
}
export default DisplayUpsetLimited