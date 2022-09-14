import React from "react";

function DisplaySentence(props) {
    if (props.sentenceProps !== "") {
        return(
            <div className="boxed">
                <div className="sentence-header"><b>Sentence Describing Region:</b></div>
                <div className="sentence-body">{props.sentenceProps}</div>
            </div>
        )
    } else {
        return(
            <div></div>
        )
    }
}
export default DisplaySentence