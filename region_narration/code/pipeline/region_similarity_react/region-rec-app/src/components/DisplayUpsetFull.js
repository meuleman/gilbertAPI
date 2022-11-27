import React from "react";
import '../App.css'

function DisplayUpsetFull(props) {
    if (props.upsetFullProps !== "") {
        /*let blob = new Blob([props.upsetFullProps], {type: 'image/svg+xml'});
        let url = URL.createObjectURL(blob);
        console.log(url)*/
        let url = props.upsetFullProps;
        return(
            /*<div className="upset">
                <svg 
                viewBox="0 0 1600 1000" 
                //preserveAspectRatio="xMinYMin meet" 
                style={{position: "absolute"}}
                //className="scaling-svg"
                //style={{paddingBottom: "92%"}}
                dangerouslySetInnerHTML={{__html: props.upsetFullProps}} />
            </div>*/
            <div className="upset">
                <img 
                    src={url}
                    //viewBox="0 0 1200 1600" 
                    width="100%"
                    height="auto"
                    alt=""
                />
            </div>
        )
    } else {
        return(
            <div className="upset"></div>
        )
    }
}
export default DisplayUpsetFull
