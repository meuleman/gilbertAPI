import React from "react";

function DisplaySpectrum(props) {
    if (props.spectrumProps !== "") {
        /*let blob = new Blob([props.spectrumProps], {type: 'image/svg+xml'});
        let url = URL.createObjectURL(blob);
        console.log(url)*/
        let url = props.spectrumProps;
        return(
            /*<div className="spectrum">
                <svg
                viewBox="0 0 1200 650"
                //preserveAspectRatio="xMinYMin meet"
                style={{position: "absolute"}}
                dangerouslySetInnerHTML={{__html: props.spectrumProps}} />
                
            </div>*/
            //<div className="container" style={{width:"80%", height:"80%", border:"1px dashed black"}}>
                //<img src=url style="width:50%"></img>
            //</div>
            <div className="spectrum">
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
            <div className="spectrum"></div>
        )
    }
}
export default DisplaySpectrum
