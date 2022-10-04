import React from "react"
import Header from "./Header"
import InputRegion from "./InputRegion"
import DisplaySentence from "./DisplaySentence"
import DisplaySpectrum from "./DisplaySpectrum"
import DisplayUpsetFull from "./DisplayUpsetFull"
import DisplayUpsetLimited from "./DisplayUpsetLimited"
import Axios from 'axios';

class RegionRecContainer extends React.Component {
    state = {
        genomicLocation: "",
        sentence: "",
        spectrum: "",
        upsetFull: "",
        upsetLimited: ""
    };
    setGenomicLocation = (location) => {
        this.setState({
            genomicLocation: location
        });
        Axios({
            method: "GET",
            url: "http://localhost:3001/",
            headers: {
                reqState: location
            }
        }).then((response) => {
            this.setState({
                sentence: response.data.sentence,
                spectrum: response.data.spectrum,
                upsetFull: response.data.upsetFull,
                upsetLimited: response.data.upsetLimited
            }, () => {
                //console.log(this.state.sentence)
                //console.log(this.state.upsetLimited)
            })
            //console.log(response.data.sentence);
        });
    }
    //showSentence = 
    render() {
        return (
            <div className="container">
                <div className="inner">
                    <Header />
                    <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
                    <DisplaySentence sentenceProps={this.state.sentence} />
                </div>
                <div className="full-width">
                    <DisplaySpectrum spectrumProps={this.state.spectrum} />
                    {/*<DisplayUpsetFull upsetFullProps={this.state.upsetFull} />*/}
                    <DisplayUpsetLimited upsetLimitedProps={this.state.upsetLimited} />
                </div>
            </div>
        )
    }
}

export default RegionRecContainer