import React from "react"
import Header from "./Header"
import InputRegion from "./InputRegion"
import DisplaySentence from "./DisplaySentence"
import DisplaySpectrum from "./DisplaySpectrum"
import DisplayUpsetFull from "./DisplayUpsetFull"
import DisplayUpsetLimited from "./DisplayUpsetLimited"
import DisplayBar from "./DisplayBar"
import Axios from 'axios';
import DisplayTabs from "./DisplayTabs"



class RegionRecContainer extends React.Component {
    state = {
        genomicLocation: "",
        sentence: "",
        spectrum: "",
        upsetFull: "",
        upsetLimited: "",
        barNormal: "",
        barCluster: "",
        barReduced: ""
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
                upsetLimited: response.data.upsetLimited,
                barNormal: response.data.barNormal,
                barCluster: response.data.barCluster,
                barReduced: response.data.barReduced
            }, () => {
            })
        });
    }
    //showSentence = 
    render() {
        return (
            <div className="container">
                <div className="region">
                    <Header />
                    <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
                    {/*<DisplaySentence sentenceProps={this.state.sentence} />*/}
                </div>
                
                {/*<div className="figures">
                    <div className="left">
                        <DisplayBar barProps={this.state.barNormal} />
                    </div>
                    <div className="right">
                        <DisplaySpectrum spectrumProps={this.state.spectrum} />
                        <DisplayUpsetFull upsetFullProps={this.state.upsetFull} />
                    </div>
                </div>*/}
                
                <div className="tab-figures">
                    <DisplayTabs tabsProps={{
                        sentence: this.state.sentence, 
                        upset: this.state.upsetFull,
                        spectrum: this.state.spectrum,
                        bar: this.state.barNormal
                    }} />
                </div>
            </div>
        )
    }
}

export default RegionRecContainer