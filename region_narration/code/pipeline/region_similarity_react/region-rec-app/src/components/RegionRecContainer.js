import React from "react"
import Header from "./Header"
import InputRegion from "./InputRegion"
import DisplaySentence from "./DisplaySentence"
import DisplaySpectrum from "./DisplaySpectrum"
import DisplayUpsetFull from "./DisplayUpsetFull"
import DisplayUpsetLimited from "./DisplayUpsetLimited"
import BarD3 from "./BarD3"
import DisplayBar from "./DisplayBar"
import Axios from 'axios';
import DisplayTabs from "./DisplayTabs"
// import rd3 from 'react-d3-library';
import sampleData from '../devmode_top_similar_regions_results_CT20221123.json';



class RegionRecContainer extends React.Component {
    devMode = true;
    state = {
        genomicLocation: "",
        sentence: "",
        spectrum: "",
        upsetFull: "",
        upsetLimited: "",
        barNormal: "",
        barCluster: "",
        barReduced: "",
        regionalVals: ""
    };
    setGenomicLocation = (location) => {
        this.setState({
            genomicLocation: location,
            sentence: "",
            spectrum: "",
            upsetFull: "",
            upsetLimited: "",
            barNormal: "",
            barCluster: "",
            barReduced: "",
            regionalVals: ""
        });
        if (this.devMode) {
            console.log(sampleData)
            this.setState({
                regionalVals: sampleData
            }, () => {
            })
        } else {
            Axios({
                method: "GET",
                url: "http://localhost:3001/",
                headers: {
                    reqState: location,
                    purpose: "getRegions"
                }
            }).then((response) => {
                this.setState({
                    regionalVals: response.data.findRegionValsObject
                }, () => {
                });
            });
            // Axios({
            //     method: "GET",
            //     url: "http://localhost:3001/",
            //     headers: {
            //         reqState: location,
            //         purpose: "functional"
            //     }
            // }).then((response) => {
            //     this.setState({
            //         sentence: response.data.sentence,
            //         spectrum: response.data.spectrum
            //         /*upsetFull: response.data.upsetFull,
            //         upsetLimited: response.data.upsetLimited,
            //         barNormal: response.data.barNormal,
            //         barCluster: response.data.barCluster,
            //         barReduced: response.data.barReduced*/
            //     }, () => {
            //     });
            // });
            // Axios({
            //     method: "GET",
            //     url: "http://localhost:3001/",
            //     headers: {
            //         reqState: location,
            //         purpose: "relations"
            //     }
            // }).then((response) => {
            //     this.setState({
            //         upsetFull: response.data.upsetFull,
            //         upsetLimited: response.data.upsetLimited
            //     }, () => {
            //     });
            // });
            // Axios({
            //     method: "GET",
            //     url: "http://localhost:3001/",
            //     headers: {
            //         reqState: location,
            //         purpose: "regions"
            //     }
            // }).then((response) => {
            //     this.setState({
            //         barNormal: response.data.barNormal,
            //         barCluster: response.data.barCluster,
            //         barReduced: response.data.barReduced
            //     }, () => {
            //     });
            // });
        }
    }
    





    render() {
        return (
            <div className="container">
                
                <div className="region">
                    <Header />
                    <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
                    <DisplaySentence sentenceProps={this.state.sentence} />
                </div>
                
                {/*<div className="figures"> //DONT USE from here
                    <div className="left">
                        <DisplayBar barProps={this.state.barNormal} />
                    </div>
                    <div className="right">
                        <DisplaySpectrum spectrumProps={this.state.spectrum} />
                        <DisplayUpsetFull upsetFullProps={this.state.upsetFull} />
                    </div>
                </div> //DONT USE to here */} 
                
                <div className="tab-figures">
                    <DisplayTabs tabsProps={{
                        //sentence: this.state.sentence, 
                        upset: this.state.upsetLimited,
                        spectrum: this.state.spectrum,
                        bar: this.state.barNormal,
                        regionalVals: this.state.regionalVals
                    }} />
                </div>
                {/*<BarD3 barProps={{regionalVals: this.state.regionalVals}} />*/}
            </div>

            
        )
    }
}

export default RegionRecContainer