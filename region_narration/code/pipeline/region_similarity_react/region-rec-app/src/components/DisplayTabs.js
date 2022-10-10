import React from 'react';
import { Tab, Tabs, TabList, TabPanel } from 'react-tabs';
import 'react-tabs/style/react-tabs.css';
import DisplaySentence from "./DisplaySentence"
import DisplaySpectrum from "./DisplaySpectrum"
import DisplayUpsetFull from "./DisplayUpsetFull"
import DisplayBar from "./DisplayBar"

class DisplayTabs extends React.Component {
    render() {
        //const sentence = this.props.tabsProps.sentence;
        const upset = this.props.tabsProps.upset;
        const spectrum = this.props.tabsProps.spectrum;
        const bar = this.props.tabsProps.bar;
        
        /*let upsetFigure;
        if (upset !== "") {
            let upsetBlob = new Blob([upset], {type: 'image/svg+xml'});
            let upsetURL = URL.createObjectURL(upsetBlob);
            upsetFigure = 
            <div className="upset">
                <img src={upsetURL} width="100%" height="auto" />
            </div>;
        } else {
            upsetFigure = <div></div>;
        }

        let spectrumFigure;
        if (spectrum !== "") {
            let spectrumBlob = new Blob([spectrum], {type: 'image/svg+xml'});
            let spectrumURL = URL.createObjectURL(spectrumBlob);
            spectrumFigure = 
            <div className="spectrum">
                <img src={spectrumURL} width="100%" height="auto" />
            </div>;
        } else {
            spectrumFigure = <div></div>;
        }

        let barFigure;
        if (bar !== "") {
            let barBlob = new Blob([bar], {type: 'image/svg+xml'});
            let barURL = URL.createObjectURL(barBlob);
            barFigure = 
            <div className="bar">
                <img src={barURL} width="100%" height="auto" />
            </div>;
        } else {
            barFigure = <div></div>;
        }*/
        
        let upsetURL = "";
        if (spectrum !== "") {
            let upsetBlob = new Blob([upset], {type: 'image/svg+xml'});
            upsetURL = URL.createObjectURL(upsetBlob);
        } 

        let spectrumURL = "";
        if (spectrum !== "") {
            let spectrumBlob = new Blob([spectrum], {type: 'image/svg+xml'});
            spectrumURL = URL.createObjectURL(spectrumBlob);
        } 

        let barURL = "";
        if (bar !== "") {
            let barBlob = new Blob([bar], {type: 'image/svg+xml'});
            barURL = URL.createObjectURL(barBlob);
        } 

        
        return(
            <Tabs>
                <TabList>
                    <Tab>Functional</Tab>
                    <Tab>Relations Between Regions</Tab>
                    <Tab>Individual Regions</Tab>
                </TabList>
                <TabPanel>
                    {/*<DisplaySentence sentenceProps={this.props.tabsProps.sentence} />*/}
                    <DisplaySpectrum spectrumProps={spectrumURL} />
                </TabPanel>
                <TabPanel>
                    <DisplayUpsetFull upsetFullProps={upsetURL} />
                </TabPanel>
                <TabPanel>
                    <DisplayBar barProps={barURL} />
                </TabPanel>
            </Tabs>
        )
    }
}

export default DisplayTabs