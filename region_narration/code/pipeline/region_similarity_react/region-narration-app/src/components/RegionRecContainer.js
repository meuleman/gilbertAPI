import React from 'react'
import Header from './Header'
import InputRegion from './InputRegion'
import DisplaySentence from './DisplaySentence'
import Axios from 'axios'
import DisplayTabs from './DisplayTabs'
import sampleData from '../devmode_top_similar_regions_results_CT20221129.json'
import devModeObj from '../devmode.json'
import SelectFactors from './SelectFactors'
import factors from '../factors.json'
import FactorHeatmap from './FactorHeatmap'
import QueryRegion from './QueryRegion'
import SimSearchHeader from './SimSearchHeader'
import SelectDetailLevel from './SelectDetailLevel'
import Spinner from './Spinner'

class RegionRecContainer extends React.Component {
  // devMode = devModeObj.devMode
  devMode = false;
  state = {
    genomicLocation: '',
    sentence: '',
    spectrum: '',
    upsetFull: '',
    upsetLimited: '',
    barNormal: '',
    barCluster: '',
    barReduced: '',
    regionalVals: '',
    fullData:'',
    selectedFactors: [...Array(factors.length).keys()],
    possibleFactors: [...Array(factors.length).keys()],
    autoFactorSelection: true,
    included: new Array((factors.length) + 1).fill(true),
    detailLevel: 1,
    spinner: false,
    devQueryHost: 'localhost',
    devQueryPort: 3001,
    prodQueryHost: 'tools.altiusinstitute.org',
    prodQueryPort: 9877,
  }

  setDetailLevel = (value) => {
    this.setState({
      detailLevel: value
    }, () => {
      console.log(this.state.detailLevel)
      if (this.state.fullData) {
        // this.setGenomicLocation(this.state.regionalVals[0].region_ind)
        this.setState({
          regionalVals: this.state.fullData[this.state.detailLevel - 1],
          sentence: this.state.fullData[this.state.detailLevel - 1][0].sentence
        }, () => {
        //   console.log(this.state.regionalVals)
        })

      }
    })
  } 

  setSelectedFactors = (factorInds) => {
    this.setState({
      selectedFactors: factorInds
    }, () => {
      if (this.state.regionalVals) {
        this.setGenomicLocation(this.state.regionalVals[0].region_ind)
      }
    })
  }

  
  setIncluded = (newIncluded) => {
    this.setState({
      included: newIncluded
    })
  }

  // setAutoFactorSelection = (value) => {
  //   console.log(value)
  //   // console.log(this.state.selectedFactors)
  //   const autoStatePromise = new Promise((resolve) => {
  //     var returnValue = this.setState({
  //       autoFactorSelection: value
  //     }, () => {
  //       return this.state.autoFactorSelection
  //     })
  //     if (returnValue === value) {
  //       resolve()
  //     }
  //   })
  //   Promise.all([autoStatePromise]).then(() => {
  //     return this.returnValue
  //   })
  // }

  setAutoFactorSelection = (value) => {
    console.log(value)
    // console.log(this.state.selectedFactors)
    this.setState({
      autoFactorSelection: value
    })
  }
  
  
  setGenomicLocation = (location) => {
    if (this.state.selectedFactors.length > 0) {
        console.log(location)
        this.setState({
        // genomicLocation: location, // so we get the actual queried location
        // sentence: '',
        spectrum: '',
        upsetFull: '',
        upsetLimited: '',
        barNormal: '',
        barCluster: '',
        barReduced: '',
        // regionalVals: ''
        spinner: true
        })
        if (this.devMode) {
        // console.log(sampleData)
        const data = sampleData
        const queryRegion = data[0]
        this.setState({
            regionalVals: data,
            sentence: queryRegion.sentence
        }, () => {
        })
        } else {
        // console.log(this.state.detailLevel)
        //   console.log(window.location.hostname)
        const queryURL = (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") ? `http://${this.state.devQueryHost}:${this.state.devQueryPort}/` : `https://${this.state.prodQueryHost}:${this.state.prodQueryPort}/`
        //   console.log(`queryURL: ${queryURL}`)
        // console.log(`location: ${location}`)
        // console.log(`selectedFactors: ${JSON.stringify(this.state.selectedFactors)}`)
        Axios({
            method: 'GET',
            url: queryURL,
            headers: {
            reqState: location,
            purpose: 'getRegions',
            factors: JSON.stringify(this.state.selectedFactors),
            }
        }).then((response) => {
            const fullData = response.data.findRegionValsObject
            const data = fullData[this.state.detailLevel - 1]
            const queryRegion = data[0]
            this.setState({
            genomicLocation: queryRegion.coordinates, // so we get the actual queried location
            regionalVals: data,
            fullData: fullData,
            sentence: queryRegion.sentence,
            spinner: false
            }, () => {
            console.log(this.state.regionalVals)
            //   console.log(this.state.fullData)
            })
        })
        }
    }
  }
  
  render () {
    return (
      <div className='container'>
        <div className='region'>
          <Header />
          <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
          <QueryRegion queryProps={this.state.genomicLocation} />
          <SelectDetailLevel detailProps={{
            setDetailLevel: this.setDetailLevel, 
            detailLevel: this.state.detailLevel,
            selectedFactors: this.state.selectedFactors
          }}/>
          <SelectFactors factorProps={{
            possibleFactors: this.state.possibleFactors, 
            setSelectedFactors: this.setSelectedFactors,
            data: this.state.regionalVals,
            selectedFactors: this.state.selectedFactors,
            setAutoFactorSelection: this.setAutoFactorSelection,
            setIncluded: this.setIncluded,
            included: this.state.included
          }}/>
          <Spinner spinnerProps={this.state.spinner}/>
          <DisplaySentence sentenceProps={this.state.sentence} />
          {/* <SimSearchHeader simSearchHeaderProps={this.state.regionalVals} /> */}
          <FactorHeatmap heatmapProps={{
            data: this.state.regionalVals, 
            search: this.setGenomicLocation,
            setAutoFactorSelection: this.setAutoFactorSelection,
            autoFactorSelection: this.state.autoFactorSelection,
            detailLevel: this.state.detailLevel
          }} />
        </div>

        {/*<div className='tab-figures'>
          <DisplayTabs tabsProps={{
            // sentence: this.state.sentence,
            upset: this.state.upsetLimited,
            spectrum: this.state.spectrum,
            bar: this.state.barNormal,
            regionalVals: this.state.regionalVals
          }}
        />
        </div>*/}
      </div>

    )
  }
}

export default RegionRecContainer