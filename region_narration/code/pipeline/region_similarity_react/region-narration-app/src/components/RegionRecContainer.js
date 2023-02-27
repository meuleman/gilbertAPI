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
    selectedFactors: [...Array(factors.length).keys()],
    possibleFactors: [...Array(factors.length).keys()],
    autoFactorSelection: true,
    included: new Array((factors.length) + 1).fill(true)
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
    console.log(location)
    this.setState({
      // genomicLocation: location, // so we get the actual queried location
      // sentence: '',
      spectrum: '',
      upsetFull: '',
      upsetLimited: '',
      barNormal: '',
      barCluster: '',
      barReduced: ''//,
      // regionalVals: ''
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
      // console.log(this.state.autoFactorSelection)
      Axios({
        method: 'GET',
        url: 'http://localhost:3001/',
        headers: {
          reqState: location,
          purpose: 'getRegions',
          factors: JSON.stringify(this.state.selectedFactors),
          autoFactorSelection: String(this.state.autoFactorSelection)
        }
      }).then((response) => {
        const data = response.data.findRegionValsObject
        const queryRegion = data[0]
        this.setState({
          genomicLocation: queryRegion.coordinates, // so we get the actual queried location
          regionalVals: data,
          sentence: queryRegion.sentence,
          selectedFactors: queryRegion.target_factors,
        }, () => {
          console.log(this.state.regionalVals)
          var newIncluded = new Array(factors.length + 1).fill(false)
          newIncluded = newIncluded.map((value, index) => 
          { if (this.state.selectedFactors.includes(index)) {
              return !value
          } else {
              return value;
          }})
          if (newIncluded.slice(0, 34).reduce((pSum, val) => pSum + val, 0) === 34) {
              newIncluded[34] = true
          }
          this.setIncluded(newIncluded)
          // return <SelectFactors preselectProps={{}}/>
        })
      })
    }
  }
  
  render () {
    return (
      <div className='container'>

        <div className='region'>
          <Header />
          <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
          <QueryRegion queryProps={this.state.genomicLocation} />
          <SelectFactors factorProps={{
            possibleFactors: this.state.possibleFactors, 
            setSelectedFactors: this.setSelectedFactors,
            data: this.state.regionalVals,
            selectedFactors: this.state.selectedFactors,
            setAutoFactorSelection: this.setAutoFactorSelection,
            setIncluded: this.setIncluded,
            included: this.state.included
          }}/>
          <DisplaySentence sentenceProps={this.state.sentence} />
          <FactorHeatmap heatmapProps={{
            data: this.state.regionalVals, 
            search: this.setGenomicLocation,
            setAutoFactorSelection: this.setAutoFactorSelection,
            autoFactorSelection: this.state.autoFactorSelection
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
