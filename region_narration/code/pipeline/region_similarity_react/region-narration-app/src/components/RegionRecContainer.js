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
    possibleFactors: [...Array(factors.length).keys()]
  }

  setSelectedFactors = (factorInds) => {
    this.setState({
      selectedFactors: factorInds
    }, () => {
      console.log(this.state.selectedFactors)
    })
  }

  setGenomicLocation = (location) => {
    this.setState({
      genomicLocation: location,
      sentence: '',
      spectrum: '',
      upsetFull: '',
      upsetLimited: '',
      barNormal: '',
      barCluster: '',
      barReduced: '',
      regionalVals: ''
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
      Axios({
        method: 'GET',
        url: 'http://localhost:3001/',
        headers: {
          reqState: location,
          purpose: 'getRegions',
          factors: JSON.stringify(this.state.selectedFactors)
        }
      }).then((response) => {
        const data = response.data.findRegionValsObject
        const queryRegion = data[0]
        this.setState({
          regionalVals: data,
          sentence: queryRegion.sentence
        }, () => {
          console.log(this.state.regionalVals)
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
          <SelectFactors factorProps={[this.state.possibleFactors, this.setSelectedFactors]}/>
          <DisplaySentence sentenceProps={this.state.sentence} />
          
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
