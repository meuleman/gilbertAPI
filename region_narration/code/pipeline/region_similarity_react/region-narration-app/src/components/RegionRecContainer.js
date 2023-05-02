import React from 'react'
import Header from './Header'
import InputRegion from './InputRegion'
import DisplaySentence from './DisplaySentence'
import Axios from 'axios'
// import DisplayTabs from './DisplayTabs'
import sampleData from '../devmode_top_similar_regions_results_CT20221129.json'
// import devModeObj from '../devmode.json'
// import SelectFactors from './SelectFactors'
import factors from '../factors.json'
import FactorHeatmap from './FactorHeatmap'
import QueryRegion from './QueryRegion'
// import SimSearchHeader from './SimSearchHeader'
import SelectDetailLevel from './SelectDetailLevel'
import Spinner from './Spinner'
// import LoadDataTest from './LoadDataTest'

class RegionRecContainer extends React.Component {
  // devMode = devModeObj.devMode
  devMode = false;
  state = {
    genomicLocation: '',
    // sentence: '',
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
    maxDetailLevel: 34,
    spinner: false,
    devQueryHost: 'localhost',
    devQueryPort: 3001,
    prodQueryHost: 'tools.altiusinstitute.org',
    prodQueryPort: 9877,
    flaskURL: "http://127.0.0.1:5435/narration",
    scales: [1000, 2000, 5000, 10000, 15000, 25000, 50000, 75000, 100000],
    currentScale: 1000,
    requiredURLParams: ['chr', 'start', 'stop']
  }

  setDetailLevel = (value) => {
    this.setState({
      detailLevel: value
    }, () => {
    //   console.log(this.state.detailLevel)
      if (this.state.fullData) {
        // this.setGenomicLocation(this.state.regionalVals[0].region_ind)
        this.setState({
          regionalVals: this.state.fullData[this.state.detailLevel - 1],
        //   sentence: this.state.fullData[this.state.detailLevel - 1][0].sentence
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

  updateURL = (coordinates, popped) => {
    const chrm = coordinates.split(":")[0]
    const start = coordinates.split(":")[1].split("-")[0]
    const end = coordinates.split(":")[1].split("-")[1]

    let hostname = "https://tools.altiusinstitute.org/ng/region-narration-app/"
    if (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") {
        hostname = "http://localhost:3000/"
    } else if (window.location.hostname === "region-narration-app.vercel.app") {
        hostname = "https://region-narration-app.vercel.app/"
    }

    const newUrl = hostname + "?chr=" + chrm + "&start=" + start + "&stop=" + end
    const historyTitle = "Genomic Narration - " + coordinates
    if (!popped) {
        window.history.pushState(newUrl, historyTitle, newUrl)
    }
    document.title = historyTitle;
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
  
  
  setGenomicLocation = (location, popped=false) => {
    if (this.state.selectedFactors.length > 0) {
        console.log(location)
        const startToEnd = location.split(':')[1]
        const startbp = startToEnd.split('-')[0]
        const endbp = startToEnd.split('-')[1]
        const queryLength = endbp - startbp

        let minDistance = this.state.scales.slice(-1)[0]
        let bestScaleInd = this.state.scales.length - 1
        this.state.scales.map((s, i) => {
            if (Math.min(Math.abs(s - queryLength), minDistance) !== minDistance) {
                minDistance = Math.abs(s - queryLength)
                bestScaleInd = i
            }
            return null
        })
        const scale = this.state.scales[bestScaleInd]

        this.setState({
            spectrum: '',
            upsetFull: '',
            upsetLimited: '',
            barNormal: '',
            barCluster: '',
            barReduced: '',
            spinner: true
        })
        if (this.devMode) {
            const data = sampleData
            // const queryRegion = data[0]
            this.setState({
                regionalVals: data,
        }, () => {
        })
        } else {
            // const queryURL = (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") ? `http://${this.state.devQueryHost}:${this.state.devQueryPort}/` : `https://${this.state.prodQueryHost}:${this.state.prodQueryPort}/`
            // let queryURL = `https://${this.state.prodQueryHost}:${this.state.prodQueryPort}/`
            let roiURL = "True"
            if (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") {
                // queryURL = `http://${this.state.devQueryHost}:${this.state.devQueryPort}/`
                roiURL = "False"
            } //else {
            //     queryURL = `https://${this.state.prodQueryHost}:${this.state.prodQueryPort}/`
            //     roiURL = "True"
            // }
            //   console.log(`queryURL: ${queryURL}`)
            // console.log(`location: ${location}`)
            // console.log(`selectedFactors: ${JSON.stringify(this.state.selectedFactors)}`)

            // // Express.js server
            // Axios({
            //     method: 'GET',
            //     url: queryURL,
            //     headers: {
            //         reqState: location,
            //         purpose: 'getRegions',
            //         factors: JSON.stringify(this.state.selectedFactors),
            //         roiURL: roiURL,
            //         scale: scale / 1000
            //     }
            // }).then((response) => {
            //     const fullData = response.data.findRegionValsObject
            //     console.log(fullData)

            //     // const startDetailLevel = Math.max(fullData[0][0].percentiles.filter((rank) => {
            //     //     return rank > Math.sqrt(0.9)
            //     // }).length, 1)
                
            //     const queryProdRanks = fullData.map((d) => {
            //         return d[0].prod_rank
            //     })
            //     let startDetailLevel = 1
            //     let maxProdRank = 0
            //     queryProdRanks.map((value, index) => {
            //         if (Math.max(maxProdRank, value) !== maxProdRank) {
            //             maxProdRank = value
            //             startDetailLevel = index + 1
            //         }
            //         return null
            //     })
                
            //     const data = fullData[startDetailLevel - 1]
            //     const queryRegion = data[0]

            //     const queryCoordinates = queryRegion.coordinates
                
            //     // push state on new query and update tab title
            //     this.updateURL(queryCoordinates, popped)

            //     this.setState({
            //         genomicLocation: queryCoordinates, // so we get the actual queried location
            //         regionalVals: data,
            //         fullData: fullData,
            //         spinner: false,
            //         maxDetailLevel: fullData.length,
            //         detailLevel: startDetailLevel,
            //         currentScale: scale
            //     }, () => {
            //         console.log(this.state.regionalVals)
            //     })
            // })
            // Flask Server
            Axios({
                method: 'POST',
                url: this.state.flaskURL,
                data: {
                    location: location,
                    factors: JSON.stringify(this.state.selectedFactors),
                    roiURL: roiURL,
                    scale: scale / 1000,
                    queryFactorThresh: 0.5,
                    numSimilarRegions: 20
                }

            }).then((response) => {
                const fullData = response.data
                console.log(fullData)

                // const startDetailLevel = Math.max(fullData[0][0].percentiles.filter((rank) => {
                //     return rank > Math.sqrt(0.9)
                // }).length, 1)
                
                const queryProdRanks = fullData.map((d) => {
                    return d[0].prod_rank
                })
                let startDetailLevel = 1
                let maxProdRank = 0
                queryProdRanks.map((value, index) => {
                    if (Math.max(maxProdRank, value) !== maxProdRank) {
                        maxProdRank = value
                        startDetailLevel = index + 1
                    }
                    return null
                })
                
                const data = fullData[startDetailLevel - 1]
                const queryRegion = data[0]

                const queryCoordinates = queryRegion.coordinates
                
                // push state on new query and update tab title
                this.updateURL(queryCoordinates, popped)

                this.setState({
                    genomicLocation: queryCoordinates, // so we get the actual queried location
                    regionalVals: data,
                    fullData: fullData,
                    spinner: false,
                    maxDetailLevel: fullData.length,
                    detailLevel: startDetailLevel,
                    currentScale: scale
                }, () => {
                    console.log(this.state.regionalVals)
                })
            })
        }
    }
  }

  parseQueryParams = (popped) => {
    const parseURL = () => {
        const urlParams = window.location.search.slice(1)
        let paramsObj = {}
        urlParams.split('&').forEach((param) => {
            const keyVal = param.split("=")
            if (keyVal[0].length > 0) {
                paramsObj[keyVal[0]] = decodeURIComponent(keyVal[1])
            }
        })
        return paramsObj
    }
    const paramsObj = parseURL();
    if (paramsObj) {
        let passed  = false;
        const paramKeys = Object.keys(paramsObj)
        let numParamDifferences = 0
        if (paramKeys.length === this.state.requiredURLParams.length) {
            numParamDifferences = paramKeys.filter((param) => {
                return !this.state.requiredURLParams.includes(param)
            }).length;
        }
        if (
            (paramKeys.length !== this.state.requiredURLParams.length) ||
            (numParamDifferences > 0)
        ) {
            if (window.location.href !== document.location.href.split("?")[0]) {
                window.location.href = document.location.href.split("?")[0];
            }
        } else {
            passed = true
        }
        if (passed) {
            const coordinates = paramsObj.chr + ":" + paramsObj.start + "-" + paramsObj.stop
            this.setGenomicLocation(coordinates, popped)
        }
    }
  }

  handlePopState = (e) => {
    this.parseQueryParams(true)
  }

  componentDidMount() {
    // component is mounting twice causing two renders and 2 pop states to be added, but only in development env?
    // https://stackoverflow.com/questions/63383473/why-componentdidmount-called-two-times#:~:text=Multiple%20componentDidMount%20calls%20may%20be,help%20detect%20unexpected%20side%20effects.
    // make sure this doesn't happen in production
    this.parseQueryParams(false)
    window.addEventListener("popstate", (e) => this.handlePopState(e));
  }

  componentWillUnmount() {
    window.removeEventListener("popstate", null);
  }
  
  
  render () {
    return (
      <div className='container'>
        <div className='region'>
          <Header />
          <InputRegion setGenomicLocationProps={this.setGenomicLocation} />
          <QueryRegion queryProps={{
            region: this.state.genomicLocation,
            scale: this.state.currentScale
        }} />
          <SelectDetailLevel detailProps={{
            setDetailLevel: this.setDetailLevel, 
            detailLevel: this.state.detailLevel,
            maxDetailLevel: this.state.maxDetailLevel
          }}/>
          {/* <SelectFactors factorProps={{
            possibleFactors: this.state.possibleFactors, 
            setSelectedFactors: this.setSelectedFactors,
            data: this.state.regionalVals,
            selectedFactors: this.state.selectedFactors,
            setAutoFactorSelection: this.setAutoFactorSelection,
            setIncluded: this.setIncluded,
            included: this.state.included
          }}/> */}
          <Spinner spinnerProps={this.state.spinner}/>
          <DisplaySentence sentenceProps={{
            // sentence: this.state.sentence, 
            data: this.state.regionalVals,
            detailLevel: this.state.detailLevel
            }} />
          {/* <SimSearchHeader simSearchHeaderProps={this.state.regionalVals} /> */}
          <FactorHeatmap heatmapProps={{
            data: this.state.regionalVals, 
            search: this.setGenomicLocation,
            setAutoFactorSelection: this.setAutoFactorSelection,
            autoFactorSelection: this.state.autoFactorSelection,
            detailLevel: this.state.detailLevel
          }} />
          {/* <LoadDataTest /> */}
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