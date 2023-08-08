import React from 'react'
import Header from './Header'
import Autocomplete from './Autocomplete/Autocomplete'
import DisplaySentence from './DisplaySentence'
import Axios from 'axios'
// import DisplayTabs from './DisplayTabs'
import sampleData from '../devmode_top_similar_regions_results_CT20221129.json'
// import devModeObj from '../devmode.json'
// import SelectFactors from './SelectFactors'
import TopFactorDropdown from './TopFactorDropdown'
import factors from '../factors.json'
import FactorHeatmap from './heatmap/FactorHeatmapContainer'
import QueryRegion from './QueryRegion'
import Order from './Order'
// import SimSearchHeader from './SimSearchHeader'
import SelectDetailLevel from './SelectDetailLevel'
import Spinner from './Spinner'
import Overlap from './Overlap'
import RegionsOnIdeogram from './RegionsOnIdeogram'
// import LoadDataTest from './LoadDataTest'

class RegionNarrationAppContainer extends React.Component {
  // devMode = devModeObj.devMode
  devMode = false;
  state = {
    genomicLocation: '',
    queriedGenomicLocation: '',
    regionalVals: '',
    fullData:'',
    selectedFactors: [...Array(factors.length).keys()],
    possibleFactors: [...Array(factors.length).keys()],
    queryBuildingSelectedFactors: [],
    factorOrder: [...Array(factors.length).keys()],
    searchedByFactor: false,
    autoFactorSelection: true,
    included: new Array((factors.length) + 1).fill(true),
    detailLevel: 1,
    maxDetailLevel: 34,
    spinner: false,
    devQueryHost: 'localhost',
    devQueryPort: 3001,
    prodQueryHost: 'tools.altiusinstitute.org',
    prodQueryPort: 9877,
    productionFlaskURL: "https://explore.altius.org:5001/narration",
    // developmentFlaskURL: "http://127.0.0.1:5435/narration",
    developmentFlaskURL: "http://127.0.0.1:5000/narration",
    productionQueryBuildingFlaskURL: "https://explore.altius.org:5001/search_by_factor",
    developmentQueryBuildingFlaskURL: "http://127.0.0.1:5000/search_by_factor",
    orders: [4, 5, 6, 7, 8, 9, 10, 11],
    orderSegmentLengths: [1048576, 262144, 65536, 16384, 4096, 1024, 256, 64],
    currentOrder: 4,
    searchedFactorInd: '',
    requiredRegionURLParams: ['chr', 'start', 'stop', 'order'],
    requiredFactorURLParams: ['factors', 'order'],
    regionMethod: 'hilbert_sfc'
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
        //   console.log(this.state.regionalVals);
        })

      }
    })
  } 

  setSelectedFactors = (factorInds) => {
    this.setState({
      selectedFactors: factorInds
    }, () => {
      if (this.state.regionalVals) {
        this.setGenomicLocation(this.state.regionalVals[0].region_ind);
      };
    });
  };

  setIncluded = (newIncluded) => {
    this.setState({
      included: newIncluded
    });
  };

  queryBuilding = (providedScale='', providedFactors='', popped=false) => {
    if ((this.state.queryBuildingSelectedFactors.length !== 0) || (providedFactors !== '')) {
      this.setState({
        spinner: true,
      });

      let scale = this.state.currentOrder
      if (providedScale !== '') {
        scale = providedScale
      }

      let factors = JSON.stringify(this.state.queryBuildingSelectedFactors)
      if (providedFactors !== '') {
        factors = JSON.stringify(providedFactors)
        this.setState({
          queryBuildingSelectedFactors: providedFactors
        })
      }

      // Flask Server
      const postBody = {
        factors: factors,
        scale: scale,
        numRegions: 20,
        regionMethod: this.state.regionMethod,
      };
      let flaskURL = this.state.productionQueryBuildingFlaskURL;
      if (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") {
        flaskURL = this.state.developmentQueryBuildingFlaskURL
      }
      Axios({
        method: 'POST',
        url: flaskURL,
        data: postBody
      }).then((response) => {
        this.setState({
          regionalVals: response.data,
          fullData: null,
          spinner: false,
          searchedByFactor: true,
          currentOrder: scale,
          factorOrder: response.data[0].factor_order
        }, () => {
          // push state on new query and update tab title
          this.updateURL(popped);
        })
      })
    }
  }

  setQueryBuildingSelectedFactors = (newSelectedFactors) => {
    this.setState({
      queryBuildingSelectedFactors: newSelectedFactors
    }, () => {
      this.queryBuilding()
    })
  };

  setAutoFactorSelection = (value) => {
    this.setState({
      autoFactorSelection: value
    });
  };

  // determines the order that results should be returned for
  determineBestOrder = (coordinates) => {
    const start = coordinates.start
    // subtract 1 in case provided region stop at end of segment
    // with floor function, this would map to the next segment over otherwise
    const stop = coordinates.stop - 1
    let orderInd = 0
    let baseline = 0  // for if there is no segment that 100% overlaps query region
    this.state.orderSegmentLengths.map((l, i) => {
      let startSegment = Math.floor(start / l)
      let stopSegment = Math.floor(stop / l)
      // if the smallest order cannot encompass entirety of query region, update baseline
      if (i === 0) {
        baseline = stopSegment - startSegment
      }
      // look for largest order that encompasses as much of the query region as possible
      if (startSegment === stopSegment - baseline) {
        orderInd = i
      }
    })
    return this.state.orders[orderInd]
  }
  
  
  setGenomicLocation = (inLocation, popped=false, providedOrder='', newQueryLocation=true, providedDetailLevel='', searchedFactorInd='') => {
    if (this.state.selectedFactors.length > 0) {
      let location = {}
      if (typeof inLocation === "string") {
        location = {
          chrom: inLocation.split(':')[0],
          start: parseInt(inLocation.split(':')[1].split('-')[0]),
          stop: parseInt(inLocation.split(':')[1].split('-')[1])
        }
      } else {
        location = inLocation
      }
      if (newQueryLocation) {
        this.setState({
          queriedGenomicLocation: location
        }, () => {
          // console.log(this.state.queriedGenomicLocation)
        })
      }

      // const startbp = location.start;
      // const stopbp = location.stop;
      // const queryLength = stopbp - startbp;

      let order = ''
      if (providedOrder === '') {
        order = this.determineBestOrder(location)
      } else {
        order = providedOrder
      }
        
      this.setState({
        spectrum: '',
        upsetFull: '',
        upsetLimited: '',
        barNormal: '',
        barCluster: '',
        barReduced: '',
        spinner: true,
        searchedFactorInd: searchedFactorInd
      });
      if (this.devMode) {
        const data = sampleData
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


          
        // Flask Server
        const postBody = {
          location: `${location.chrom}:${location.start}-${location.stop}`,
          factors: JSON.stringify(this.state.selectedFactors),
          roiURL: roiURL,
          scale: order,
          queryFactorThresh: 0.5,
          numSimilarRegions: 20,
          regionMethod: this.state.regionMethod,
        };
        // console.log(`flaskURL: ${this.state.productionFlaskURL}`);
        let flaskURL = this.state.productionFlaskURL;
        if (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") {
          flaskURL = this.state.developmentFlaskURL
        }
        Axios({
          method: 'POST',
          url: flaskURL,
          data: postBody
        }).then((response) => {
          const fullData = response.data;
          // console.log(fullData);
          
          let startDetailLevel = 1;
          if (providedDetailLevel === '') {
            const queryProdRanks = fullData.map((d) => {
              return d[0].prod_rank
            });
            
            let maxProdRank = 0;
            queryProdRanks.map((value, index) => {
              if (Math.max(maxProdRank, value) !== maxProdRank) {
                maxProdRank = value
                startDetailLevel = index + 1
              }
              return null
            });
          } else {
            startDetailLevel = providedDetailLevel
          }
                  
          const data = fullData[startDetailLevel - 1];
          const queryRegion = data[0];

          const returnedCoordinates = queryRegion.coordinates;

          const returnedCoordinatesObj = {
            chrom: returnedCoordinates.split(':')[0],
            start: parseInt(returnedCoordinates.split(':')[1].split('-')[0]),
            stop: parseInt(returnedCoordinates.split(':')[1].split('-')[1])
          }

          this.setState({
            genomicLocation: returnedCoordinatesObj, // so we get the actual queried/returned location
            regionalVals: data,
            fullData: fullData,
            spinner: false,
            maxDetailLevel: fullData.length,
            detailLevel: startDetailLevel,
            currentOrder: order,
            searchedByFactor: false,
          }, () => {
            // push state on new query and update tab title
            this.updateURL(popped);
          })
      })
      .catch((err) => {
        console.error(`error:     ${JSON.stringify(err)}`);
        console.error(`post body: ${JSON.stringify(postBody)}`);
        alert('Query Failed: Try another region.');
        this.setState({
          spinner: false
        })
		  });
      }
    }
  }

  updateURL = (popped) => {
    let hostname = "https://explore.altius.org/ng/"
    if (window.location.hostname === "localhost" || window.location.hostname === "127.0.0.1") {
      const localPort = window.location.port;
      hostname = "http://localhost:" + localPort + "/";
    } else if (window.location.hostname === "region-narration-app.vercel.app") {
      hostname = "https://region-narration-app.vercel.app/";
    }

    const searchedByFactor = this.state.searchedByFactor
    let newUrl = ''
    let historyTitle = "Genomic Narration - "

    if (searchedByFactor) {
      const order = this.state.currentOrder;
      const factors = this.state.queryBuildingSelectedFactors;
      newUrl = hostname + "?order=" + order + "&factors=" + factors;
      historyTitle += "Factor Search"
    } else {
      const coordinates = this.state.genomicLocation
      const chrm = coordinates.chrom
      const start = coordinates.start
      const end = coordinates.stop
      const order = this.state.currentOrder
      newUrl = hostname + "?order=" + order + "&chr=" + chrm + "&start=" + start + "&stop=" + end;
      historyTitle += coordinates;
    }

    if (!popped) {
      try {
        window.history.pushState(newUrl, historyTitle, newUrl);
      }
      catch (err) {
        console.log(`Error: ${JSON.stringify(err)}`);
	    }
    }
    document.title = historyTitle;
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
      const paramKeys = Object.keys(paramsObj)

      const regionalURL = JSON.stringify(paramKeys.sort()) === JSON.stringify(this.state.requiredRegionURLParams.sort())
      const factorURL = JSON.stringify(paramKeys.sort()) === JSON.stringify(this.state.requiredFactorURLParams.sort())

      if (regionalURL) {
        const coordinates = paramsObj.chr + ":" + paramsObj.start + "-" + paramsObj.stop
        this.setGenomicLocation(coordinates, popped)
      } else if (factorURL) {
        const order = parseInt(paramsObj.order)
        const factors = paramsObj.factors.split(',').map(i => parseInt(i))
        this.queryBuilding(order, factors, popped)
      } else if ((paramKeys.length > 0) && (window.location.href !== document.location.href.split("?")[0])) {
        window.location.href = document.location.href.split("?")[0];
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
          <div className='search-options'>
            <Autocomplete onChangeLocation={this.setGenomicLocation} />
            <TopFactorDropdown topFactorDropdownProps={{
              // scale: this.state.currentOrder,
              // search: this.setGenomicLocation,
              setQueryBuildingSelectedFactors: this.setQueryBuildingSelectedFactors,
              queryBuildingSelectedFactors: this.state.queryBuildingSelectedFactors,
              factorOrder: this.state.factorOrder
            }}/>
          </div>
        </div>
        <div className='results' id='results'>
            <QueryRegion 
              queryProps={{
                region: this.state.genomicLocation,
                scale: this.state.currentOrder,
                spinner: this.state.spinner,
                searchedByFactor: this.state.searchedByFactor,
              }} 
            />
            <Order 
              orderProps={{
                region: this.state.genomicLocation,
                searchedByFactor: this.state.searchedByFactor,
                setGenomicLocation: this.setGenomicLocation,
                queryBuilding: this.queryBuilding,
                currentScale: this.state.currentOrder,
                queriedLocation: this.state.queriedGenomicLocation,
                returnedLocation: this.state.genomicLocation,
                possibleScales: this.state.orders,
                searchedByFactor: this.state.searchedByFactor,
            }}
            />
            <Overlap overlapProps={{
              queriedLocation: this.state.queriedGenomicLocation,
              returnedLocation: this.state.genomicLocation,
              searchedByFactor: this.state.searchedByFactor,
            }}/>
            <SelectDetailLevel 
              detailProps={{
                setDetailLevel: this.setDetailLevel, 
                detailLevel: this.state.detailLevel,
                maxDetailLevel: this.state.maxDetailLevel,
                returnedLocation: this.state.genomicLocation,
                searchedByFactor: this.state.searchedByFactor,
              }}
            />
            <Spinner spinnerProps={this.state.spinner}/>
            <DisplaySentence sentenceProps={{
              // sentence: this.state.sentence, 
              data: this.state.regionalVals,
              detailLevel: this.state.detailLevel,
              searchedByFactor: this.state.searchedByFactor,
              }} />
            {/* <SimSearchHeader simSearchHeaderProps={this.state.regionalVals} /> */}
            <div className='figures'>
              <RegionsOnIdeogram ideogramProps={{
                data: this.state.regionalVals,
                searchedByFactor: this.state.searchedByFactor
              }}/>
              <FactorHeatmap heatmapProps={{
                data: this.state.regionalVals, 
                search: this.setGenomicLocation,
                setAutoFactorSelection: this.setAutoFactorSelection,
                autoFactorSelection: this.state.autoFactorSelection,
                detailLevel: this.state.detailLevel,
                scale: this.state.currentOrder,
                searchedByFactor: this.state.searchedByFactor,
              }}/>
            </div>
          </div>
        </div>

    )
  }
}

export default RegionNarrationAppContainer
