const express = require('express')
const app = express()
const port = 3001
const cors = require('cors')
// const { compileFunction } = require("vm");
// const { resolve } = require("path");
// const cluster = require('cluster');
// import { Buffer } from 'node:buffer';
app.use(cors())

let sentence = ''
let spectrum = ''
let upsetBuffArr = []
let upsetFull = ''
let upsetLimited = ''
let barBuffArr = []
let barNormal = ''
let barCluster = ''
let barReduced = ''
let findRegionsBuffArr = []
let findRegionValsObject = ''

app.get('/', (request, response) => {
  const { spawn } = require('child_process')
  if (request.headers.purpose === 'getRegions') {
    // get region information
    const pythonFindRegions = spawn('bash', ['/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_top_similar_regions_CT20230110.sh', request.headers.reqstate, 'chromatin'])
    // const pythonFindRegions = spawn('bash', ['/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_top_similar_regions_CT20221128.sh', request.headers.reqstate, 'chromatin'])
    // const pythonFindRegions = spawn('bash', ['/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/pipeline/rank_based/region_comparison_top_similar_regions_CT20221215.sh', request.headers.reqstate, 'chromatin'])
    pythonFindRegions.stdout.on('data', (data) => {
      // findRegionValsObject = JSON.parse(data)
      findRegionsBuffArr.push(data)
    })
    pythonFindRegions.stderr.on('data', (data) => {
      console.error(`Bar stderr: ${data}`)
    })
    pythonFindRegions.on('close', () => {
      findRegionValsObject = JSON.parse(Buffer.concat(findRegionsBuffArr))
      // console.log(findRegionValsObject)
      // console.log(findRegionValsObject)
      response.status(200).send(
        JSON.stringify({
          findRegionValsObject
        })
      )
      findRegionValsObject = ''
      findRegionsBuffArr = []
    })
  } else if (request.headers.purpose === 'functional') {
    // const sentencePromise = new Promise((resolve, reject) => {
    //   // Sentence generation
    //   const pythonSentence = spawn('bash', ['/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_sentences_react_CT20221010.sh', request.headers.reqstate])
    //   pythonSentence.stdout.on('data', (data) => {
    //     sentence += data.toString()
    //   })
    //   pythonSentence.stderr.on('data', (data) => {
    //     console.error(`stderr: ${data}`)
    //   })
    //   pythonSentence.on('close', () => {
    //     resolve()
    //   })
    // })

  //   const spectrumPromise = new Promise((resolve, reject) => {
  //     // Spectrum plot generation
  //     const pythonSpectrum = spawn('bash', ['/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_spectrum_react_CT20220928.sh', request.headers.reqstate])
  //     pythonSpectrum.stdout.on('data', (data) => {
  //       spectrum += data.toString()
  //     })
  //     pythonSpectrum.stderr.on('data', (data) => {
  //       console.error(`spectrum stderr: ${data}`)
  //     })
  //     pythonSpectrum.on('close', () => {
  //       resolve()
  //     })
  //   })

  //   Promise.all([spectrumPromise]).then(() => { // sentencePromise, 
  //     response.status(200).send(
  //       JSON.stringify({
  //         //sentence,
  //         spectrum
  //       })
  //     )
  //     //sentence = ''
  //     spectrum = ''
  //   })
  // } else if (request.headers.purpose === 'relations') {
  //   // Upset plot generation
  //   const pythonUpset = spawn('bash', ['/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_upset_react_CT20220928.sh', request.headers.reqstate])
  //   pythonUpset.stdout.on('data', (data) => {
  //     upsetBuffArr.push(data)
  //   })
  //   pythonUpset.stderr.on('data', (data) => {
  //     console.error(`Upset stderr: ${data}`)
  //   })
  //   pythonUpset.on('close', () => {
  //     const upsetBuffObj = JSON.parse(Buffer.concat(upsetBuffArr))
  //     upsetFull += upsetBuffObj.fullUpset
  //     upsetLimited += upsetBuffObj.limitedUpset
  //     response.status(200).send(
  //       JSON.stringify({
  //         upsetFull,
  //         upsetLimited
  //       })
  //     )
  //     upsetFull = ''
  //     upsetLimited = ''
  //     upsetBuffArr = []
  //   })
  // } else if (request.headers.purpose === 'regions') {
  //   // Bar plot generation
  //   const pythonBar = spawn('bash', ['/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_bar_react_CT20220928.sh', request.headers.reqstate])
  //   pythonBar.stdout.on('data', (data) => {
  //     barBuffArr.push(data)
  //   })
  //   pythonBar.stderr.on('data', (data) => {
  //     console.error(`Bar stderr: ${data}`)
  //   })
  //   pythonBar.on('close', () => {
  //     const barBuffObj = JSON.parse(Buffer.concat(barBuffArr))
  //     barNormal += barBuffObj.normalBar
  //     barCluster += barBuffObj.clusterBar
  //     barReduced += barBuffObj.reducedBar
  //     response.status(200).send(
  //       JSON.stringify({
  //         barNormal,
  //         barCluster,
  //         barReduced
  //       })
  //     )
  //     barNormal = ''
  //     barCluster = ''
  //     barReduced = ''
  //     barBuffArr = []
  //   })
  };

  /* const upsetPromise = new Promise((resolve, reject) => {
        //Upset plot generation
        const pythonUpset = spawn('bash', ["/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_upset_react_CT20220928.sh", request.headers.reqstate]);
        pythonUpset.stdout.on('data', (data) => {
            upsetBuffArr.push(data)
        });
        pythonUpset.stderr.on('data', (data) => {
            console.error(`Upset stderr: ${data}`);
        });
        pythonUpset.on("close", () => {
            var upsetBuffObj =  JSON.parse(Buffer.concat(upsetBuffArr));
            upsetFull += upsetBuffObj.fullUpset;
            upsetLimited += upsetBuffObj.limitedUpset;
            resolve();
        });
    }); */

  /* const barPromise = new Promise((resolve, reject) => {
        //Bar plot generation
        const pythonBar = spawn('bash', ["/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_bar_react_CT20220928.sh", request.headers.reqstate]);
        pythonBar.stdout.on('data', (data) => {
            barBuffArr.push(data)
        });
        pythonBar.stderr.on('data', (data) => {
            console.error(`Bar stderr: ${data}`);
        });
        pythonBar.on("close", () => {
            var barBuffObj =  JSON.parse(Buffer.concat(barBuffArr));
            barNormal += barBuffObj.normalBar;
            barCluster += barBuffObj.clusterBar;
            barReduced += barBuffObj.reducedBar;
            response.status(200).write(
                JSON.stringify({
                    barNormal: barNormal,
                    barCluster: barCluster,
                    barReduced: barReduced
                })
            )
            resolve();
        });
    }); */
  /* Promise.all([sentencePromise, spectrumPromise, upsetPromise, barPromise]).then(() => {// upsetPromise, barPromise
        response.status(200).send(
            JSON.stringify({
                sentence: sentence,
                spectrum: spectrum,
                upsetFull: upsetFull,
                upsetLimited: upsetLimited,
                barNormal: barNormal,
                barCluster: barCluster,
                barReduced: barReduced
            })
        )
        //response.status(200).end()
        sentence = '';
        spectrum = '';
        upsetFull = '';
        upsetLimited = '';
        barNormal = '';
        barCluster = '';
        barReduced = '';
        barBuffArr = [];
        upsetBuffArr = [];
    }); */
})

app.listen(port, () => {
  console.log('Ready to rumble')
})
