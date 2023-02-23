const express = require('express')
const app = express()
const port = 3001
const cors = require('cors')
// const { compileFunction } = require("vm");
// const { resolve } = require("path");
// const cluster = require('cluster');
// import { Buffer } from 'node:buffer';
app.use(cors())

// let sentence = ''
// let spectrum = ''
// let upsetBuffArr = []
// let upsetFull = ''
// let upsetLimited = ''
// let barBuffArr = []
// let barNormal = ''
// let barCluster = ''
// let barReduced = ''
let findRegionsBuffArr = []
let findRegionValsObject = ''

app.get('/', (request, response) => {
  const { spawn } = require('child_process')
  if (request.headers.purpose === 'getRegions') {
    // get region information
    const pythonFindRegions = spawn(
      'bash', [
        '/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch/similar_regions_CT20230222.sh', 
        request.headers.reqstate, 
        request.headers.factors
      ])
    pythonFindRegions.stdout.on('data', (data) => {
      findRegionsBuffArr.push(data)
    })
    pythonFindRegions.stderr.on('data', (data) => {
      console.error(`time: ${data}`)
    })
    pythonFindRegions.on('close', () => {
      findRegionValsObject = JSON.parse(Buffer.concat(findRegionsBuffArr))
      response.status(200).send(
        JSON.stringify({
          findRegionValsObject
        })
      )
      findRegionValsObject = ''
      findRegionsBuffArr = []
    })
  } else if (request.headers.purpose === 'functional') {
  }
})

app.listen(port, () => {
  console.log('Ready to rumble')
})
