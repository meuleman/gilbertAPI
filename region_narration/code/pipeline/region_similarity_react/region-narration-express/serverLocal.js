const express = require('express')
const app = express()
const port = process.env.PORT || 3001
const cors = require('cors')
const { inherits } = require('util')
// const { compileFunction } = require("vm");
// const { resolve } = require("path");
// const cluster = require('cluster');
// import { Buffer } from 'node:buffer';
app.use(cors())

let findRegionsBuffArr = []
let findRegionValsObject = ''

app.get('/', (request, response) => {
  const { spawn } = require('child_process')
  // console.log(JSON.stringify(request.headers))
  if (request.headers.purpose === 'getRegions') {
    // get region information
    // console.log(request.headers.factors)
    const pythonFindRegions = spawn(
      'bash', [
        '/Users/ctrader/Documents/AG/NarrativeGenomics/region_narration/code/simsearch/similar_regions.sh',
        request.headers.reqstate, 
        request.headers.factors,
        request.headers.roiurl
      ], {
        stdio: 'pipe',
        shell: true,
      })
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
