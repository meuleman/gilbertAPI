const express = require('express')
const https = require('https')
const fs = require('fs')
const normalizePort = require('normalize-port')
const nocache = require('nocache')
const morgan = require('morgan')

/**
 * Specify port for running service
 */

const app = module.exports = express();
const port = normalizePort(process.env.PORT || 3001)
app.set("port", port)

/**
 * Allow overriding of locations of similar regions script and other parameters
 */

const similarRegionsShellScript = process.env.SIMILAR_REGIONS_SH || '/home/ctrader/AG/NarrativeGenomics/region_narration/code/simsearch/similar_regions.sh'
const alternateDataDir = process.env.ALTERNATE_DATA_DIR || '/home/ctrader/AG/NarrativeGenomics/region_narration/data'
const alternateSimsearchDir = process.env.ALTERNATE_SIMSEARCH_DIR || '/home/ctrader/AG/NarrativeGenomics/region_narration/code/simsearch'

// const cors = require('cors')
// const { inherits } = require('util')
// const { compileFunction } = require("vm");
// const { resolve } = require("path");
// const cluster = require('cluster');
// import { Buffer } from 'node:buffer';
// app.use(cors())

let findRegionsBuffArr = []
let findRegionValsObject = ''

const privateKeyFn = process.env.SSLPRIVATEKEY || "/var/www/auth/altiusinstitute.org.key"
const certificateFn = process.env.SSLCERTIFICATE || "/var/www/auth/altiusinstitute.org.withIntermediateChain.crt"

const privateKey = fs.readFileSync(privateKeyFn)
const certificate = fs.readFileSync(certificateFn)

const options = {
  key: privateKey,
  cert: certificate,
}
const server = https.createServer(options, app)
server.listen(port)
server.on("error", onError)
server.on("listening", onListening)

/**
 * Event listener for HTTP server "error" event.
 */

function onError(error) {
  if (error.syscall !== "listen") {
    throw error;
  }

  const bind = typeof port === "string"
    ? "Pipe " + port
    : "Port " + port;

  // handle specific listen errors with friendly messages
  switch (error.code) {
    case "EACCES":
      console.error(`${bind} requires elevated privileges`);
      process.exit(1);
      break;
    case "EADDRINUSE":
      console.error(`${bind} is already in use`);
      process.exit(1);
      break;
    default:
      throw error;
  }
}

/**
 * Event listener for HTTP server "listening" event.
 */

function onListening() {
  const addr = server.address();
  const bind = typeof addr === "string"
    ? "pipe " + addr
    : "port " + addr.port;
  console.log(`Listening on ${bind}`);
}

/**
 * Allow CORS
 */

function cors(req, res, next) {
  res.set("Access-Control-Allow-Origin", "*");
  res.set("Access-Control-Allow-Methods", "GET");
  res.set("Access-Control-Allow-Headers", "X-Requested-With, Content-Type, reqstate, purpose, factors, roiurl");
  res.set("Access-Control-Allow-Credentials", true);
  return next();
}

/**
 * Response, CORS, cache policy and logging
 */

app.use(cors);
app.use(nocache());
app.use(morgan("combined"));

app.get("/favicon.ico", (req, res) => {
  res.sendStatus(404);
});

app.get('/', (request, response) => {
  const { spawn } = require('child_process')
  // console.log(JSON.stringify(request.headers))
  if (request.headers.purpose === 'getRegions') {
    // get region information
    // console.log(request.headers.factors)
    const similarRegionsCall = [
      similarRegionsShellScript,
      request.headers.reqstate, 
      request.headers.factors,
      request.headers.roiurl,
      alternateDataDir,
      alternateSimsearchDir,
    ]
    const pythonFindRegions = spawn(
      'bash', 
      similarRegionsCall, 
      {
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
