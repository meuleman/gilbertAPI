const express = require("express")
const app = express()
const port = 3001
const cors = require("cors");
const { compileFunction } = require("vm");
app.use(cors());

var sentence = '';
var spectrum = '';
var upsetFull = '';
var upsetLimited = '';
var completionNum = 0;

app.get("/", (request, response) => {
    const { spawn } = require("child_process");
    const pythonSentence = spawn('bash', ["/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_sentences_react_CT20220830.sh", request.headers.reqstate]);
    const pythonSpectrum = spawn('bash', ["/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_spectrum_react_CT20220831.sh", request.headers.reqstate]);
    const pythonUpset = spawn('bash', ["/Users/ctrader/Documents/region_similarity_pipeline/pipeline_analysis_scripts/region_comparison_pipeline_analysis_upset_react_CT20220831.sh", request.headers.reqstate]);
    
    //Sentence generation
    pythonSentence.stdout.on('data', (data) => {
        //console.log("data = " + data.toString());
        sentence += data.toString()
    });
    pythonSentence.stderr.on('data', (data) => {
        console.error(`stderr: ${data}`);
    });
    pythonSentence.on("close", () => {
        completionNum += 1;
        //console.log("sentence = " + sentence);
        //sentence = '';
    });

    //Spectrum plot generation
    pythonSpectrum.stdout.on('data', (data) => {
        //console.log("spectrum = " + data.toString());
        spectrum += data.toString()
    });
    pythonSpectrum.stderr.on('data', (data) => {
        console.error(`spectrum stderr: ${data}`);
    });
    pythonSpectrum.on("close", () => {
        completionNum += 1;
        //console.log("spectrum = " + spectrum);
        response.status(200).send(
            JSON.stringify({
                sentence: sentence,
                spectrum: spectrum,
                upsetFull: upsetFull,
                upsetLimited: upsetLimited
            })
        )
        sentence = '';
        spectrum = '';
        upsetFull = '';
        upsetLimited = '';
    });
    //response.end()
    
    //Upset plot generation
    pythonUpset.stdout.on('data', (data) => {
        //console.log("Upset = " + data.toString().split('split_figures_here')[1]);
        //upsetFull += data.toString().split('split_figures_here')[0]
        // upsetLimited += data.toString().split('split_figures_here')[1]
        console.log("upset = " + data)
        
    });
    pythonUpset.stderr.on('data', (data) => {
        console.error(`Upset stderr: ${data}`);
        //console.log("upset = " + "f\napp s")
    });
    pythonUpset.on("close", () => {
        completionNum += 1;
        //console.log("Upset = " + upsetLimited);
        
    });
})

app.listen(port, () => {
    console.log("Ready to receive response")
})