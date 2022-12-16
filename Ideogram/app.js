console.log("THIS IS IN NODE.JS");
const PORT_NUMBER = 3000;

// Requiring express for routing                                                                                                                                                                                        
const express = require('express');

// Creating app from express                                                                                                                                                                                            
const app = express();

// Requiring in-built file system                                                                                                                                                                                       
const fs = require('fs');
/*var fs = require('fs');
var lines = fs.readFileSync('EpilogosScores.txt').toString().split("\r\n");
var scores = [];
for(i in lines) {
  scores.push(lines[i]);
}
console.log(scores);

for (let index = 0; index < scores.length; index++) {
  var tabs = //how to split a single string in scores array by tabs to get individual characteristics
  scores.push ({
    
  })
}
console.log(scores);*/

app.use(express.static(__dirname + 'public'));
app.use('/js', express.static(__dirname, + '/js'));
app.use('/css', express.static(__dirname, + '/css'));

// GET request to the root of the app                                                                                                                                                                                   
app.get('/', function(req, res) {
  res.sendFile(__dirname + '/index.html');
})

// Creating server at port 3000 (PORT_NUMBER)                                                                                                                                                                           
app.listen(PORT_NUMBER, function(req, res) {
  console.log(`Server started on port ${PORT_NUMBER}; visit http://localhost:${PORT_NUMBER}`);
})

