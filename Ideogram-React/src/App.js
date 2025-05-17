/* Functionality: File upload -> plotted on ideogram (Target gene shown in dropdown + related genes plotted)
   Need to be added: click through ideogram to change target gene

   Notes for future iterations to improve upon:
   1. Find ideogram endpoints when element on dropdown is pressed so that dropdown is created faster with larger files
   2. Epilogos State Score integration
   3. Gene names added next to annotations (Space could be an issue)
   4. Manual input option along with file upload option
   5. Add time feature (how long it takes for dropdown to be populated)
   6. Custom annotation options (rectanlges/circles instead of triangles)

*/

import React, { Component, useState } from "react";
import Ideogram from 'ideogram';
import axios from "axios";

import * as Constants from "./Constants.js";
import './App.css';

class App extends Component {
  
  constructor(props) {
    super(props);

    this.state = {
      genes: {},
      targetGeneNames: [],
      targetGeneSelectKey: 0,
      selectedGene: null,
      geneAnnotations: [],
      ideoKey: 0,
      ideoOrientation: 'vertical',
      ideoAssembly: 'GRCh38',
      ideoAnnotations: [],
    };

    this.uploadInputRef = React.createRef();
  }

  componentDidMount() {}

  componentWillUnmount() {}

  handleUploadInputChange() {
    const self = this;

    async function queryAnnotationService(genesObj) {
      const newTargetGeneSelectKey = self.state.targetGeneSelectKey + 1;
      const newIdeoKey = self.state.ideoKey + 1;

      const targetGenes = Object.keys(genesObj);
      let allTargetedGenesSet = new Set();
      targetGenes.forEach((targetGene) => {
        const targetedGenes = genesObj[targetGene];
        targetedGenes.forEach((targetedGene) => {
          allTargetedGenesSet.add(targetedGene);
        })
      });
      const allTargetedGenes = Array.from(allTargetedGenesSet);

      // helper function for generating service URLs
      const annotationServiceUrlForGeneNameAndAssembly = (geneName, assembly) => {
        return Constants.annotationService.scheme 
          + "://" 
          + Constants.annotationService.host 
          + ":" 
          + Constants.annotationService.port 
          + "/sets?q=" 
          + geneName 
          + "&assembly=" 
          + assembly;
      }
      
      try {
        const annotationServicePromiseResponses = await Promise.allSettled(allTargetedGenes.map(geneName => {
          const annotationServiceUrl = annotationServiceUrlForGeneNameAndAssembly(geneName, Constants.assemblyToOlderName[self.state.ideoAssembly]);
          return axios.get(annotationServiceUrl)
            .then((res) => {
              if (res.data.hits) {
                const desiredHitObjects = res.data.hits[geneName];
                const desiredHit = desiredHitObjects.filter(hit => hit.name === geneName);
                return desiredHit;
              }
            });
        }));
        
        /**
         The annotationServicePromiseResponses variable is an array of resolved 
         Promises. Resolved here means either "fulfilled" or "rejected".
         
         If the Promise is fulfilled, it will have a value key that we can use to 
         build the "geneAnnotations" state object -- the gene name, and its 
         corresponding position on the genome.
         
         Note that not every Promise will get fulfilled, as we might not get an 
         annotation object back from the service for a gene name that has no 
         corresponding entry in the annotation database. The gene name might be
         misspelled, or there just isn't anything in our database for it, which
         can happen. 
         
         In that case, we don't add it to our list of selectable genes. We just
         skip over it and consider the rest.
         */

        const fulfilledPromises = annotationServicePromiseResponses.filter(r => r.status === 'fulfilled');

        const newGenes = fulfilledPromises.map(r => {
          const annotation = r.value[0];
          return annotation.name;
        });

        const newSelectedTargetGene = targetGenes[0];

        const newGeneAnnotations = fulfilledPromises.map(r => {
          const annotation = r.value[0];
          const ideoFormattedAnnotation = {
            name: annotation.name,
            chr: annotation.chrom.replace('chr', ''), // strip 'chr' prefix from chromosome name, if present
            start: annotation.start,
            stop: annotation.stop,
          }
          return ideoFormattedAnnotation;
        });

        /**
         Once we have built the "geneAnnotations" array, we can simply pick the
         first object in this array and set that to the gene labeled in the 
         ideogram. Note that this should be an array, as this is required by the 
         ideogram library.
         */

        const targetedGenesForSelectedTargetGene = genesObj[newSelectedTargetGene];

        const newIdeoAnnotations = newGeneAnnotations.filter((ideoFormattedAnnotation) => {
          const annotationName = ideoFormattedAnnotation.name;
          if (targetedGenesForSelectedTargetGene.includes(annotationName)) {
            return ideoFormattedAnnotation;
          }
        });

        /**
         Finally, we update application state.
         */

        const newState = {
          genes: genesObj,
          targetGeneNames: targetGenes,
          targetGeneSelectKey: newTargetGeneSelectKey,
          selectedGene: newSelectedTargetGene,
          geneAnnotations: newGeneAnnotations,
          ideoKey: newIdeoKey,
          ideoAnnotations: newIdeoAnnotations,
        };

        self.setState(newState);
      }
      catch (err) {
        /**
         Something could go wrong; our annotation service might be down, or something
         else. Reporting the error can help with debugging. But we should still update 
         the state of the web application and the ideogram with some default values,
         basically "resetting" the ideogram.
         */
        console.log(`Error: ${err}`);
        const defaultState = {
          genes: {},
          targetGeneNames: [],
          targetGeneSelectKey: newTargetGeneSelectKey,
          selectedGene: null,
          geneAnnotations: [],
          ideoKey: newIdeoKey,
          ideoAnnotations: [],
        };
        self.setState(defaultState);
      }
    }

    const fileInput = document.getElementById('uploadInput');
    const reader = new FileReader();
    reader.onload = (e) => {
      const data = e.target.result;
      const newGeneLines = data.split('\n').filter((e) => e);
      const genesObj = {};
      
      newGeneLines.forEach(line => {
        const elements = line.split('\t');
        const targetGene = elements[0];
        const relatedGenesStr = elements[1];
        const relatedGenes = relatedGenesStr.split(/[,\s]+/);
        genesObj[targetGene] = relatedGenes;
      });
      
      /** 
       Once genes are read in, we query our annotations service to turn
       those gene names into positions on the genome, so that we can build 
       an "ideoAnnotations" object-array that we can pass to the AppIdeogram
       component. We also update the "ideoKey" property, which will update
       the AppIdeogram component upon change-of-state.
      */

      queryAnnotationService(genesObj);
    }
    for (let file of fileInput.files) {
      reader.readAsText(file);
    }
  }

  handleGeneSelectChange(newTargetGene) {
    /**
     All we have to do is pick the item from the "geneAnnotations" array with
     the matching name, and use that to change "selectedGene" and the "ideoAnnotations"
     object-array. We also update "ideoKey" to force the ideogram to redraw.
     */
    const newSelectedTargetGene = newTargetGene;
    const newIdeoKey = this.state.ideoKey + 1;
    const genesObj = this.state.genes;
    const targetedGenesForSelectedTargetGene = genesObj[newSelectedTargetGene];
    const geneAnnotations = this.state.geneAnnotations;
    const newIdeoAnnotations = geneAnnotations.filter((ideoFormattedAnnotation) => {
      const annotationName = ideoFormattedAnnotation.name;
      if (targetedGenesForSelectedTargetGene.includes(annotationName)) {
        return ideoFormattedAnnotation;
      }
    });
    const newStateChanges = {
      selectedGene: newSelectedTargetGene,
      ideoKey: newIdeoKey,
      ideoAnnotations: newIdeoAnnotations,
    };
    const newState = {...this.state, ...newStateChanges};
    this.setState(newState, () => {
      document.activeElement.blur(); // remove the blue outline around the select drop-down menu
    });
  }

  render() {
    return (
      <div className="box">
        <div className="row header">
          <div>
            Ideogram Prototype | React
          </div>
        </div>
        <div className="row content">
          <form name="uploadForm">
            <div>
              <input 
                id="uploadInput" 
                type="file" 
                ref={(component) => this.uploadInputRef = component}
               
                onChange={(e) => this.handleUploadInputChange(e.target.result)} />
            </div>
          </form>
          <hr />
          {
            (this.state.targetGeneNames.length > 0) 
              ? 
              <select key={this.state.targetGeneSelectKey} onChange={(e) => this.handleGeneSelectChange(e.target.value)}>
                {this.state.targetGeneNames.map((targetGeneName) => {
                  return (
                    <option value={targetGeneName} key={`gene-${targetGeneName}`}>{targetGeneName}</option>
                  )
                })}
              </select>
              : 
              <div className='noGenesLabel'>
                Please choose a text file to add its genes to a pull-down menu
              </div>
          }
          {
            (this.state.selectedGene) 
              ?
              <div>
                <hr />
                <div className='geneLabelParent'>
                  Selected/Target gene: <span className='geneLabel'>{this.state.selectedGene}</span>
                </div>
              </div>
              : 
              <div />
          }
          <hr />
          <AppIdeogram
            key={this.state.ideoKey}
            assembly={this.state.ideoAssembly}
            orientation={this.state.ideoOrientation}
            annotations={this.state.ideoAnnotations} />
        </div>
        <div className="row footer">
          {this.state.ideoAnnotations.length} genes found
        </div>
      </div>
    );
  }
}

class AppIdeogram extends Component {
  
  componentDidMount() {
    return new Ideogram({
      organism: 'human',
      assembly: this.props.assembly,
      container: '#ideo-container',
      orientation: this.props.orientation,
      annotations: this.props.annotations,
    });
  }

  render() {
    return (
      <div id="ideo-container" className="ideogramContainer"></div>
    );
  }
}

export default App;
