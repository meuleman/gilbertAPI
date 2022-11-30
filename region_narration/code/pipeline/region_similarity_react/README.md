Steps for setting up proof-of-concept / development environment:
1) Clone https://github.com/meuleman/NarrativeGenomics.git
2) Make sure you have npm installed. I think the js packages needed will be pulled through cloning the repo, but if not, also install d3, Axios, and react-tabs.
3) Ensure that in NarrativeGenomics/region_narration/code/pipeline/region_similarity_react/region-rec-app/src/components/RegionRecContainer.js, devMode is set to true.
4) Go to NarrativeGenomics/region_narration/code/pipeline/region_similarity_react/region-narration-app/ and run npm start
5) Open your browser and enter any genomic coordinate into the search bar (ie chr1:198636000-198636001 or chr2:1-2). The results for chr1:198636000-198636001 will display.

The app will run on port 3000
