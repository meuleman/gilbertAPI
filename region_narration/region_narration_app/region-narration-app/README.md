Setup:
1) conda create -n genomic_narration_env --file ../region-narration-flask/requirements.txt --channel default --channel bioconda --channel conda-forge 
2) conda activate genomic_narration_env
3) npm install

To run the local development version on port ####:
1) conda activate genomic_narration_env
2) PORT=#### npm run localdev
The frontend will run on port #### and the server will run on port 5435. 

By instead running:
1) conda activate genomic_narration_env
2) npm run localdev
The frontend will run on port 3000. 
