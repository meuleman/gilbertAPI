To run the local development version:
1) npm install (if needed)
2) conda create -n genomic_narration_env --file ../region-narration-flask/requirements.txt  --channel default --channel bioconda (if needed)
3) conda activate genomic_narration_env
4) npm run localdev

The app will run on port 3000 and the server will run on port 5435.

To change the port the app runs on:
4) PORT=#### npm run localdev
