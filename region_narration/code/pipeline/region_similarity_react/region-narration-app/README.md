To run the local development version:
1) npm install (if needed)
2) conda create -n genomic_narration_env python=3.11 --no-default-packages --yes (if needed)
3) conda activate genomic_narration_env
4) conda install $(cat ../region-narration-flask/requirements.txt | cut -d= -f1 | sed 1,3d) (if needed)
5) npm run localdev

The app will run on port 3000 and the server will run on port 5435.

To change the port the app runs on:
5) PORT=#### npm run localdev
