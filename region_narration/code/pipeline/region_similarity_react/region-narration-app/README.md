To run the development version:
1) npm install (if needed)
2) conda create -n genomic_narration_env python=3.8 --no-default-packages --yes (if needed)
3) conda activate genomic_narration_env
4) conda install $(cat ../region-narration-flask/requirements.txt | cut -d= -f1 | sed 1,3d) (if needed)
5) npm run localdev

The app will run on port 3000 and the server will run on port 5435





To run the production version:

1) npm install (if needed)
2) PUBLIC_URL=. npm run build
3) Open the site at: https://tools.altiusinstitute.org/ng/region-narration-app/

Using PUBLIC_URL should allow the site to be served from a subdirectory of tools0. Otherwise, the site would need to be served from the top-most level of tools.altiusinstitute.org.
