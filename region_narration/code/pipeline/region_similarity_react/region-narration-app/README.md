To run the development version:
1) npm install (if needed)
2) conda create -n genomic_narration_env --file ../region-narration-flask/requirements.txt (if needed)
3) conda activate genomic_narration_env
4) npm run localdev

The app will run on port 3000 and the server will run on port 5435





To run the production version:

1) npm install (if needed)
2) PUBLIC_URL=. npm run build
3) Open the site at: https://tools.altiusinstitute.org/ng/region-narration-app/

Using PUBLIC_URL should allow the site to be served from a subdirectory of tools0. Otherwise, the site would need to be served from the top-most level of tools.altiusinstitute.org.
