cleaner_score_fusion_VIZBI.ipynb
- A notebook with an example used in the 2022 VIZBI poster for score fusion between different data types. This notebook is gene-centric and does not use preset gene intervals like the other notebooks do. This notebook was used to create the original idea of score fusion / combined score with the same formula used in the STRING database: combinedscore=1∏(1−individualscore)

create_combined_score.ipynb
- This notebook is used only to combined individual data modalities 

individual_data_percentile_calculation.ipynb
- This notebook was used to test using percentiles as a universal scoring method for any data modality. 

region_recipe.ipynb
- A notebook with the intention to create a 'recipe' that can be applied to any region of any scale which provides all information known about a specific region (TF motif, DHS enrichment, epilogos chromatin states, STRING associations, etc.). This notebook is currently incomplete and only has 25kb region samples. 

STRING_subscore_analyses.ipynb
- A general guide for understanding how the STRING database was designed as well as how the combined score/associations were created. Plotting for getting fraction of combined score each individual score contributes is included.

utilizing_STRING.ipynb
- The goal of this notebook is so that given a random genomic coordinate, we can find the closest gene and give its associated genes in either assemblies (hg38 and hg19) and provide an ROIurl to visualize these genes.


All *.py files are dataframes that were created by having to run code through parallel such as getting distribution or enrichments for all 100k 25kb regions. The only exceptions are: combined_subscores.py and split_on_column_via_pandas.py

combined_subscores.py
- This is the original STRING database file indicating how they calculated subscores and how they combined scores.

split_on_column_via_pandas.py
- This python file is written by Alex to split dataframes more easily and saves them as individual bed files. This file in particular separates the DHS dataframe into separate component bed files. 