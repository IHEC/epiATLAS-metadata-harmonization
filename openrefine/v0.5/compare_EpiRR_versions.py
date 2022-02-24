import os

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

old = pd.read_json('merged/EpiAtlas_EpiRR_metadata_all.merged_minimal_sorted_old.json', orient='records')
old.index = old.EpiRR.str.split('.').str[0]
new = pd.read_json('merged/EpiAtlas_EpiRR_metadata_all.merged_minimal_sorted.json', orient='records')
new.index = new.EpiRR.str.split('.').str[0]

shared_ids = old.index.intersection(new.index)
old[~old.index.isin(shared_ids)].value_counts(['project', 'cell_type'])
new[~new.index.isin(shared_ids)].value_counts(['project', 'biomaterial_type'])

diff_tbl = old[old.index.isin(shared_ids)].compare(new[new.index.isin(shared_ids)])
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('diff.json', indent=True)

old_v0_5 = pd.read_csv('openrefine/v0.5/IHEC_metadata_harmonization.v0.5.intermediate_old.csv')
old_v0_5.index = old_v0_5.EpiRR.str.split('.').str[0]
old_v0_5.drop(columns=['Unnamed: 0'], inplace=True)
new_v0_5 = pd.read_csv('openrefine/v0.5/IHEC_metadata_harmonization.v0.5.intermediate.csv')
new_v0_5.index = new_v0_5.EpiRR.str.split('.').str[0]
new_v0_5.drop(columns=['Unnamed: 0'], inplace=True)

shared_ids_v0_5 = old_v0_5.index.intersection(new_v0_5.index)
old_v0_5[~old_v0_5.index.isin(shared_ids_v0_5)].value_counts(['project', 'cell_type'])
new_v0_5[~new_v0_5.index.isin(shared_ids_v0_5)].value_counts(['project', 'biomaterial_type'])

diff_tbl_v0_5 = old_v0_5[old_v0_5.index.isin(shared_ids_v0_5)].compare(new_v0_5[new_v0_5.index.isin(shared_ids_v0_5)])
diff_tbl_v0_5.apply(lambda x: [x.dropna()], axis=1).to_json('diff_v0.5.json', indent=True)

diff_blueprint = old_v0_5[old_v0_5.index.isin(shared_ids_v0_5) & (old_v0_5.project == 'BLUEPRINT')].compare(
    new_v0_5[new_v0_5.index.isin(shared_ids_v0_5) & (new_v0_5.project == 'BLUEPRINT')])
