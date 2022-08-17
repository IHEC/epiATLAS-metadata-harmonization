import os
import pandas as pd

# make sure the working directory when running this file is the project root of the git repository
os.chdir('../../')

old = pd.read_csv('openrefine/v0.8/IHEC_metadata_harmonization.v0.8.csv')
old.index = old.EpiRR
new = pd.read_csv('openrefine/v0.9/IHEC_metadata_harmonization.v0.9.csv')
new.index = new.EpiRR

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v0.8', 'other': 'v0.9'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v0.9/diff_v0.8_v0.9.json', indent=True)
