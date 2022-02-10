import os.path
import re
from subprocess import run

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')
mypath = './openrefine/v0.5'

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.4/IHEC_metadata_harmonization.v0.4.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]

# create project with old version
run([openrefine_client, '--create', initial_csv], check=True)

# apply rules
run([openrefine_client, '--apply', 'openrefine/v0.5/Solving_Markers_HealthStatus_CellType_GFrosi_PE_v3.json',
     old_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.5/disease_ontology_term_changes.json',
     old_project_name], check=True)

# export to intermediate version, because we need to do some cleaning in python now
intermediate_csv = "./openrefine/v0.5/IHEC_metadata_harmonization.v0.5.intermediate.csv"
intermediate_project_name = os.path.splitext(os.path.basename(intermediate_csv))[0]

run([openrefine_client, '--export', f'--output={intermediate_csv}',
     old_project_name], check=True)

# read csv to pandas
intermediate_tbl = pd.read_csv(intermediate_csv)
separator = '::'

# first clean disease ontology column
pattern = r'code=(C[0-9]*)'  # extract ncim term from url
ncim_ids = intermediate_tbl.disease_ontology_term.str.extractall(pattern)
extracted_ncim = ncim_ids.groupby(level=0).agg(lambda x: separator.join(sorted(['ncim:' + i for i in set(x)])))
intermediate_tbl.loc[extracted_ncim.index, 'disease_ontology_term'] = extracted_ncim[0]

# now clean sample ontology column
biomaterial2onto = {'primary cell': r'cl:[0-9]{7}', 'primary cell culture': r'cl:[0-9]{7}',
                    'primary tissue': r'uberon:[0-9]{7}', 'cell line': r'efo:[0-9]{7}'}
# make sure ontologies are lowercase
for index, row in intermediate_tbl.iterrows():
    if not pd.isna(row.sample_ontology_term):
        sample_ontology_split = re.findall(biomaterial2onto[row.biomaterial_type], row.sample_ontology_term, re.I)
        intermediate_tbl.loc[index, 'sample_ontology_term'] = separator.join(sample_ontology_split).lower()

intermediate_tbl.to_csv(intermediate_csv)

# run([openrefine_client, '--create', intermediate_csv], check=True)
