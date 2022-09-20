import os.path
from subprocess import run

import numpy as np
import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.9/IHEC_metadata_harmonization.v0.9.extended.csv'  # csv to build project from

v10_intermediate_tbl = pd.read_csv(initial_csv)

disease_higher_tbl = pd.read_csv(
    './openrefine/v0.10/internal--IHEC_metadata_harmonization.v0.9.extended - Comp_health_disease.csv',
    usecols=['Sample_disease_high_level', 'Sample_disease_intermediate_level', 'disease', 'donor_health_status'])
disease_higher_tbl.replace('-blank-', np.nan, inplace=True)

disease_higher_tbl.rename(columns={'Sample_disease_high_level': 'disease_high_order_manual',
                                   'Sample_disease_intermediate_level': 'disease_intermediate_order_manual'},
                          inplace=True)

v10_merged = pd.merge(v10_intermediate_tbl, disease_higher_tbl.drop_duplicates(), 'outer',
                      on=['disease', 'donor_health_status'], validate='many_to_one')

assert (len(v10_intermediate_tbl) == len(v10_merged))

v10_merged.rename(columns={'age': 'donor_age'}, inplace=True)

v10_intermediate_csv = './openrefine/v0.10/IHEC_metadata_harmonization.v0.10.intermediate.csv'
v10_merged.to_csv(v10_intermediate_csv, index=False)

# create project with intermediate version
run([openrefine_client, '--create', v10_intermediate_csv], check=True)
intermediate_project_name = os.path.splitext(os.path.basename(v10_intermediate_csv))[0]

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', 'openrefine/v0.10/biomaterial_type_ENCODE_Roadmap.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v0.10/minor_disease_issues.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v0.10/sample_ontology_high-level_manual.json',
     intermediate_project_name],
    check=True)

v10_extended_csv = './openrefine/v0.10/IHEC_metadata_harmonization.v0.10.extended.csv'
run([openrefine_client, '--export', f'--output={v10_extended_csv}', intermediate_project_name], check=True)

v10_extended = pd.read_csv(v10_extended_csv)
final_csv = './openrefine/v0.10/IHEC_metadata_harmonization.v0.10.csv'
v10_extended[['EpiRR', 'EpiRR_status', 'project', 'biomaterial_type', 'cell_type', 'line', 'tissue_type',
              'sample_ontology_curie', 'sample_ontology_term_high_order_manual', 'markers', 'disease',
              'disease_ontology_curie', 'disease_high_order_manual', 'disease_intermediate_order_manual', 'donor_id',
              'donor_age', 'donor_age_unit', 'donor_life_stage', 'sex', 'donor_health_status',
              'donor_health_status_ontology_curie', 'health_state']].to_csv(final_csv, index=False)

old = pd.read_csv(initial_csv)
old.rename(columns={'age': 'donor_age'}, inplace=True)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
old.sort_index(1, inplace=True)
new = pd.read_csv(v10_extended_csv)
new.index = new.EpiRR
new.drop(columns=['disease_high_order_manual', 'disease_intermediate_order_manual'], inplace=True)
new.sort_index(0, inplace=True)
new.sort_index(1, inplace=True)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v0.9', 'other': 'v0.10'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v0.10/diff_v0.9_v0.10.json', indent=True)
