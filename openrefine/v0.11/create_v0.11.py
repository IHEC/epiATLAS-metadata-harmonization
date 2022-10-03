import os.path
import urllib
from subprocess import run

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.10/IHEC_metadata_harmonization.v0.10.extended.csv'  # csv to build project from

# create project with intermediate version
intermediate_project_name = os.path.splitext(os.path.basename(initial_csv))[0]
# run([openrefine_client, '--delete', intermediate_project_name], check=True)
run([openrefine_client, '--create', initial_csv], check=True)

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', 'openrefine/v0.11/fixing_inconsistencies.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v0.11/fix_other_and_blank.json', intermediate_project_name],
    check=True)

v0_11_extended_intermediate_csv = './openrefine/v0.11/IHEC_metadata_harmonization.v0.11.extended.intermediate.csv'
run([openrefine_client, '--export', f'--output={v0_11_extended_intermediate_csv}', intermediate_project_name],
    check=True)

v0_11_extended_intermediate = pd.read_csv(v0_11_extended_intermediate_csv)
epirr_all = urllib.request.urlopen('https://www.ebi.ac.uk/vg/epirr/view/all?format=json').read()
epirr_all_dt = pd.read_json(epirr_all)
epirr_all_dt.rename(columns={'type': 'donor_type'}, inplace=True)
v0_11_extended_intermediate_merged = pd.merge(
    v0_11_extended_intermediate,
    epirr_all_dt[['full_accession', 'donor_type']],
    how="inner",
    left_on='EpiRR',
    right_on='full_accession',
    validate='one_to_one'
)
v0_11_extended_intermediate_merged.drop(columns='full_accession', inplace=True)

assert (len(v0_11_extended_intermediate) == len(v0_11_extended_intermediate_merged))
automatic_higher_level = v0_11_extended_intermediate_merged.columns.str.endswith(
    'order') | v0_11_extended_intermediate_merged.columns.str.endswith('unique')
v0_11_extended_intermediate_merged.loc[:, ~automatic_higher_level].to_csv(v0_11_extended_intermediate_csv, index=False)

run(['Rscript', 'add_higher_order_v0.11.R'], check=True, cwd='./openrefine/v0.11')

v0_11_extended_csv = './openrefine/v0.11/IHEC_metadata_harmonization.v0.11.extended.csv'  # csv to build project from
v0_11_extended = pd.read_csv(v0_11_extended_csv)
v0_11_extended.sort_values(by='EpiRR', inplace=True)

renaming_dict = {'EpiRR': 'EpiRR',
                 'EpiRR_status': 'EpiRR_status',
                 'project': 'project',
                 'biomaterial_type': 'harm_biomaterial_type',
                 'line': 'harm_line',
                 'markers': 'harm_markers',
                 'cell_type': 'harm_cell_type',
                 'tissue_type': 'harm_tissue_type',
                 'sample_ontology_curie': 'harm_sample_ontology_curie',
                 'disease': 'harm_disease',
                 'disease_ontology_curie': 'harm_disease_ontology_curie',
                 'donor_age': 'harm_donor_age',
                 'donor_age_unit': 'harm_donor_age_unit',
                 'donor_health_status': 'harm_donor_health_status',
                 'donor_health_status_ontology_curie': 'harm_donor_health_status_ontology_curie',
                 'donor_id': 'harm_donor_id',
                 'donor_life_stage': 'harm_donor_life_stage',
                 'health_state': 'harm_donor_life_status',
                 'sex': 'harm_donor_sex',
                 'sample_ontology_term_high_order_manual': 'harm_sample_ontology_intermediate',
                 'disease_high_order_manual': 'harm_disease_high',
                 'disease_intermediate_order_manual': 'harm_disease_intermediate',
                 'donor_type': 'donor_type'}

v0_11_extended.rename(columns=renaming_dict, inplace=True)
v0_11_extended = v0_11_extended[['EpiRR', 'project', 'harm_biomaterial_type', 'harm_sample_ontology_intermediate', 'harm_disease_high', 'harm_disease_intermediate',
     'EpiRR_status', 'harm_cell_type', 'harm_line', 'harm_tissue_type', 'harm_sample_ontology_curie', 'harm_markers',
     'sample_ontology', 'sample_ontology_term', 'sample_ontology_term_high_order_JeffreyHyacinthe', 'sample_ontology_term_high_order_JonathanSteif', 'sample_ontology_term_intermediate_order_unique', 'sample_ontology_term_high_order_unique', 'sample_ontology_term_intermediate_order', 'sample_ontology_term_high_order',
     'harm_disease', 'harm_disease_ontology_curie', 'disease_ontology_curie_ncit', 'disease_ontology_term_intermediate_order_unique', 'disease_ontology_term_high_order_unique', 'disease_ontology_term_intermediate_order', 'disease_ontology_term_high_order',
     'donor_type', 'harm_donor_id', 'harm_donor_age', 'harm_donor_age_unit', 'harm_donor_life_stage', 'harm_donor_sex',
     'harm_donor_health_status', 'harm_donor_health_status_ontology_curie', 'donor_health_status_ontology_curie_ncit', 'donor_health_status_ontology_term_intermediate_order_unique', 'donor_health_status_ontology_term_high_order_unique', 'donor_health_status_ontology_term_intermediate_order', 'donor_health_status_ontology_term_high_order',
     'harm_donor_life_status']]

v0_11_extended.to_csv(v0_11_extended_csv, index=False)

final_csv = './openrefine/v0.11/IHEC_metadata_harmonization.v0.11.csv'

v0_11_extended.loc[:, v0_11_extended.columns.isin(renaming_dict.values())].to_csv(final_csv, index=False)

old = pd.read_csv(initial_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
old.sort_index(1, inplace=True)
new = pd.read_csv(v0_11_extended_csv)
new.index = new.EpiRR
new.rename(columns={v: k for k, v in renaming_dict.items()}, inplace=True)
new.drop(columns=['donor_type'], inplace=True)
new.sort_index(0, inplace=True)
new.sort_index(1, inplace=True)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v0.10', 'other': 'v0.11'}, inplace=True)
diff_tbl.rename(columns={k: k + ':' + v for k, v in renaming_dict.items()}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v0.11/diff_v0.10_v0.11.json', indent=True)
