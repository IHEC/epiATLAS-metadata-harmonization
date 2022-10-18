import os.path
from subprocess import run

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.11/IHEC_metadata_harmonization.v0.11.extended.csv'  # csv to build project from

# create project with intermediate version
intermediate_project_name = os.path.splitext(os.path.basename(initial_csv))[0]
run([openrefine_client, '--delete', intermediate_project_name], check=True)
run([openrefine_client, '--create', initial_csv], check=True)

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', 'openrefine/v1.0/amed-crest_update.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v1.0/sex_tissue_donorid_AMED_CREST_updated_OpenRefine.json',
     intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v1.0/ENCODE_disease.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v1.0/additional_fixes.json', intermediate_project_name],
    check=True)

v1_0_extended_intermediate_csv = './openrefine/v1.0/IHEC_metadata_harmonization.v1.0.extended.intermediate.csv'
run([openrefine_client, '--export', f'--output={v1_0_extended_intermediate_csv}', intermediate_project_name],
    check=True)

# v1_0_extended_intermediate = pd.read_csv(v1_0_extended_intermediate_csv)

f = open('./openrefine/v1.0/Routput.txt', 'w')
run(['Rscript', 'add_higher_order_v1.0.R'], check=True, cwd='./openrefine/v1.0', stdout=f)

v1_0_extended_csv = './openrefine/v1.0/IHEC_metadata_harmonization.v1.0.extended.csv'  # csv to build project from
v1_0_extended = pd.read_csv(v1_0_extended_csv)
v1_0_extended.sort_values(by='EpiRR', inplace=True)
#
old_names = ['EpiRR', 'project', 'harm_biomaterial_type',
             'harm_sample_ontology_intermediate', 'harm_disease_high',
             'harm_disease_intermediate', 'EpiRR_status', 'harm_cell_type',
             'harm_line', 'harm_tissue_type', 'harm_sample_ontology_curie',
             'harm_markers', 'sample_ontology_term_high_order_JeffreyHyacinthe',
             'sample_ontology_term_high_order_JonathanSteif', 'harm_disease',
             'harm_disease_ontology_curie',
             'automated_harm_disease_ontology_curie_ncit', 'donor_type',
             'harm_donor_id', 'harm_donor_age', 'harm_donor_age_unit',
             'harm_donor_life_stage', 'harm_donor_sex', 'harm_donor_health_status',
             'harm_donor_health_status_ontology_curie',
             'automated_harm_donor_health_status_ontology_curie_ncit',
             'harm_donor_life_status', 'automated_sample_ontology',
             'automated_sample_ontology_term',
             'automated_harm_sample_ontology_term_intermediate_order_unique',
             'automated_harm_sample_ontology_term_high_order_unique',
             'automated_harm_sample_ontology_term_intermediate_order',
             'automated_harm_sample_ontology_term_high_order',
             'automated_harm_disease_ontology_term_intermediate_order_unique',
             'automated_harm_disease_ontology_term_high_order_unique',
             'automated_harm_disease_ontology_term_intermediate_order',
             'automated_harm_disease_ontology_term_high_order',
             'automated_harm_donor_health_status_ontology_term_intermediate_order_unique',
             'automated_harm_donor_health_status_ontology_term_high_order_unique',
             'automated_harm_donor_health_status_ontology_term_intermediate_order',
             'automated_harm_donor_health_status_ontology_term_high_order']

new_names = ['EpiRR', 'project', 'harmonized_biomaterial_type',
             'harmonized_sample_ontology_intermediate', 'harmonized_sample_disease_high',
             'harmonized_sample_disease_intermediate', 'harmonized_EpiRR_status', 'harmonized_cell_type',
             'harmonized_cell_line', 'harmonized_tissue_type', 'harmonized_sample_ontology_curie',
             'harmonized_cell_markers', 'sample_ontology_term_high_order_JeffreyHyacinthe',
             'sample_ontology_term_high_order_JonathanSteif', 'harmonized_sample_disease',
             'harmonized_sample_disease_ontology_curie',
             'automated_harmonized_sample_disease_ontology_curie_ncit',
             'harmonized_donor_type',
             'harmonized_donor_id', 'harmonized_donor_age', 'harmonized_donor_age_unit',
             'harmonized_donor_life_stage', 'harmonized_donor_sex', 'harmonized_donor_health_status',
             'harmonized_donor_health_status_ontology_curie',
             'automated_harmonized_donor_health_status_ontology_curie_ncit',
             'harmonized_donor_life_status', 'automated_harmonized_sample_ontology',
             'automated_harmonized_sample_ontology_term',
             'automated_harmonized_sample_ontology_term_intermediate_order_unique',
             'automated_harmonized_sample_ontology_term_high_order_unique',
             'automated_harmonized_sample_ontology_term_intermediate_order',
             'automated_harmonized_sample_ontology_term_high_order',
             'automated_harmonized_sample_disease_ontology_term_intermediate_order_unique',
             'automated_harmonized_sample_disease_ontology_term_high_order_unique',
             'automated_harmonized_sample_disease_ontology_term_intermediate_order',
             'automated_harmonized_sample_disease_ontology_term_high_order',
             'automated_harmonized_donor_health_status_ontology_term_intermediate_order_unique',
             'automated_harmonized_donor_health_status_ontology_term_high_order_unique',
             'automated_harmonized_donor_health_status_ontology_term_intermediate_order',
             'automated_harmonized_donor_health_status_ontology_term_high_order']

v1_0_extended.rename(columns=dict(zip(old_names, new_names)), inplace=True)
v1_0_extended = v1_0_extended[['EpiRR', 'project', 'harmonized_biomaterial_type',
                               'harmonized_sample_ontology_intermediate', 'harmonized_sample_disease_high',
                               'harmonized_sample_disease_intermediate',
                               'harmonized_EpiRR_status',
                               'harmonized_cell_type',
                               'harmonized_cell_line', 'harmonized_tissue_type', 'harmonized_sample_ontology_curie',
                               'harmonized_cell_markers', 'automated_harmonized_sample_ontology',
                               'automated_harmonized_sample_ontology_term',
                               'sample_ontology_term_high_order_JeffreyHyacinthe',
                               'sample_ontology_term_high_order_JonathanSteif',
                               'automated_harmonized_sample_ontology_term_intermediate_order_unique',
                               'automated_harmonized_sample_ontology_term_high_order_unique',
                               'automated_harmonized_sample_ontology_term_intermediate_order',
                               'automated_harmonized_sample_ontology_term_high_order',
                               'harmonized_sample_disease',
                               'harmonized_sample_disease_ontology_curie',
                               'automated_harmonized_sample_disease_ontology_curie_ncit',
                               'automated_harmonized_sample_disease_ontology_term_intermediate_order_unique',
                               'automated_harmonized_sample_disease_ontology_term_high_order_unique',
                               'automated_harmonized_sample_disease_ontology_term_intermediate_order',
                               'automated_harmonized_sample_disease_ontology_term_high_order',
                               'harmonized_donor_type',
                               'harmonized_donor_id', 'harmonized_donor_age', 'harmonized_donor_age_unit',
                               'harmonized_donor_life_stage', 'harmonized_donor_sex',
                               'harmonized_donor_health_status',
                               'harmonized_donor_health_status_ontology_curie',
                               'automated_harmonized_donor_health_status_ontology_curie_ncit',
                               'automated_harmonized_donor_health_status_ontology_term_intermediate_order_unique',
                               'automated_harmonized_donor_health_status_ontology_term_high_order_unique',
                               'automated_harmonized_donor_health_status_ontology_term_intermediate_order',
                               'automated_harmonized_donor_health_status_ontology_term_high_order',
                               'harmonized_donor_life_status']]

v1_0_extended.to_csv(v1_0_extended_csv, index=False)

final_csv = './openrefine/v1.0/IHEC_metadata_harmonization.v1.0.csv'

v1_0_extended.loc[:,
~(v1_0_extended.columns.str.startswith('automated') | v1_0_extended.columns.str.contains('order'))].to_csv(final_csv,
                                                                                                           index=False)
old2new_cols = dict(zip([
    'EpiRR',
    'project',
    'harm_biomaterial_type',
    'harm_sample_ontology_intermediate',
    'harm_disease_high',
    'harm_disease_intermediate',
    'EpiRR_status',
    'harm_cell_type',
    'harm_line',
    'harm_tissue_type',
    'harm_sample_ontology_curie',
    'harm_markers',
    'sample_ontology',
    'sample_ontology_term',
    'sample_ontology_term_high_order_JeffreyHyacinthe',
    'sample_ontology_term_high_order_JonathanSteif',
    'sample_ontology_term_intermediate_order_unique',
    'sample_ontology_term_high_order_unique',
    'sample_ontology_term_intermediate_order',
    'sample_ontology_term_high_order',
    'harm_disease',
    'harm_disease_ontology_curie',
    'disease_ontology_curie_ncit',
    'disease_ontology_term_intermediate_order_unique',
    'disease_ontology_term_high_order_unique',
    'disease_ontology_term_intermediate_order',
    'disease_ontology_term_high_order',
    'donor_type',
    'harm_donor_id',
    'harm_donor_age',
    'harm_donor_age_unit',
    'harm_donor_life_stage',
    'harm_donor_sex',
    'harm_donor_health_status',
    'harm_donor_health_status_ontology_curie',
    'donor_health_status_ontology_curie_ncit',
    'donor_health_status_ontology_term_intermediate_order_unique',
    'donor_health_status_ontology_term_high_order_unique',
    'donor_health_status_ontology_term_intermediate_order',
    'donor_health_status_ontology_term_high_order',
    'harm_donor_life_status'
],
    [
    'EpiRR',
    'project',
    'harmonized_biomaterial_type',
    'harmonized_sample_ontology_intermediate',
    'harmonized_sample_disease_high',
    'harmonized_sample_disease_intermediate',
    'harmonized_EpiRR_status',
    'harmonized_cell_type',
    'harmonized_cell_line',
    'harmonized_tissue_type',
    'harmonized_sample_ontology_curie',
    'harmonized_cell_markers',
    'automated_harmonized_sample_ontology',
    'automated_harmonized_sample_ontology_term',
    'sample_ontology_term_high_order_JeffreyHyacinthe',
    'sample_ontology_term_high_order_JonathanSteif',
    'automated_harmonized_sample_ontology_term_intermediate_order_unique',
    'automated_harmonized_sample_ontology_term_high_order_unique',
    'automated_harmonized_sample_ontology_term_intermediate_order',
    'automated_harmonized_sample_ontology_term_high_order',
    'harmonized_sample_disease',
    'harmonized_sample_disease_ontology_curie',
    'automated_harmonized_sample_disease_ontology_curie_ncit',
    'automated_harmonized_sample_disease_ontology_term_intermediate_order_unique',
    'automated_harmonized_sample_disease_ontology_term_high_order_unique',
    'automated_harmonized_sample_disease_ontology_term_intermediate_order',
    'automated_harmonized_sample_disease_ontology_term_high_order',
    'harmonized_donor_type',
    'harmonized_donor_id',
    'harmonized_donor_age',
    'harmonized_donor_age_unit',
    'harmonized_donor_life_stage',
    'harmonized_donor_sex',
    'harmonized_donor_health_status',
    'harmonized_donor_health_status_ontology_curie',
    'automated_harmonized_donor_health_status_ontology_curie_ncit',
    'automated_harmonized_donor_health_status_ontology_term_intermediate_order_unique',
    'automated_harmonized_donor_health_status_ontology_term_high_order_unique',
    'automated_harmonized_donor_health_status_ontology_term_intermediate_order',
    'automated_harmonized_donor_health_status_ontology_term_high_order',
    'harmonized_donor_life_status'
    ]))

old = pd.read_csv(initial_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
old.sort_index(1, inplace=True)
new = pd.read_csv(v1_0_extended_csv)
new.index = new.EpiRR
new.rename(columns={v: k for k, v in old2new_cols.items()}, inplace=True)
new.sort_index(0, inplace=True)
new.sort_index(1, inplace=True)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v0.11', 'other': 'v1.0'}, inplace=True)
diff_tbl.rename(columns={k: k + ':' + v for k, v in old2new_cols.items()}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v1.0/diff_v0.11_v1.0.json', indent=True)
