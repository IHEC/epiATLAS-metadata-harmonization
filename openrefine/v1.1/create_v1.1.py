import os.path
from subprocess import run

import numpy as np
import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v1.0/IHEC_metadata_harmonization.v1.0.extended.csv'  # csv to build project from

# remove non-reprocessed epigenomes
reprocessed_metadata = pd.read_csv('./openrefine/v1.1/epiatlas_metadata.csv')
v1_0_df = pd.read_csv(initial_csv)
v1_0_df['epirr_id_without_version'] = v1_0_df['EpiRR'].str.split('\\.', 1, expand=True)[0]
v1_0_df = v1_0_df[v1_0_df['epirr_id_without_version'].isin(reprocessed_metadata['epirr_id_without_version'].unique())]

# remove pruned epigenomes due to coverage
with open('./openrefine/v1.1/pruned_cov.txt') as file:
    pruned_epirrs = [line.rstrip() for line in file]
v1_0_df = v1_0_df[~v1_0_df['epirr_id_without_version'].isin(pruned_epirrs)]

v1_1_extended_intermediate_csv = './openrefine/v1.1/IHEC_metadata_harmonization.v1.1.extended.intermediate.csv'
v1_0_df.to_csv(v1_1_extended_intermediate_csv, index=False)

# create project with intermediate version
intermediate_project_name = os.path.splitext(os.path.basename(v1_1_extended_intermediate_csv))[0]
run([openrefine_client, '--delete', intermediate_project_name], check=False)
run([openrefine_client, '--create', v1_1_extended_intermediate_csv], check=True)

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', 'openrefine/v1.1/new_col_and_fix_some_harmonized_sample_ontology_intermediate.json',
     intermediate_project_name],
    check=True)

run([openrefine_client, '--export', f'--output={v1_1_extended_intermediate_csv}', intermediate_project_name],
    check=True)

v1_1_extended_df = pd.read_csv(v1_1_extended_intermediate_csv)
v1_1_extended_df.insert(v1_1_extended_df.columns.get_loc("harmonized_donor_life_stage"),
                        "automated_harmonized_donor_age_in_years",
                        v1_1_extended_df['harmonized_donor_age'].str.split('-').apply(
                            lambda l: np.nan if l[0] == 'unknown' else np.mean([float(i.rstrip("+")) for i in l])))
v1_1_extended_df.loc[
    v1_1_extended_df['harmonized_donor_age_unit'] == 'week', 'automated_harmonized_donor_age_in_years'] /= 52

v1_1_extended_df.loc[
    v1_1_extended_df['harmonized_donor_age_unit'] == 'day', 'automated_harmonized_donor_age_in_years'] /= 365

v1_1_extended_csv = './openrefine/v1.1/IHEC_metadata_harmonization.v1.1.extended.csv'
v1_1_extended_df.to_csv(v1_1_extended_csv, index=False)

# v1_0_extended_intermediate =

f = open('./openrefine/v1.1/Routput.txt', 'w')
run(['Rscript', 'add_higher_order_v1.1.R'], check=True, cwd='./openrefine/v1.1', stdout=f)

v1_1_extended_df = pd.read_csv(v1_1_extended_csv)
v1_1_extended_df = v1_1_extended_df[["EpiRR",
                                     "project",
                                     "harmonized_biomaterial_type",
                                     "harmonized_sample_ontology_intermediate",
                                     "harmonized_sample_disease_high",
                                     "harmonized_sample_disease_intermediate",
                                     "harmonized_EpiRR_status",
                                     "harmonized_cell_type",
                                     "harmonized_cell_line",
                                     "harmonized_tissue_type",
                                     "harmonized_sample_ontology_curie",
                                     "harmonized_cell_markers",
                                     "automated_harmonized_sample_ontology",
                                     "automated_harmonized_sample_ontology_term",
                                     "sample_ontology_term_high_order_JeffreyHyacinthe",
                                     "sample_ontology_term_high_order_JonathanSteif",
                                     "automated_harmonized_sample_ontology_term_intermediate_order_unique",
                                     "automated_harmonized_sample_ontology_term_high_order_unique",
                                     "automated_harmonized_sample_ontology_term_intermediate_order",
                                     "automated_harmonized_sample_ontology_term_high_order",
                                     "harmonized_sample_disease",
                                     "harmonized_sample_disease_ontology_curie",
                                     "automated_harmonized_sample_disease_ontology_curie_ncit",
                                     "automated_harmonized_sample_disease_ontology_term_intermediate_order_unique",
                                     "automated_harmonized_sample_disease_ontology_term_high_order_unique",
                                     "automated_harmonized_sample_disease_ontology_term_intermediate_order",
                                     "automated_harmonized_sample_disease_ontology_term_high_order",
                                     "harmonized_donor_type",
                                     "harmonized_donor_id",
                                     "harmonized_donor_age",
                                     "harmonized_donor_age_unit",
                                     "automated_harmonized_donor_age_in_years",
                                     "harmonized_donor_life_stage",
                                     "harmonized_donor_sex",
                                     "harmonized_donor_health_status",
                                     "harmonized_donor_health_status_ontology_curie",
                                     "automated_harmonized_donor_health_status_ontology_curie_ncit",
                                     "automated_harmonized_donor_health_status_ontology_term_intermediate_order_unique",
                                     "automated_harmonized_donor_health_status_ontology_term_high_order_unique",
                                     "automated_harmonized_donor_health_status_ontology_term_intermediate_order",
                                     "automated_harmonized_donor_health_status_ontology_term_high_order",
                                     "epirr_id_without_version"]]

v1_1_extended_df.to_csv(v1_1_extended_csv, index=False)
final_csv = './openrefine/v1.1/IHEC_metadata_harmonization.v1.1.csv'
v1_1_extended_df.loc[:,
~(v1_1_extended_df.columns.str.startswith('automated') | v1_1_extended_df.columns.str.contains('order'))].to_csv(
    final_csv,
    index=False)

old = pd.read_csv(initial_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
old.sort_index(1, inplace=True)
new = pd.read_csv(v1_1_extended_csv)
new.index = new.EpiRR
new.sort_index(0, inplace=True)
new.sort_index(1, inplace=True)
assert new.index.isin(old.index).all()
old = old[old.index.isin(new.index)]
assert old.columns.difference(new.columns) == "harmonized_donor_life_status"
assert (new.columns.difference(old.columns) == ["automated_harmonized_donor_age_in_years",
                                                "epirr_id_without_version"]).all()
shared_cols = old.columns.intersection(new.columns)
assert shared_cols.size == 40

diff_tbl = old[shared_cols].compare(new[shared_cols])
diff_tbl.rename(columns={'self': 'v1.0', 'other': 'v1.1'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v1.1/diff_v1.0_v1.1.json', indent=True)
