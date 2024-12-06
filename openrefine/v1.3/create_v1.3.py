import os.path
from subprocess import run

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
v1_2_csv = './openrefine/v1.2/IHEC_metadata_harmonization.v1.2.extended.csv'  # csv to build project from

# create project from v1.2
intermediate_project_name = os.path.splitext(os.path.basename(v1_2_csv))[0]
run([openrefine_client, '--delete', intermediate_project_name], check=False)
run([openrefine_client, '--create', v1_2_csv], check=True)

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script
run([openrefine_client, '--apply', './openrefine/v1.3/life_stage_fixes.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', './openrefine/v1.3/fig1_high_order_fixes.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', './openrefine/v1.3/sample_ontology_intermediate_fixes.json', intermediate_project_name],
    check=True)



v1_3_intermediate_csv = './openrefine/v1.3/IHEC_metadata_harmonization.v1.3.extended.intermediate.csv'
run([openrefine_client, '--export', f'--output={v1_3_intermediate_csv}', intermediate_project_name],
    check=True)

v1_3_df_intermediate = pd.read_csv(v1_3_intermediate_csv)

# add columns that are not corrected by EpiClass
v1_2_intermediate_csv = './openrefine/v1.2/IHEC_metadata_harmonization.v1.2.extended.intermediate.csv'
v1_2_df_intermediate = pd.read_csv(v1_2_intermediate_csv)
assert v1_3_df_intermediate["EpiRR"].equals(v1_2_df_intermediate["EpiRR"])
for col in ['harmonized_donor_life_stage', 'harmonized_donor_sex']:
    v1_3_df_intermediate.insert(v1_3_df_intermediate.columns.get_loc(col) + 1, col + "_uncorrected",
                                v1_2_df_intermediate[col])

# load diff between v1.1 and v1.2
import json

diff_v1_1_v1_2 = json.load(open('./openrefine/v1.2/diff_v1.1_v1.2.json', 'r'))

# read v1.1/pruned_cov.txt into a list in one line
with open('./openrefine/v1.1/pruned_cov.txt') as file:
    pruned_epirrs = [line.rstrip() for line in file]

# get the epirr ids where harmonized_donor_life_stage_uncorrected is not equal to harmonized_donor_life_stage from v1_3_df_intermediate
epirr_ids_life_stage_diff = v1_3_df_intermediate.loc[
    v1_3_df_intermediate["harmonized_donor_life_stage_uncorrected"] != v1_3_df_intermediate[
        "harmonized_donor_life_stage"],
    "epirr_id_without_version"
]
epirr_ids_life_stage_json_diff = [entry.split('.')[0] for entry, diff in diff_v1_1_v1_2.items() if
                                  "('harmonized_donor_life_stage', 'v1.1')" in diff[0].keys()]
reverted_ids_life_stage = [entry['engineConfig']['facets'][0]['query'] for entry in
                           json.load(open('./openrefine/v1.3/life_stage_fixes.json', 'r'))]

# assert that the epirr_ids_life_stage_diff without the pruned_epirrs are the same as the epirr_ids_life_stage_json_diff without the reverted ids
assert set(epirr_ids_life_stage_diff).difference(pruned_epirrs) == set(epirr_ids_life_stage_json_diff).difference(
    set(reverted_ids_life_stage))

# get the epirr ids where harmonized_donor_sex_uncorrected is not equal to harmonized_donor_sex from v1_3_df_intermediate
epirr_ids_sex_diff = v1_3_df_intermediate.loc[
    v1_3_df_intermediate["harmonized_donor_sex_uncorrected"] != v1_3_df_intermediate["harmonized_donor_sex"],
    "epirr_id_without_version"
]
epirr_ids_sex_json_diff = [entry.split('.')[0] for entry, diff in diff_v1_1_v1_2.items() if
                           "('harmonized_donor_sex', 'v1.1')" in diff[0].keys()]

# assert that the epirr_ids_sex_diff without the pruned_epirrs are the same as the epirr_ids_sex_json_diff
assert set(epirr_ids_sex_diff).difference(pruned_epirrs) == set(epirr_ids_sex_json_diff)

### Add automated experiments from experiment metadata
reprocessed_metadata = pd.read_csv('./openrefine/v1.1/epiatlas_metadata.csv')
# make reprocessed_metadata a wide table with epirr_id_without_version as index and assay_type and experiment_type as columns
reprocessed_metadata_wide = reprocessed_metadata.pivot(index='epirr_id_without_version',
                                                       columns=['assay_type', 'experiment_type'],
                                                       values='uuid')

# merge the to column indexes into one
reprocessed_metadata_wide.columns = ['_'.join(col) for col in reprocessed_metadata_wide.columns]
# reorder reprocessed_metadata_wide columns
reprocessed_metadata_wide = reprocessed_metadata_wide[['ChIP-Seq_H3K27ac', 'ChIP-Seq_H3K27me3', 'ChIP-Seq_H3K36me3',
                                                       'ChIP-Seq_H3K4me1', 'ChIP-Seq_H3K4me3', 'ChIP-Seq_H3K9me3',
                                                       'WGBS_standard', 'WGBS_PBAT',
                                                       'RNA-Seq_mRNA-Seq', 'RNA-Seq_total-RNA-Seq']]

# first, all of the epigenomes are set to partial
reprocessed_metadata_wide['epiATLAS_status'] = 'Partial'
# now, if there is at least one histone mark AND RNA-Seq, this means there is imputed data for the other marks (also for WGBS)
# i.e., the status is set to full_imputed
chip_n = reprocessed_metadata_wide.filter(like='ChIP-Seq', axis=1).notnull().sum(axis=1)
have_rna = reprocessed_metadata_wide.filter(like='RNA-Seq', axis=1).notnull().sum(axis=1) >= 1
have_wgbs = reprocessed_metadata_wide.filter(like='WGBS', axis=1).notnull().sum(axis=1) >= 1
reprocessed_metadata_wide.loc[(chip_n >= 1) & have_rna, 'epiATLAS_status'] = 'Complete_imputed'
reprocessed_metadata_wide.loc[(chip_n == 6) & have_rna & have_wgbs, 'epiATLAS_status'] = 'Complete'

# load imputed metadata file
imputed_metadata = pd.read_csv('./openrefine/v1.2/epiatlas_imputed_metadata.csv')
imputed_metadata.data_file_path = 'imputed'
# remove "." and suffix from epirr_id column
imputed_metadata['epirr_id_without_version'] = imputed_metadata['epirr_id'].str.split('\\.', 1, expand=True)[0]
# rename DNAMethyl to WGBS in column "mark"
imputed_metadata.loc[imputed_metadata['mark'] == 'DNAMethyl', 'mark'] = 'WGBS_standard'
# add prefix ChIP-Seq for histone marks
imputed_metadata.loc[imputed_metadata['mark'] != 'WGBS_standard', 'mark'] = 'ChIP-Seq_' + imputed_metadata.loc[
    imputed_metadata['mark'] != 'DNAMethyl', 'mark']
# pivot table of imputed_metadata
imputed_metadata_wide = imputed_metadata.pivot(index='epirr_id_without_version', columns='mark',
                                               values='data_file_path')
reprocessed_metadata_wide.fillna(imputed_metadata_wide, inplace=True)

# make index to column in reprocessed_metadata_wide
reprocessed_metadata_wide.reset_index(inplace=True)

fully_imputed = (
        (reprocessed_metadata_wide.filter(like='ChIP-Seq', axis=1).notnull().sum(axis=1) == 6) &
        (reprocessed_metadata_wide.filter(like='RNA-Seq', axis=1).notnull().sum(axis=1) >= 1) &
        (reprocessed_metadata_wide.filter(like='WGBS', axis=1).notnull().sum(axis=1) >= 1))

assert reprocessed_metadata_wide.loc[
    reprocessed_metadata_wide['epiATLAS_status'].isin(
        ['Complete', 'Complete_imputed']), 'epirr_id_without_version'].equals(
    reprocessed_metadata_wide.loc[fully_imputed, 'epirr_id_without_version'])

assert reprocessed_metadata_wide.loc[
    reprocessed_metadata_wide['epiATLAS_status'] == 'Partial', 'epirr_id_without_version'].equals(
    reprocessed_metadata_wide.loc[~fully_imputed, 'epirr_id_without_version'])

# merge reprocessed_metadata_wide with v1_2_df_intermediate, validate one_to_one
v1_3_df_intermediate = v1_3_df_intermediate.merge(
    reprocessed_metadata_wide.rename(columns={'epiATLAS_status': 'epiATLAS_status_new'}), 'left',
    validate='one_to_one')

assert v1_3_df_intermediate.epiATLAS_status.equals(v1_3_df_intermediate.epiATLAS_status_new)
v1_3_df_intermediate.drop(columns='epiATLAS_status_new', inplace=True)

# double check that the histone columns are the same
for hm in ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']:
    # make sure all columns ending with hm have the same values
    assert v1_3_df_intermediate[f'automated_experiments_{hm}'].equals(v1_3_df_intermediate[f'ChIP-Seq_{hm}'])
# drop all columns starting with automated_experiments_
v1_3_df_intermediate.drop(columns=v1_3_df_intermediate.filter(like='automated_experiments_').columns, inplace=True)
# add prefix automated_experiments_ to columns
v1_3_df_intermediate.rename(columns={col: f'automated_experiments_{col}' for col in reprocessed_metadata_wide.drop(
    columns=['epirr_id_without_version', 'epiATLAS_status']).columns}, inplace=True)

### add colors to harmonized_sample_ontology_intermediate
# load the json with the colors
color_json = json.load(open('./openrefine/v1.3/IHEC_EpiATLAS_IA_colors_Apl01_2024.json', 'r'))
# get the entry where the key is 'harmonized_sample_ontology_intermediate'
intermediate_colors = {}
for color_dict in color_json:
    if color_dict.get('harmonized_sample_ontology_intermediate'):
        intermediate_colors.update(color_dict['harmonized_sample_ontology_intermediate'][0])
assert intermediate_colors
# add a column called harmonized_sample_ontology_intermediate_color to v1_3_df_intermediate right after harmonized_sample_ontology_intermediate
v1_3_df_intermediate.insert(v1_3_df_intermediate.columns.get_loc('harmonized_sample_ontology_intermediate') + 1,
                            'harmonized_sample_ontology_intermediate_color',
                            v1_3_df_intermediate['harmonized_sample_ontology_intermediate'].map(intermediate_colors))

# get the entry where the key is 'harmonized_sample_ontology_intermediate'
fig1_colors = {}
for color_dict in color_json:
    if color_dict.get('fig1_ontology_intermediate_merged'):
        fig1_colors.update(color_dict['fig1_ontology_intermediate_merged'][0])
assert fig1_colors
# add a column called harmonized_sample_ontology_intermediate_color to v1_3_df_intermediate right after harmonized_sample_ontology_intermediate
v1_3_df_intermediate.insert(v1_3_df_intermediate.columns.get_loc('harmonized_sample_ontology_term_high_order_fig1') + 1,
                            'harmonized_sample_ontology_term_high_order_fig1_color',
                            v1_3_df_intermediate['harmonized_sample_ontology_term_high_order_fig1'].map(fig1_colors))

# Read the file with the Manual_Groups
manual_groups = pd.read_csv('./openrefine/v1.3/EpiAtlas_Manual_Groups+mod_HSOI.tsv', sep='\t')

# check that order of epirr_id_without_version is the same in manual_groups and v1_3_df_intermediate
manual_groups = manual_groups.set_index('epirr_id_without_version').loc[v1_3_df_intermediate['epirr_id_without_version']].reset_index()
assert manual_groups['epirr_id_without_version'].equals(v1_3_df_intermediate['epirr_id_without_version'])

# insert Final_Manual_Group_Label from manual_groups into v1_3_df_intermediate after harmonized_sample_ontology_intermediate_color
v1_3_df_intermediate.insert(v1_3_df_intermediate.columns.get_loc('harmonized_biomaterial_type') + 1,
                            'harmonized_sample_label',
                            manual_groups['Final_Manual_Group_Label'])

### Write to csv
v1_3_extended_csv = './openrefine/v1.3/IHEC_metadata_harmonization.v1.3.extended.csv'
v1_3_df_intermediate.to_csv(v1_3_extended_csv, index=False)

# for this step,
with open('./openrefine/v1.3/Routput.txt', 'w') as f:
    run(['conda', 'run', '-n', 'v1_2_r_env', 'Rscript', 'add_higher_order_v1.3.R'], check=True, cwd='./openrefine/v1.3',
        stdout=f)


v1_3_extended_df = pd.read_csv(v1_3_extended_csv)

v1_3_extended_df = v1_3_extended_df[["EpiRR",
                                     "project",
                                     "harmonized_biomaterial_type",
                                     "harmonized_sample_label",
                                     "harmonized_sample_ontology_intermediate",
                                     "harmonized_sample_ontology_intermediate_color",
                                     "harmonized_sample_disease_high",
                                     "harmonized_sample_disease_intermediate",
                                     "harmonized_EpiRR_status",
                                     "epiATLAS_status",
                                     "harmonized_cell_type",
                                     "harmonized_cell_line",
                                     "harmonized_tissue_type",
                                     "harmonized_sample_ontology_curie",
                                     "harmonized_cell_markers",
                                     "automated_harmonized_sample_ontology",
                                     "automated_harmonized_sample_ontology_term",
                                     "harmonized_sample_ontology_term_high_order_fig1",
                                     "harmonized_sample_ontology_term_high_order_fig1_color",
                                     "harmonized_sample_organ_system_order_AnetaMikulasova",
                                     "harmonized_sample_organ_order_AnetaMikulasova",
                                     "harmonized_sample_organ_part_or_lineage_order_AnetaMikulasova",
                                     "harmonized_sample_cell_order_AnetaMikulasova",
                                     "harmonized_sample_cell_2_order_AnetaMikulasova",
                                     "harmonized_sample_cell_3_order_AnetaMikulasova",
                                     "harmonized_sample_cancer_type_order_AnetaMikulasova",
                                     "harmonized_sample_cancer_subtype_order_AnetaMikulasova",
                                     "harmonized_sample_disease",
                                     "harmonized_sample_disease_ontology_curie",
                                     "automated_harmonized_sample_disease_ontology_curie_ncit",
                                     "harmonized_donor_type",
                                     "harmonized_donor_id",
                                     "harmonized_donor_age",
                                     "harmonized_donor_age_unit",
                                     "automated_harmonized_donor_age_in_years",
                                     "harmonized_donor_life_stage",
                                     "harmonized_donor_life_stage_uncorrected",
                                     "harmonized_donor_sex",
                                     "harmonized_donor_sex_uncorrected",
                                     "harmonized_donor_health_status",
                                     "harmonized_donor_health_status_ontology_curie",
                                     "automated_harmonized_donor_health_status_ontology_curie_ncit",
                                     'automated_experiments_ChIP-Seq_H3K27ac',
                                     'automated_experiments_ChIP-Seq_H3K27me3',
                                     'automated_experiments_ChIP-Seq_H3K36me3',
                                     'automated_experiments_ChIP-Seq_H3K4me1',
                                     'automated_experiments_ChIP-Seq_H3K4me3',
                                     'automated_experiments_ChIP-Seq_H3K9me3',
                                     'automated_experiments_WGBS_standard',
                                     'automated_experiments_WGBS_PBAT',
                                     'automated_experiments_RNA-Seq_mRNA-Seq',
                                     'automated_experiments_RNA-Seq_total-RNA-Seq',
                                     "epirr_id_without_version"]]

# order the rows by: harmonized_sample_ontology_term_high_order_fig1, harmonized_sample_ontology_intermediate, harmonized_sample_label, harmonized_sample_disease_high, harmonized_sample_disease_intermediate, harmonized_donor_sex, automated_harmonized_donor_age_in_years, EpiRR
# Define the columns to sort by
sort_columns = [
    "harmonized_sample_ontology_term_high_order_fig1",
    "harmonized_sample_ontology_intermediate",
    "harmonized_sample_label",
    "harmonized_sample_disease_high",
    "harmonized_sample_disease_intermediate",
    "harmonized_donor_sex",
    "automated_harmonized_donor_age_in_years",
    "EpiRR"
]

# Sort the DataFrame
v1_3_extended_df.sort_values(
    by=sort_columns,
    key=lambda col: col if col.name == "automated_harmonized_donor_age_in_years" else col.str.casefold(),
    inplace=True
)
v1_3_extended_df.to_csv(v1_3_extended_csv, index=False)

v1_3_csv = './openrefine/v1.3/IHEC_metadata_harmonization.v1.3.csv'
v1_3_extended_df.loc[:,
~(v1_3_extended_df.columns.str.startswith('automated')
  | v1_3_extended_df.columns.str.endswith('uncorrected')
  | v1_3_extended_df.columns.str.endswith('color')
  | v1_3_extended_df.columns.str.contains('order'))].to_csv(
    v1_3_csv, index=False)

### Compare v1.2 and v1.3

old = pd.read_csv(v1_2_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
new = pd.read_csv(v1_3_extended_csv)
new.index = new.EpiRR
new.sort_index(0, inplace=True)
assert old.index.equals(new.index)

old.columns = old.columns.str.replace('H3K', 'ChIP-Seq_H3K')
# subset old for columns that do not contain WGBS or RNA-Seq
old.drop(columns=old.columns[old.columns.str.contains('WGBS') | old.columns.str.contains('RNA-Seq')], inplace=True)
old.sort_index(1, inplace=True)
new.drop(columns=new.columns[
    new.columns.str.contains('WGBS') | new.columns.str.contains('RNA-Seq') | new.columns.str.endswith('_uncorrected') |
                 new.columns.str.endswith('_color') | new.columns.str.endswith('harmonized_sample_label')],
         inplace=True)
new.sort_index(1, inplace=True)
assert old.columns.equals(new.columns)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v1.2', 'other': 'v1.3'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v1.3/diff_v1.2_v1.3.json', indent=True)
# print column in markdown format

na_props = pd.DataFrame.from_dict({col: {'Column': col, 'Examples': '', 'Explanation': '',
                                         '# Not Null (%)': f'{v1_3_extended_df[col].notnull().sum()} ({v1_3_extended_df[col].notnull().mean() * 100:.1f}%)'}
                                   for col in v1_3_extended_df.columns}, orient='index')
na_props.to_markdown('openrefine/v1.3/README.md', index=False, tablefmt='github', mode='a')
