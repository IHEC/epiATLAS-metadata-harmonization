import os.path
from subprocess import run

import numpy as np
import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
v1_0_csv = './openrefine/v1.0/IHEC_metadata_harmonization.v1.0.extended.csv'  # csv to build project from

# remove non-reprocessed epigenomes
reprocessed_metadata = pd.read_csv('./openrefine/v1.1/epiatlas_metadata.csv')
v1_0_df = pd.read_csv(v1_0_csv)
v1_0_df['epirr_id_without_version'] = v1_0_df['EpiRR'].str.split('\\.', 1, expand=True)[0]
v1_0_df = v1_0_df[v1_0_df['epirr_id_without_version'].isin(reprocessed_metadata['epirr_id_without_version'].unique())]

# remove pruned epigenomes due to coverage
with open('./openrefine/v1.1/pruned_cov.txt') as file:
    pruned_epirrs = [line.rstrip() for line in file]
pruned_df = v1_0_df[v1_0_df['epirr_id_without_version'].isin(pruned_epirrs)]

v1_1_csv = './openrefine/v1.1/IHEC_metadata_harmonization.v1.1.extended.csv'
v1_1_df = pd.read_csv(v1_1_csv)
assert pruned_df.columns.difference(v1_1_df.columns) == pd.Index(['harmonized_donor_life_status'], dtype='object')
pruned_df = pruned_df.drop(columns='harmonized_donor_life_status')
pruned_df.insert(pruned_df.columns.get_loc("harmonized_donor_life_stage"),
                 "automated_harmonized_donor_age_in_years",
                 pruned_df['harmonized_donor_age'].str.split('-').apply(
                     lambda l: np.nan if l[0] == 'unknown' else np.mean([float(i.rstrip("+")) for i in l])))
pruned_df.loc[
    pruned_df['harmonized_donor_age_unit'] == 'week', 'automated_harmonized_donor_age_in_years'] /= 52

pruned_df.loc[
    pruned_df['harmonized_donor_age_unit'] == 'day', 'automated_harmonized_donor_age_in_years'] /= 365

v1_2_df_intermediate = pd.concat([v1_1_df, pruned_df])

aneta_assignments = pd.read_excel(
    './openrefine/v1.2/IHEC_metadata_harmonization.v1.1.extended_editAMv2_shared231107.xlsx')
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    aneta_assignments[['EpiRR'] + list(aneta_assignments.columns.difference(v1_2_df_intermediate.columns))], 'left',
    validate='one_to_one')

martin_tcell = pd.read_excel('./openrefine/v1.2/IHEC_TCell_Groups_MH_09302023.xlsx', sheet_name="Tcell_Grouped").dropna(
    how='all')
martin_tcell.rename(columns={'Group Label': 'sample_ontology_term_high_order_MartinHirst'}, inplace=True)
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    martin_tcell[['EpiRR'] + list(martin_tcell.columns.difference(v1_2_df_intermediate.columns))], 'left',
    validate='one_to_one')

martin_blood_cancer = pd.read_excel('./openrefine/v1.2/IHEC_Blood_Cancer_Groups_MH_09302023.xlsx',
                                    sheet_name="Groupings").dropna(how='all')
martin_blood_cancer.rename(columns={'Groupings': 'sample_ontology_term_high_order_MartinHirst'}, inplace=True)
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    martin_blood_cancer[['EpiRR'] + list(martin_blood_cancer.columns.difference(v1_2_df_intermediate.columns))], 'left',
    validate='one_to_one')

assert v1_2_df_intermediate.shape[0] == pruned_df.shape[0] + v1_1_df.shape[0]

v1_2_intermediate_csv = './openrefine/v1.2/IHEC_metadata_harmonization.v1.2.extended.intermediate.csv'

v1_2_df_intermediate.to_csv(v1_2_intermediate_csv, index=False)
# create project with intermediate version
intermediate_project_name = os.path.splitext(os.path.basename(v1_2_intermediate_csv))[0]
run([openrefine_client, '--delete', intermediate_project_name], check=False)
run([openrefine_client, '--create', v1_2_intermediate_csv], check=True)

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', './openrefine/v1.2/aneta_comments.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--export', f'--output={v1_2_intermediate_csv}', intermediate_project_name],
    check=True)

v1_2_df_intermediate = pd.read_csv(v1_2_intermediate_csv)

deselect_columns = ['CANCER_SUBTYPE', 'CANCER_TYPE', 'CELL',
                    'CELL_2', 'CELL_3', 'CL_ID', 'COMMENT', 'DONOR_DISEASE_CATEGORY',
                    'DONOR_DISEASE_TYPE', 'DONOR_HEALTH_STATUS', 'DONOR_ID', 'DONOR_ID_N',
                    'DONOR_RELATIVES_N', 'DONOR_SAMPLE_N', 'DONOR_SAMPLE_N_DE-RELATIVED',
                    'NOTE_V2', 'ORGAN', 'ORGAN_PART_OR_LINEAGE', 'ORGAN_SYSTEM',
                    'RELATIVES', 'RELATIVES_CATEGORY', 'RELATIVES_NOTE',
                    'SAMPLE_HEALTH_STATUS', 'sample_ontology_term_high_order_MartinHirst']
# remove columns from v1_2_df_intermediate
v1_2_df_intermediate.drop(columns=deselect_columns, inplace=True)

### add Aneta Mikulasova's assignments
aneta_assignments_new = pd.read_csv('./openrefine/v1.2/IHEC_metadata_harmonization.v1.2.preliminary_AM_231218.csv')
aneta_columns = ['EpiRR', 'ORGAN_SYSTEM', 'ORGAN_PART_OR_LINEAGE', 'ORGAN', 'CELL', 'CELL_2', 'CELL_3', 'CANCER_TYPE',
                 'CANCER_SUBTYPE']

# merge aneta_assignments_new with only aneta_columns with v1_2_df_intermediate, validate one_to_one
v1_2_df_intermediate = v1_2_df_intermediate.merge(aneta_assignments_new[aneta_columns], 'left', validate='one_to_one')
# remove Epirr column from aneta_columns
aneta_columns.remove('EpiRR')
# add suffix to column names aneta_columns in v1_2_df_intermediate and make column to lower case
v1_2_df_intermediate.rename(columns={col: f'sample_{col.lower()}_order_AnetaMikulasova' for col in aneta_columns},
                            inplace=True)

### add EpiClass predictions
epiclass_predictions = pd.read_csv('./openrefine/v1.2/EpiAtlas_sex+life_stage_analyses - v1.2_metadata_integration.csv')
# merge EpiClass predictions with v1_2_df_intermediate, validate one_to_one
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    epiclass_predictions[['epirr_noV', 'EpiClass_pred_Sex', 'EpiClass_pred_Life_stage']], 'left',
    left_on='epirr_id_without_version',
    right_on='epirr_noV', validate='one_to_one')
# drop column epirr_noV
v1_2_df_intermediate.drop(columns='epirr_noV', inplace=True)
# rename EpiClass predictions columns
v1_2_df_intermediate.rename(
    columns={'EpiClass_pred_Sex': 'EpiClass_donor_sex', 'EpiClass_pred_Life_stage': 'EpiClass_donor_life_stage'},
    inplace=True)

### add manual groupings
# read manual groupings
manual_groupings = pd.read_csv('./openrefine/v1.2/IHEC Manual Grouping - Complete_Sample_metadata_harm.extended.csv')
# merge manual_groupings with v1_2_df_intermediate, validate one_to_one
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    manual_groupings[['epirr_id_without_version', 'Final_Manual_Group_Label'
                      # , 'cluster'
                      ]], 'left',
    validate='one_to_one')
# rename Final_Manual_Group_Label column
v1_2_df_intermediate.rename(columns={'Final_Manual_Group_Label': 'harmonized_sample_group_label'}, inplace=True)


### Add automated experiments from experiment metadata
# make reprocessed_metadata a wide table with epirr_id_without_version as index and experiment_type as columns
reprocessed_metadata_wide = reprocessed_metadata.pivot(index='epirr_id_without_version', columns='experiment_type',
                                                       values='uuid')
# merge columns PBAT and standard
reprocessed_metadata_wide['WGBS'] = reprocessed_metadata_wide[['PBAT', 'standard']].apply(
    lambda x: '::'.join(x.dropna().astype(str)) if not x.dropna().empty else np.nan, axis=1)
assert not reprocessed_metadata_wide['WGBS'].str.contains('::').any()
# drop columns PBAT and standard
# reprocessed_metadata_wide.drop(columns=['PBAT', 'standard'], inplace=True)
# join the columns total-RNA-Seq and mRNA-Seq with the separator ::
reprocessed_metadata_wide['RNA-Seq'] = reprocessed_metadata_wide[['total-RNA-Seq', 'mRNA-Seq']].apply(
    lambda x: '::'.join(x.dropna().astype(str)) if not x.dropna().empty else np.nan, axis=1)
# drop columns total-RNA-Seq and mRNA-Seq
# reprocessed_metadata_wide.drop(columns=['total-RNA-Seq', 'mRNA-Seq'], inplace=True)

histone_assays = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
all_assays = histone_assays + ['WGBS', 'RNA-Seq']
# count number of NaN values in each row
reprocessed_metadata_wide['n_histones'] = reprocessed_metadata_wide[histone_assays].notnull().sum(axis=1)
# first, all of the epigenomes are set to partial
reprocessed_metadata_wide['epiATLAS_status'] = 'partial'
# now, if there is ar least one histone mark AND RNA-Seq, this means there is imputed data for the other marks (also for WGBS)
# i.e., the status is set to full_imputed
reprocessed_metadata_wide.loc[(reprocessed_metadata_wide['n_histones'] >= 1) & (
    ~reprocessed_metadata_wide['RNA-Seq'].isnull()), 'epiATLAS_status'] = 'full_imputed'
reprocessed_metadata_wide.loc[
    reprocessed_metadata_wide[all_assays].isnull().sum(axis=1) == 0, 'epiATLAS_status'] = 'full'

# load imputed metadata file
imputed_metadata = pd.read_csv('./openrefine/v1.2/epiatlas_imputed_metadata.csv')
imputed_metadata.data_file_path = 'imputed'
# remove "." and suffix from epirr_id column
imputed_metadata['epirr_id_without_version'] = imputed_metadata['epirr_id'].str.split('\\.', 1, expand=True)[0]
# rename DNAMethyl to WGBS in column "mark"
imputed_metadata.loc[imputed_metadata['mark'] == 'DNAMethyl', 'mark'] = 'WGBS'
# pivot table of imputed_metadata
imputed_metadata_wide = imputed_metadata.pivot(index='epirr_id_without_version', columns='mark',
                                               values='data_file_path')
reprocessed_metadata_wide.fillna(imputed_metadata_wide, inplace=True)

# make index to column in reprocessed_metadata_wide
reprocessed_metadata_wide.reset_index(inplace=True)
# imputed_metadata_wide.reset_index(inplace=True)
#
# # merge reprocessed_metadata_wide with imputed_metadata_wide on epirr_id_without_version, validate one_to_one
# merged_metadata_wide = reprocessed_metadata_wide.merge(imputed_metadata_wide, 'left', on='epirr_id_without_version',
#                                                        validate='one_to_one')
# columns_to_combine = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3', 'WGBS']
#
# # Step 2: Combine values for each specified column with _x and _y suffixes
# for col in columns_to_combine:
#     x_col = f'{col}_x'
#     y_col = f'{col}_y'
#     # Here we concatenate non-null values into a single string with a separator "::"
#     # Adjust the lambda function for different combination logic as needed
#     merged_metadata_wide[col] = merged_metadata_wide[[x_col, y_col]].apply(
#         lambda row: '::'.join(row.dropna().astype(str)) if not row.dropna().empty else np.nan, axis=1)
#     # Drop the original _x and _y columns if they are no longer needed
#     merged_metadata_wide.drop(columns=[x_col, y_col], inplace=True)

reprocessed_metadata_wide['fully_imputed'] = reprocessed_metadata_wide[all_assays].notnull().sum(axis=1) == 8

assert reprocessed_metadata_wide.loc[
    reprocessed_metadata_wide['epiATLAS_status'].isin(['full', 'full_imputed']), 'epirr_id_without_version'].equals(
    reprocessed_metadata_wide.loc[reprocessed_metadata_wide['fully_imputed'], 'epirr_id_without_version'])

assert reprocessed_metadata_wide.loc[
    reprocessed_metadata_wide['epiATLAS_status'] == 'partial', 'epirr_id_without_version'].equals(
    reprocessed_metadata_wide.loc[~reprocessed_metadata_wide['fully_imputed'], 'epirr_id_without_version'])

# merge reprocessed_metadata_wide with v1_2_df_intermediate, validate one_to_one
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    reprocessed_metadata_wide[all_assays + ['epiATLAS_status', 'epirr_id_without_version']], 'left',
    validate='one_to_one')
# add prefix automated_experiments_ to columns H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3, WGBS, RNA-Seq
v1_2_df_intermediate.rename(columns={col: f'automated_experiments_{col}' for col in all_assays},
                            inplace=True)

### Add Fig1 ontology high order terms
# read IHEC_metadata_harmonization.v1.1.ontology_colors.csv
ontology_colors = pd.read_csv('./openrefine/v1.2/IHEC_metadata_harmonization.v1.1.ontology_colors.csv')
# add EpiRR_no_decimal
ontology_colors['EpiRR_no_decimal'] = ontology_colors['EpiRR'].str.split('\\.', 1, expand=True)[0]
# read pruned_cov_merged_ontology.csv
pruned_cov_merged_ontology = pd.read_csv('./openrefine/v1.2/pruned_cov_merged_ontology.csv')
# stack ontology_colors and pruned_cov_merged_ontology
stacked_ontology = pd.concat([ontology_colors, pruned_cov_merged_ontology])
# rename fig1_ontology_intermediate_merged column to harmonized_sample_intermediate_merged
stacked_ontology.rename(columns={'fig1_ontology_intermediate_merged': 'harmonized_sample_ontology_term_high_order_fig1',
                                 'EpiRR_no_decimal': 'epirr_id_without_version'}, inplace=True)
# merge stacked_ontology with v1_2_df_intermediate , validate one_to_one
v1_2_df_intermediate = v1_2_df_intermediate.merge(
    stacked_ontology[['epirr_id_without_version', 'harmonized_sample_ontology_term_high_order_fig1']], 'left',
    on='epirr_id_without_version', validate='one_to_one')

v1_2_extended_csv = './openrefine/v1.2/IHEC_metadata_harmonization.v1.2.extended.csv'

### Skip automated checks for now
v1_2_df_intermediate.to_csv(v1_2_extended_csv, index=False)

# for this step,
with open('./openrefine/v1.2/Routput.txt', 'w') as f:
    run(['conda', 'run', '-n', 'v1_2_r_env', 'Rscript', 'add_higher_order_v1.2.R'], check=True, cwd='./openrefine/v1.2', stdout=f)

v1_2_extended_df = pd.read_csv(v1_2_extended_csv)

v1_2_extended_df = v1_2_extended_df[["EpiRR",
                                         "project",
                                         "harmonized_biomaterial_type",
                                         "harmonized_sample_ontology_intermediate",
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
                                         "sample_ontology_term_high_order_JeffreyHyacinthe",
                                         "sample_ontology_term_high_order_JonathanSteif",
                                         "sample_organ_system_order_AnetaMikulasova",
                                         "sample_organ_part_or_lineage_order_AnetaMikulasova",
                                         "sample_organ_order_AnetaMikulasova",
                                         "sample_cell_order_AnetaMikulasova",
                                         "sample_cell_2_order_AnetaMikulasova",
                                         "sample_cell_3_order_AnetaMikulasova",
                                         "sample_cancer_type_order_AnetaMikulasova",
                                         "sample_cancer_subtype_order_AnetaMikulasova",
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
                                         "automated_experiments_H3K27ac",
                                         "automated_experiments_H3K27me3",
                                         "automated_experiments_H3K36me3",
                                         "automated_experiments_H3K4me1",
                                         "automated_experiments_H3K4me3",
                                         "automated_experiments_H3K9me3",
                                         "automated_experiments_WGBS",
                                         "automated_experiments_RNA-Seq",
                                         "epirr_id_without_version"]]

v1_2_extended_df.to_csv(v1_2_extended_csv, index=False)
v1_2_csv = './openrefine/v1.2/IHEC_metadata_harmonization.v1.2.csv'
v1_2_extended_df.loc[:,
~(v1_2_extended_df.columns.str.startswith('automated') | v1_2_extended_df.columns.str.contains('order'))].to_csv(
    v1_2_csv, index=False)

### Compare v1.1 and v1.2

old = pd.read_csv(v1_1_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
old.sort_index(1, inplace=True)
new = pd.read_csv(v1_2_extended_csv)
new.index = new.EpiRR
new.sort_index(0, inplace=True)
new.sort_index(1, inplace=True)
assert old.index.isin(new.index).all()
new = new[new.index.isin(old.index)]
assert old.columns.difference(new.columns).empty
assert (new.columns.difference(old.columns) == ['automated_experiments_H3K27ac', 'automated_experiments_H3K27me3',
                                                'automated_experiments_H3K36me3', 'automated_experiments_H3K4me1',
                                                'automated_experiments_H3K4me3', 'automated_experiments_H3K9me3',
                                                'automated_experiments_RNA-Seq', 'automated_experiments_WGBS',
                                                'epiATLAS_status', 'harmonized_sample_ontology_term_high_order_fig1',
                                                'sample_cancer_subtype_order_AnetaMikulasova',
                                                'sample_cancer_type_order_AnetaMikulasova',
                                                'sample_cell_2_order_AnetaMikulasova',
                                                'sample_cell_3_order_AnetaMikulasova',
                                                'sample_cell_order_AnetaMikulasova',
                                                'sample_organ_order_AnetaMikulasova',
                                                'sample_organ_part_or_lineage_order_AnetaMikulasova',
                                                'sample_organ_system_order_AnetaMikulasova']).all()
shared_cols = old.columns.intersection(new.columns)
assert shared_cols.size == 42

diff_tbl = old[shared_cols].compare(new[shared_cols])
diff_tbl.rename(columns={'self': 'v1.1', 'other': 'v1.2'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v1.2/diff_v1.1_v1.2.json', indent=True)
