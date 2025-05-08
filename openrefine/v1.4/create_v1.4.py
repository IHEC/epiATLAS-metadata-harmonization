import os

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

old_version = 'v1.3'
new_version = 'v1.4'
old_extended_csv = f'./openrefine/{old_version}/IHEC_sample_metadata_harmonization.{old_version}_extended.csv'  # csv to build project from
new_extended_csv = old_extended_csv.replace(old_version, new_version)  # csv to build project from
new_csv = new_extended_csv.replace("_extended", "")

new_extended_df = pd.read_csv(old_extended_csv)

# some epirrs have to be reverted to "unknown"
revert_lifestage_file = f"./openrefine/{new_version}/list_EpiRR_to_revert_life-stage_to_unknown.txt"
# read all lines from the file and remove the new line character and add the to a set called revert_epirrs
with open(revert_lifestage_file, 'r') as f:
    revert_epirrs = {line.strip() for line in f.readlines()}

# revert the life stage of the epirrs in the list
new_extended_df.loc[
    new_extended_df.epirr_id_without_version.isin(revert_epirrs), 'harmonized_donor_life_stage'] = 'unknown'

# write the DataFrame to a csv file
new_extended_df.to_csv(new_extended_csv, index=False)

updated_column_mapping = {
    "EpiRR_ordering": "EpiRR Ordering",
    "EpiRR": "EpiRR",
    "harmonized_sample_disease_high": "Biospecimen Disease",
    "harmonized_sample_ontology_term_high_order_fig1": "Broad Biospecimen Label",
    "harmonized_sample_ontology_term_high_order_fig1_color": "Broad Colour",
    "harmonized_sample_ontology_intermediate": "Intermediate Biospecimen Label",
    "harmonized_sample_ontology_intermediate_color": "Intermediate Colour",
    "harmonized_sample_label": "Curated Biospecimen Label"
}

# write the DataFrame to a csv file using the names in updated_column_mapping and only those columns
new_extended_df.loc[:, list(updated_column_mapping.keys())].rename(columns=updated_column_mapping).to_csv(new_csv,
                                                                                                          index=False)

### Compare v1.2 and v1.3

old = pd.read_csv(old_extended_csv)
old.index = old.EpiRR
old.sort_index(0, inplace=True)
new = pd.read_csv(new_extended_csv)
new.index = new.EpiRR
new.sort_index(0, inplace=True)
assert old.index.equals(new.index)

old.sort_index(1, inplace=True)
new.sort_index(1, inplace=True)
assert old.columns.equals(new.columns)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': old_version, 'other': new_version}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json(
    f'openrefine/{new_version}/diff_{old_version}_{new_version}.json', indent=True)
# print column in markdown format

na_props = pd.DataFrame.from_dict({col: {'Column': col, 'Examples': '', 'Explanation': '',
                                         '# Not Null (%)': f'{new_extended_df[col].notnull().sum()} ({new_extended_df[col].notnull().mean() * 100:.1f}%)'}
                                   for col in new_extended_df.columns}, orient='index')
na_props.to_markdown(f'openrefine/{new_version}/README.md', index=False, tablefmt='github', mode='a')
