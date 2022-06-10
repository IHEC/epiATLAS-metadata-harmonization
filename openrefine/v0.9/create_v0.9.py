import json
import os.path
import warnings
from subprocess import run
import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# first get the ncit thesaurus
from typing import Dict, List, Tuple
import xml.etree.ElementTree as et

tree = et.parse('ontologies/Thesaurus.owl')
root = tree.getroot()
classes: List[et.Element] = root.findall('{http://www.w3.org/2002/07/owl#}Class')
names: List[str] = []
name2ncim: Dict[str, List[str]] = dict()
name2ncit: Dict[str, List[str]] = dict()
ncit2name: Dict[str, str] = dict()
ncit2ncim: Dict[str, str] = dict()
ncim2name: Dict[str, str] = dict()

for currClass in classes:
    ncitElement = currClass.find('{http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#}NHC0')
    ncit = ncitElement.text
    nameElement = currClass.find('{http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#}P108')
    name = None
    if nameElement is not None:
        name = nameElement.text
        names.append(name)
        ncit2name[ncit] = name
        if name in name2ncit:
            name2ncit[name].append(ncit)
        else:
            name2ncit[name] = [ncit]
    ncimElement = currClass.find('{http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#}P207')
    if ncimElement is None:
        ncimElement = currClass.find('{http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#}P208')
    if ncimElement is not None:
        ncim = ncimElement.text
        ncit2ncim[ncit] = ncim
        if nameElement is not None:
            ncim2name[ncim] = name
            if name in name2ncim:
                name2ncim[name].append(ncim)
            else:
                name2ncim[name] = [ncim]

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.8/IHEC_metadata_harmonization.v0.8.csv'  # csv to build project from

# read csv to pandas
working_tbl = pd.read_csv(initial_csv)
separator = '::'

disagree: Dict[Tuple[str, str], Dict[str, str]] = {}
mapping_problems: Dict[Tuple[str, str], Dict[str, str]] = {}
free_text_cols: List[str] = ['donor_health_status', 'disease']
onto_cols: List[str] = ['donor_health_status_ontology_curie', 'disease_ontology_curie']
both_missing: List[str] = []

# check disease and donor_health_status
for index, row in working_tbl.iterrows():
    for free_text_col, onto_col in zip(free_text_cols, onto_cols):
        free_text = row[free_text_col]
        onto = row[onto_col]
        if not pd.isna(onto):  # id exists
            cuis = [curie.split(':')[1] for curie in
                    onto.split(separator)]  # split by separator in case we have multiple
            cuis = [ncit2ncim[cui] if cui in ncit2ncim else cui for cui in cuis]
            ontos_ncim = sorted([ncim2name[cui] for cui in cuis if cui in ncim2name])  # map assuming its a ncim cui
            ontos_ncit = sorted([ncit2name[ncit] for ncit in cuis if ncit in ncit2ncim])
            if ontos_ncit:
                warnings.warn(f'{row.EpiRR}: NCIT ID(s) could not be mapped to NCIM')
            if len(ontos_ncim) != len(cuis):
                # print(f'{row.EpiRR} potential problem when mapping ncim to name')
                mapping_problems[(row.EpiRR, free_text_col)] = {'free_text': free_text, 'ontology_curie': onto,
                                                                'description': 'ncim not found'}
                continue
            if not pd.isna(free_text):  # free text exists
                free_texts = sorted(free_text.split(separator))
                if all([a.casefold() == b.casefold() for a, b in zip(ontos_ncim, free_texts)]):
                    working_tbl.loc[index, free_text_col] = separator.join(ontos_ncim)
                else:
                    disagree[(row.EpiRR, free_text_col)] = {'free_text': free_text, 'ontology_curie': onto,
                                                            'ontology_term': separator.join(ontos_ncim)}
                    # print(f'{row.EpiRR} onto and free text disagree')
            else:  # since free text is empty we put in the name from the ontos
                working_tbl.loc[index, free_text_col] = separator.join(ontos_ncim)
        elif not pd.isna(free_text):  # we have no onto_id but we have a free_text
            free_texts = sorted(free_text.split(separator))
            ncims = sorted([sorted(name2ncim[free_text]) for free_text in free_texts if free_text in name2ncim])
            name_not_available = sorted([free_text for free_text in free_texts if free_text not in name2ncim])
            if name_not_available:
                mapping_problems[(row.EpiRR, free_text_col)] = {'free_text': free_text, 'ontology_curie': onto,
                                                                'description': 'name not found'}
                continue
            working_tbl.loc[index, onto_col] = separator.join(
                [separator.join(['NCIM:' + ncim for ncim in ncim_name]) for ncim_name in ncims])
        else:
            both_missing.append(row.EpiRR)

disagree_df = pd.DataFrame.from_dict(disagree, orient='index')
disagree_df.index.names = ['EpiRR', 'field']
agg_disagree = disagree_df.reset_index().groupby(['field', 'free_text', 'ontology_curie', 'ontology_term'])[
    'EpiRR'].apply('|'.join).reset_index()
agg_disagree.to_csv('./openrefine/v0.9/disagree.csv', index=False)

mapping_problems_df = pd.DataFrame.from_dict(mapping_problems, orient='index')
mapping_problems_df.index.names = ['EpiRR', 'field']
agg_problems = mapping_problems_df.reset_index().groupby(['field', 'free_text', 'ontology_curie', 'description'])[
    'EpiRR'].apply('|'.join).reset_index()
agg_problems.to_csv('./openrefine/v0.9/mapping_problems.csv', index=False)

initial_tbl = pd.read_csv(initial_csv)
diff_tbl = ~((initial_tbl == working_tbl) | ((initial_tbl != initial_tbl) & (working_tbl != working_tbl)))
diff_tbl.sum().to_csv('./openrefine/v0.9/changed_entries.csv', header=False)

intermediate_csv = './openrefine/v0.9/IHEC_metadata_harmonization.v0.9.intermediate.csv'
working_tbl.to_csv(intermediate_csv, index=False)

# create project with intermediate version
run([openrefine_client, '--create', intermediate_csv], check=True)
intermediate_project_name = os.path.splitext(os.path.basename(intermediate_csv))[0]

# here we manually solve some mapping issues and conflicts and the resulting json is then used in this script

run([openrefine_client, '--apply', 'openrefine/v0.9/mapping_and_conflict_fixes.json', intermediate_project_name],
    check=True)
run([openrefine_client, '--apply', 'openrefine/v0.9/DEEP_disease_donor_health_status.json', intermediate_project_name],
    check=True)

v09_csv = './openrefine/v0.9/IHEC_metadata_harmonization.v0.9.csv'

run([openrefine_client, '--export', f'--output={v09_csv}', intermediate_project_name], check=True)
