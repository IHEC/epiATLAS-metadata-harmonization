import difflib
import json
import os.path
import typing
from subprocess import run

import pandas as pd
from pronto import Ontology

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.5/IHEC_metadata_harmonization.v0.5.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]

# create project with old version
run([openrefine_client, '--create', initial_csv], check=True)

# apply rules
run([openrefine_client, '--apply', 'openrefine/v0.6/metadata_jamboree_15022022_solving_github_issues_in_v0.5.json',
     old_project_name], check=True)

# export to intermediate version, because we need to do some cleaning in python now
intermediate_csv = "./openrefine/v0.6/IHEC_metadata_harmonization.v0.6.intermediate.csv"

run([openrefine_client, '--export', f'--output={intermediate_csv}',
     old_project_name], check=True)

# read csv to pandas
intermediate_tbl = pd.read_csv(intermediate_csv)
separator = '::'


# use pronto to check the corresponding ontologies this will generate warnings!
def get_term2id(ontology: Ontology) -> typing.Tuple[typing.Dict[str, str], typing.List[str]]:
    terms = []
    term2id = {}
    for term in ontology.terms():
        if term.name:
            terms.append(term.name)
            term2id[term.name] = term.id
    return term2id, terms


# make efo separately because pronto had trouble parsing it:
efo_terms: list = []
efo = {}
efo_term2id = {}
with open('ontologies/efo.obo') as efo_file:
    efo_file = efo_file.read()
line_iterator = iter(efo_file.splitlines())
line = next(line_iterator, None)
while line is not None:
    if line == '[Term]':
        curr_id = next(line_iterator)[4:]
        if curr_id.startswith('EFO'):
            curr_name = next(line_iterator)[6:]
            efo_terms.append(curr_name)
            efo[curr_id] = curr_name
            efo_term2id[curr_name] = curr_id
    line = next(line_iterator, None)

cl = Ontology('ontologies/cl-basic.obo')
cl_term2id, cl_terms = get_term2id(cl)
cl_dict = {'ontology': cl, 'terms': cl_terms, 'terms_lower': [t.lower() for t in cl_terms], 'term2id': cl_term2id,
           'column': 'cell_type'}
uberon = Ontology('ontologies/uberon-basic.obo')
uberon_term2id, uberon_terms = get_term2id(uberon)
biomaterial2onto = {'primary cell': cl_dict,
                    'primary cell culture': cl_dict,
                    'primary tissue': {'ontology': uberon, 'terms': uberon_terms,
                                       'terms_lower': [t.lower() for t in uberon_terms], 'term2id': uberon_term2id,
                                       'column': 'tissue_type'},
                    'cell line': {'ontology': efo, 'terms': efo_terms, 'terms_lower': [t.lower() for t in efo_terms],
                                  'term2id': efo_term2id, 'column': 'line'}
                    }

unmatched: typing.Dict[str, dict] = {}
fill_sample_ontology_term: dict = {}

# check sample ontologies
for index, row in intermediate_tbl.iterrows():
    if row.EpiRR.startswith('IHECRE00000922'):
        # this entry was processed seperately because there are multiple tissues it came from
        continue
    is_efo = row.biomaterial_type == 'cell line'
    onto_dict: dict = biomaterial2onto[row.biomaterial_type]
    name_column = onto_dict['column']
    onto_id: str = row.sample_ontology_term
    current_name = row[name_column]
    if not pd.isna(onto_id):  # there exists an ontology curie
        curr_term_name = onto_dict['ontology'][onto_id]
        if not is_efo:
            curr_term_name = curr_term_name.name
        if pd.isna(current_name):
            intermediate_tbl.loc[index, name_column] = curr_term_name
        elif current_name.casefold() != curr_term_name.casefold():
            unmatched[row.EpiRR] = {'current': current_name, 'term_name': curr_term_name, 'ontology_curie': onto_id,
                                    'similarity': round(difflib.SequenceMatcher(None, current_name.casefold(),
                                                                                curr_term_name.casefold()).ratio(), 2)}
        else:
            # this means the term are equal.
            # We still replace the name by the one from the ontology in case there are capitalization problems
            # e.g. CD4-positive or cd4-positive
            intermediate_tbl.loc[index, name_column] = curr_term_name

    else:
        if pd.isna(current_name):
            fill_sample_ontology_term[row.EpiRR] = []
        else:
            try:
                intermediate_tbl.loc[index, 'sample_ontology_term'] = onto_dict['term2id'][current_name]
            except KeyError:
                try:
                    # lowercase instead of correct capitalization
                    index_lower = onto_dict['terms_lower'].index(current_name)
                    correct_name = onto_dict['terms'][index_lower]
                    intermediate_tbl.loc[index, name_column] = correct_name
                    intermediate_tbl.loc[index, 'sample_ontology_term'] = onto_dict['term2id'][correct_name]
                except ValueError:
                    # no ontology curie -> check if the current term name has similar entries in the ontology
                    fill_sample_ontology_term[row.EpiRR] = {'current': current_name,
                                                            'similiar_terms': [
                                                                {'ontology_curie': x,
                                                                 'term_name': onto_dict['term2id'][x]}
                                                                for x in
                                                                difflib.get_close_matches(current_name,
                                                                                          onto_dict['terms'])]}


# we can fix overall 745 falsy entries using this
print(f'we can fix overall {len(pd.read_csv(intermediate_csv).compare(intermediate_tbl))} falsy entries using this')

# there are still 127 na's in sample_ontology_term, 17 of these also do not have a name (e.g. cell_type)
print(
    f'overall {len(fill_sample_ontology_term)} nas in sample_ontology_term, {len([(k, v) for k, v in fill_sample_ontology_term.items() if v == []])} of these also do not have a name (e.g. cell_type)')

# 541 entries have conflicting information in the columns
print(f'{len(unmatched)} entries have conflicting information in the columns')

with open('openrefine/v0.6/fill_sample_ontology_term.json', 'w') as f:
    json.dump(fill_sample_ontology_term, f, indent=True)

with open('openrefine/v0.6/conflicts_sample_ontology_term.json', 'w') as f:
    json.dump(unmatched, f, indent=True)

# fix case for donor_id:
original_table = pd.read_csv('raw/EpiAtlas_EpiRR_metadata_all.csv')
with open('config/merges.json') as rule_file:
    rules = json.load(rule_file)
cols_to_use = next(
    (rule['strategy']['merge_delimited'] for rule in rules if rule['harmonized'] == 'donor_id'))
donor_id_case_sensitive = original_table[cols_to_use].apply(
    lambda x: separator.join(sorted(set(x.dropna().astype(str)))), axis=1)

intermediate_tbl.donor_id = donor_id_case_sensitive
intermediate_tbl.to_csv(intermediate_csv)

run([openrefine_client, '--create', intermediate_csv], check=True)
intermediate_project_name = os.path.splitext(os.path.basename(intermediate_csv))[0]

# apply rules
run([openrefine_client, '--apply', 'openrefine/v0.6/fix_donor_id.json',
     intermediate_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.6/remove_columns.json',
     intermediate_project_name], check=True)

v06_csv = './openrefine/v0.6/IHEC_metadata_harmonization.v0.6.csv'
run([openrefine_client, '--export', f'--output={v06_csv}',
     intermediate_project_name], check=True)
