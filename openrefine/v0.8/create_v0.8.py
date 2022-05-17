import os.path
from subprocess import run

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

import pandas as pd
import numpy as np
import sys
from pronto import Ontology


def create_agreement_cols(df_to_work_filled_merge):
    """Receives a df and save a new one as
    a csv file including the Diff columns
    (agree/disagree with ontology description)"""

    df_to_work_filled_merge['Diff_health'] = np.where(
        df_to_work_filled_merge['donor_health_status'] == df_to_work_filled_merge['Descrip_onto_health_status'],
        'agree', 'disagree')
    df_to_work_filled_merge['Diff_health_merge'] = np.where(
        df_to_work_filled_merge['donor_health_status_merge'] == df_to_work_filled_merge['Descrip_onto_health_status'],
        'agree', 'disagree')
    df_to_work_filled_merge['Diff_disease'] = np.where(
        df_to_work_filled_merge['disease'] == df_to_work_filled_merge['Descrip_onto_disease'], 'agree', 'disagree')

    # reorder columns
    cols = ['EpiRR', 'EpiRR_status', 'age', 'biomaterial_type', 'cell_type',
            'disease_ontology_term', 'donor_age_unit',
            'donor_id', 'donor_life_stage', 'health_state', 'line', 'markers',
            'project', 'sample_ontology_term', 'sex', 'taxon_id', 'tissue_type',
            'donor_health_status_merge', 'donor_health_status', 'Descrip_onto_health_status',
            'Diff_health', 'Diff_health_merge', 'donor_health_status_ontology_uri',
            'donor_health_status_ontology_curie', 'disease', 'Descrip_onto_disease', 'Diff_disease',
            'disease_ontology_uri', 'disease_ontology_curie', 'Descrip_all_terms']

    df_to_work_filled_merge = df_to_work_filled_merge[cols]
    df_to_work_filled_merge.to_csv('openrefine/v0.8/IHEC_diff_onto_disease_health_ncitowl_2.csv', index=False)  # saving file to explore


def fill_health_disease(epirr_amed_113, df_to_work):
    df_to_work_filled_merge = fill_amed_crest_merge(epirr_amed_113, df_to_work)

    # Fill na disease cells with descrp_disease
    df_to_work_filled_merge['disease'] = np.where(
        (pd.isna(df_to_work_filled_merge['disease']) & pd.notna(df_to_work_filled_merge['Descrip_onto_disease'])),
        df_to_work_filled_merge['Descrip_onto_disease'], df_to_work_filled_merge['disease'])
    # Fill na heath status cells with descrp_health
    df_to_work_filled_merge['donor_health_status'] = np.where((pd.isna(
        df_to_work_filled_merge['donor_health_status']) & pd.notna(
        df_to_work_filled_merge['Descrip_onto_health_status'])), df_to_work_filled_merge['Descrip_onto_health_status'],
                                                              df_to_work_filled_merge['donor_health_status'])
    # Fill na health status cells with descrp_all when disease is na
    df_to_work_filled_merge['donor_health_status'] = np.where(
        (pd.isna(df_to_work_filled_merge['donor_health_status']) & pd.notna(df_to_work_filled_merge['disease'])),
        df_to_work_filled_merge['Descrip_all_terms'], df_to_work_filled_merge['donor_health_status'])

    return df_to_work_filled_merge


def fill_amed_crest_merge(epirr_amed_113, df_to_work):
    # whole dict
    dict_merge = dict(zip(df_to_work['EpiRR'], df_to_work['donor_health_status_merge']))
    # list to filter dict (113 - AMED-CREST)
    list_to_dict = [line.strip() for line in epirr_amed_113]
    # filtering dict merge to map health column with new info
    dict_to_map = {my_key: dict_merge[my_key] for my_key in list_to_dict}
    # map dict
    df_to_work['donor_health_status'] = df_to_work['EpiRR'].map(dict_to_map)

    return df_to_work


def map_term_ncit(df_to_work, ncit_obo, ncit_dat):
    """Receives a df with the desired columns
    and two dictionaries (from .obo and .dat
    files). Return a df including the description
    NCI terms columns."""

    dict_terms_nci = create_ncit_obo_dict(ncit_obo)
    dict_dat = create_dict_dat(ncit_dat)

    # list ontologies in disease_ontology_term merged
    dis_ont = df_to_work['disease_ontology_term'].str.split(':').str[-1].tolist()  # problem merged ::
    # list ontologies in disease_ontology_uri
    dis_ont_uri = df_to_work['disease_ontology_uri'].str.split('code=').str[-1].str.split('&').str[0].to_list()
    df_to_work['disease_ontology_uri'] = dis_ont_uri  # reassigning values
    # list ontologies in donor_health_status_uri
    dhealth_ont_uri = df_to_work['donor_health_status_ontology_uri'].str.split('code=').str[-1].str.split('&').str[
        0].tolist()
    df_to_work['donor_health_status_ontology_uri'] = dhealth_ont_uri

    # list to be columns
    health_ont_desc = []
    disease_ont_desc = []
    general_ont_desc = []  # to create as well

    health_ont_desc = create_description_col(dis_ont, dhealth_ont_uri, dict_terms_nci, dict_dat, health_ont_desc)
    disease_ont_desc = create_description_col(dis_ont, dis_ont_uri, dict_terms_nci, dict_dat, disease_ont_desc)

    # list with all descriptions available for all terms in disease_ontology_term
    for ele in dis_ont:
        if pd.notna(ele):

            if ele not in dict_terms_nci.keys():
                if ele == 'C0277545':
                    general_ont_desc.append('Disease type AND/OR category unknown')
                else:
                    term = dict_dat.get(ele)  # getting value (ncit term) for CUI
                    general_ont_desc.append(dict_terms_nci.get(term, "No description available"))
            else:
                general_ont_desc.append(dict_terms_nci.get(ele, "No description available"))
        else:
            general_ont_desc.append(np.nan)

    # creating columns
    df_to_work['Descrip_onto_health_status'] = health_ont_desc
    df_to_work['Descrip_onto_disease'] = disease_ont_desc
    df_to_work['Descrip_all_terms'] = general_ont_desc

    return df_to_work


def create_description_col(list_ont, list_desired_col, dict_terms_nci, dict_dat, list_to_col):
    """Inputs: list ont, list_desired_col, dicts
       Output: list_to_col. The two input lists have
       the terms related with disease_ont_term and health
       or disease uri colums. The dicts were generated from
       .obo and .dat files. The output list constains the
       description of each term"""

    for d_ont, d_uri in zip(list_ont, list_desired_col):
        if pd.notna(d_ont) and pd.notna(d_uri):

            if d_uri not in dict_terms_nci.keys():
                if d_uri == 'C0277545':
                    list_to_col.append('Disease type AND/OR category unknown')
                else:
                    term = dict_dat.get(d_uri)  # getting value (ncit term) for CUI in .dat dict
                    list_to_col.append(dict_terms_nci.get(term, "No description available"))
            else:
                list_to_col.append(
                    dict_terms_nci.get(d_uri, "No description available"))  # not working to replace None
        else:
            list_to_col.append(np.nan)

    return list_to_col


def merge_col_ori_version(df_raw, df_v7):  # ok
    """Receives two dfs and returns a
    merged df (e.g IHEC metadata raw
    and v7)"""

    df_ori_col = get_desired_cols(df_raw)  # v01 - correct
    df_v7.rename({'donor_health_status': 'donor_health_status_merge'}, axis=1, inplace=True)  # v07
    df_to_work = df_v7.merge(df_ori_col, on='EpiRR', how='left')  # merging ori columns

    return df_to_work


def get_desired_cols(df_raw):  # original version
    """Receives a df and returns
    a df with specific
    columns"""

    return df_raw[['EpiRR', 'donor_health_status', 'disease', 'donor_health_status_ontology_curie',
                   'donor_health_status_ontology_uri', 'disease_ontology_curie',
                   'disease_ontology_uri']]


def create_dict_dat(ncit_dat):
    """Receives a .dat files and
    returns a dict (Keys=CUI term;
    values=NCIT term to map on
    obo dict)"""

    dict_dat = {}
    for line in ncit_dat:
        dict_dat[line.strip().split('|')[1]] = line.strip().split('|')[0]  # keys CUI; val NCIT term

    return dict_dat


def create_ncit_obo_dict(ncit_obo):
    """Receives a .obo file and
    returns a dict (Keys=terms;
    values=description_term)"""

    dict_terms_nci = {}

    for term in ncit_obo.terms():  # from obo file
        if term.name:
            dict_terms_nci[term.id.split('#')[-1]] = term.name

    return dict_terms_nci


def generate_obo_file(owl_file):
    """Receives a OWL file and
    returns a .obo file"""

    with open("ontologies/ncit_2204.obo", "wb") as f:
        owl_file.dump(f, format="obo")



print('Starting script....')

if not os.path.exists('ontologies/ncit_2204.obo'):
    # make sure to gunzip the owl and obo file, otherwise a memory intensive computation is run
    ncit_owl = Ontology("ontologies/Thesaurus.owl") # from https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus_22.04d.OWL.zip
    generate_obo_file(ncit_owl)

ncit_obo = Ontology('ontologies/ncit_2204.obo')
ncit_dat = open('ontologies/nci_code_cui_map_201706.dat.txt', 'r')  # file .dat downloaded from ftp server (NCI) https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/nci_code_cui_map/nci_code_cui_map_201706.dat
df_raw = pd.read_csv('raw/EpiAtlas_EpiRR_metadata_all.csv')  # version 1 (EpiRR updated - no merged columns)
df_v7 = pd.read_csv('openrefine/v0.7/IHEC_metadata_harmonization.v0.7.csv')  # version 7 (last version)
epirr_amed_113 = open('openrefine/v0.8/list_EpiRR_tofill_healthstatus.txt', 'r')
df_to_work = merge_col_ori_version(df_raw, df_v7)
df_to_work_1 = map_term_ncit(df_to_work, ncit_obo, ncit_dat)
df_to_work_filled_merge = fill_health_disease(epirr_amed_113, df_to_work_1)

print('Saving df...')
create_agreement_cols(df_to_work_filled_merge)
print('Finished!')


# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.8/IHEC_diff_onto_disease_health_ncitowl_2.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]

# create project with old version
run([openrefine_client, '--create', initial_csv], check=True)

# apply rules
run([openrefine_client, '--apply', 'openrefine/v0.8/01-IHEC_diff_onto_disease_health_ncitowl_2.json', old_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.8/02-resolving_github_issues.json', old_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.8/03-filter_and_rename_rows_and_cols.json', old_project_name], check=True)

v08_csv = './openrefine/v0.8/IHEC_metadata_harmonization.v0.8.csv'

run([openrefine_client, '--export', f'--output={v08_csv}', old_project_name], check=True)
