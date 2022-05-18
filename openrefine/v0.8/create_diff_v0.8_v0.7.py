import os
import pandas as pd

# make sure the working directory when running this file is the project root of the git repository
os.chdir('../../')

old = pd.read_csv('openrefine/v0.7/IHEC_metadata_harmonization.v0.7.csv')
old.index = old.EpiRR
new = pd.read_csv('openrefine/v0.8/IHEC_metadata_harmonization.v0.8.csv')
new.index = new.EpiRR

removed_entries = old.EpiRR[~old.EpiRR.isin(new.EpiRR)]
removed_entries.to_csv('openrefine/v0.8/removed_entries.csv', index=False)

old = old[~old.EpiRR.isin(removed_entries)]
old.drop(columns='taxon_id', inplace=True)
old['disease_ontology_term'] = old['disease_ontology_term'].str.upper()
old['donor_health_status2'] = old['donor_health_status']
old['disease_ontology_term2'] = old['disease_ontology_term']

old.rename(columns={'donor_health_status': 'v0.7:donor_health_status:v0.8:donor_health_status',
                    'disease_ontology_term': 'v0.7:disease_ontology_term:v0.8:donor_health_status_ontology_curie',
                    'donor_health_status2': 'v0.7:donor_health_status:v0.8:disease',
                    'disease_ontology_term2': 'v0.7:disease_ontology_term:v0.8:disease_ontology_curie',
                    'sample_ontology_term': 'v0.7:sample_ontology_term:v0.8:sample_ontology_curie'}, inplace=True)
old.sort_index(1, inplace=True)

new.rename(columns={'donor_health_status': 'v0.7:donor_health_status:v0.8:donor_health_status',
                    'disease': 'v0.7:donor_health_status:v0.8:disease',
                    'donor_health_status_ontology_curie': 'v0.7:disease_ontology_term:v0.8:donor_health_status_ontology_curie',
                    'disease_ontology_curie': 'v0.7:disease_ontology_term:v0.8:disease_ontology_curie',
                    'sample_ontology_curie': 'v0.7:sample_ontology_term:v0.8:sample_ontology_curie'}, inplace=True)
new.sort_index(1, inplace=True)

diff_tbl = old.compare(new)
diff_tbl.rename(columns={'self': 'v0.7', 'other': 'v0.8'}, inplace=True)
diff_tbl.apply(lambda x: [x.dropna()], axis=1).to_json('openrefine/v0.8/diff_v0.7_v0.8.json', indent=True)
