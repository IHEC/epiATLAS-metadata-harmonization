import json

import pandas as pd

conflicts = pd.read_json('../v0.6/conflicts_sample_ontology_term.json', orient='records').T
actual_columns = conflicts.columns.to_list()
agg_conflicts = conflicts.reset_index().groupby(actual_columns).agg({'index': ' ,'.join}).reset_index().sort_values(
    'index')
agg_conflicts.to_json('conflicts_aggregated.json', orient='records', indent=True)
agg_conflicts.to_csv('conflicts_aggregated.csv', index=False)

to_fill = json.load(open('../v0.6/fill_sample_ontology_term.json'))
tmp_list = []
for k, v in to_fill.items():
    tmp_dict = {'index': k}
    if 'current' in v:
        tmp_dict['current'] = v['current']
    if 'similiar_terms' in v:
        for i, d in enumerate(v['similiar_terms']):
            tmp_dict.update({k_i + str(i + 1): v_i for k_i, v_i in d.items()})
    tmp_list.append(tmp_dict)
fill_df = pd.DataFrame.from_dict(tmp_list)
agg_fill = fill_df.groupby(fill_df.columns.difference({'index'}).to_list(), dropna=False).agg(
    {'index': ' ,'.join}).reset_index().sort_values('index')
agg_fill.to_json('fill_aggregated.json', orient='records', indent=True)
agg_fill.to_csv('fill_aggregated.csv', index=False)
