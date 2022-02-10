import json
import os.path
from copy import deepcopy
from os import listdir
from os.path import isfile, join, splitext
from subprocess import run

import pandas as pd

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')
mypath = './openrefine/v0.1'  # path to the rule files
example_filename = 'example.json'  # file containing the example rules in the OpenRefine project
delimiter = '::'

rule_files = [f for f in listdir(mypath) if
              isfile(join(mypath, f)) and f != example_filename and (
                      f.endswith('.txt') or f.endswith('.json') or f.endswith('JSON'))]

all_rules = {f: json.load(open(join(mypath, f), 'r')) for f in rule_files}
example_rules = json.load(open(join(mypath, example_filename), 'r'))

# for this to work properly the rule files have to be named the following way:
# all edited columns have to at the beginning of the filename, separated by the delimiter (here '::')
# the columns are followed by the project or author if applicable
file2columns = {file: set(splitext(file)[0].split(sep='-', maxsplit=1)[0].split(delimiter)) for file in rule_files}


def sort_merged(my_str: str, sep: str) -> str:
    return sep.join(sorted(my_str.split(sep)))


# check column dependencies for each file, sort merged fields in rules, and count #edits
dependencies = {}
edit_count = 0
for file, rules in all_rules.items():
    edit_count += len(rules)
    if rules[:len(example_rules)] == example_rules:
        edit_count -= len(example_rules)
    dependencies[file] = set()
    rules_sorted_merge = []
    for rule in rules:
        if 'columnName' in rule:
            if rule['columnName'] not in file2columns[file] and rule not in example_rules:
                # check whether this file is editing another column than specified in the filename
                print(f'ATTENTION: file {file} changes column {rule["columnName"]}')

        # because the merged fields were not sorted in the version, for which rules were created first,
        # we need to sort all merged fields
        rule_sorted = deepcopy(rule)
        conf = rule_sorted['engineConfig']
        if 'facets' in conf:
            for facet in conf['facets']:
                dependencies[file] |= {facet['columnName']}
                if 'selection' in facet:
                    for sel in facet['selection']:
                        if delimiter in sel['v']['v'] or delimiter in sel['v']['l']:
                            sel['v']['v'] = sort_merged(sel['v']['v'], delimiter)
                            sel['v']['l'] = sort_merged(sel['v']['l'], delimiter)
        if 'edits' in rule_sorted:
            for edit in rule_sorted['edits']:
                if any([(delimiter in s) for s in edit['from']]) or delimiter in edit['to']:
                    edit['from'] = [sort_merged(s, delimiter) for s in edit['from']]
                    edit['to'] = sort_merged(edit['to'], delimiter)

        rules_sorted_merge.append(rule_sorted)

    with open(join(mypath, f'{file}.merged_sorted'), 'w') as f:
        json.dump(rules_sorted_merge, f, indent=2)
    dependencies[file] -= file2columns[file]
print(f'{edit_count=}')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
raw_file_sorted = 'merged/EpiAtlas_EpiRR_metadata_all.merged_minimal_sorted.json'
raw_file_csv_sorted = raw_file_sorted + '.csv'
df = pd.read_json(raw_file_sorted)
df.to_csv(raw_file_csv_sorted, index=False)
openrefine_name_sorted = os.path.splitext(os.path.basename(raw_file_csv_sorted))[0]
run([openrefine_client, '--create', raw_file_csv_sorted],
    check=True)

remaining_rules = set(rule_files)
order = []
while remaining_rules:
    last_remaining_rules = deepcopy(remaining_rules)
    for file in sorted(remaining_rules):
        not_a_dependency = True
        for other_file in sorted(remaining_rules):
            if file != other_file:  # otherwise continue
                if bool(file2columns[file] & dependencies[other_file]):
                    not_a_dependency = False
                    break
        if not_a_dependency:
            order.append(file)
            run([openrefine_client, '--apply', (join(mypath, file) + '.merged_sorted'),
                 openrefine_name_sorted],
                check=True)
    remaining_rules -= set(order)
    assert last_remaining_rules != remaining_rules
print(f'order of rules applied: {order}')

run([openrefine_client, '--export', f'--output={join(mypath, "IHEC_metadata_harmonization.v0.1.csv")}',
     openrefine_name_sorted], check=True)

run([openrefine_client, '--export', f'--output={join(mypath, "IHEC_metadata_harmonization.v0.1.xls")}',
     openrefine_name_sorted], check=True)
