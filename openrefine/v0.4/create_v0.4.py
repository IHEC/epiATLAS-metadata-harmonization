import json
import os.path
from copy import deepcopy
from os import listdir
from os.path import join, isfile
from subprocess import run

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')
mypath = './openrefine/v0.4'

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.3/IHEC_metadata_harmonization.v0.3.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]

run([openrefine_client, '--create', initial_csv], check=True)

run([openrefine_client, '--apply', 'openrefine/v0.4/metadata_jamboree_20012022_solving_github_issues_in_v0.3.json', old_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.4/cleanToBeUpdated.json', old_project_name], check=True)
run([openrefine_client, '--apply', 'openrefine/v0.4/removeIDs.json', old_project_name], check=True)

new_csv = "./openrefine/v0.4/IHEC_metadata_harmonization.v0.4.csv"
run([openrefine_client, '--export', f'--output={new_csv}',
     old_project_name], check=True)

run([openrefine_client, '--create', new_csv], check=True)
