import os.path
from subprocess import run

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.6/IHEC_metadata_harmonization.v0.6.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]

# create project with old version
run([openrefine_client, '--create', initial_csv], check=True)

# apply rules
run([openrefine_client, '--apply', 'openrefine/v0.7/solving_DEEP.json', old_project_name], check=True)

v07_csv = './openrefine/v0.7/IHEC_metadata_harmonization.v0.7.csv'

run([openrefine_client, '--export', f'--output={v07_csv}', old_project_name], check=True)
