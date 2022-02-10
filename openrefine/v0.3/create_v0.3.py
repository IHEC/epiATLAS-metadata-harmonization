import os.path
from subprocess import run

# make sure the working directory when running this file is the project root of the git project
os.chdir('../../')
rule_file = './openrefine/v0.3/metadata_jamboree_02122021_solving_github_issues_in_v0.2.json'  # path to the rule file

# create openrefine project and apply rules - OPENREFINE SERVER HAS TO BE RUNNING
# creating openrefine projects via the openrefine-client needs a csv as input in order to work properly
openrefine_client = './openrefine/openrefine-client_0-3-10_linux'  # path to the openrefine executable
initial_csv = './openrefine/v0.1/IHEC_metadata_harmonization.v0.1.csv'  # csv to build project from

old_project_name = os.path.splitext(os.path.basename(initial_csv))[0]
run([openrefine_client, '--create', initial_csv], check=True)
run([openrefine_client, '--apply', rule_file, old_project_name], check=True)

new_csv = "./openrefine/v0.3/IHEC_metadata_harmonization.v0.3.csv"
run([openrefine_client, '--export', f'--output={new_csv}',
     old_project_name], check=True)

# run([openrefine_client, '--create', new_csv], check=True)
