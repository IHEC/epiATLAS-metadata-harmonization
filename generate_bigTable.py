#!/usr/bin/env python3

"""
This script is used to generate metadata for EpiAtlas data
based on EpiRR registeries available from https://www.ebi.ac.uk/vg/epirr/view/all?format=json.
It generates two files:
1. EpiAtlas_EpiRR.txt (long format)
2. EpiAtlas_EpiRR_metadata_all.csv (wide format: reshaped from EpiAtlas_EpiRR.txt file)
NaN was produced by pivot_table when reshaping long to wide foramt (i.e. the value of a specific metadata item was not reported originally by the project).
Na means the value of a specific metadata item was reported by the project as "NA"
for questions please contact Abdulrahman Salhab: abdulrahman.salhab@uni-saarland.de
"""


import urllib.request, json
import pandas as pd

url="https://www.ebi.ac.uk/vg/epirr/view/all?format=json"
response = urllib.request.urlopen(url)
data = json.loads(response.read())


d = []

fo = open('raw/EpiAtlas_EpiRR.txt', 'w')
print('EpiRR'+'\t'+'EpiRR_status'+'\t'+'project'+'\t'+'metadata'+'\t'+'value', file=fo)

for idx in range(0,len(data)):
    url=data[idx]["_links"]["self"]
    response = urllib.request.urlopen(url)
    url_json = json.loads(response.read())
    for key,value in url_json["meta_data"].items():
        print(url_json["full_accession"]+'\t'+url_json["status"]+'\t'+url_json["project"]+'\t'+key+'\t'+value, file=fo)
        d.append([url_json["full_accession"],url_json["status"],url_json["project"],key,value])

fo.close()


df = pd.DataFrame(d)
df.columns = ['EpiRR', 'EpiRR_status', 'project', 'metadata', 'value']
df2 = df.pivot_table(index=["EpiRR","EpiRR_status","project"], columns='metadata', values='value', aggfunc='first')
df2.to_csv(r'raw/EpiAtlas_EpiRR_metadata_all.csv', na_rep="NaN")

