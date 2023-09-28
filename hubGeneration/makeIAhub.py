#!/usr/bin/python3
# Convert the sample metadata input .csv to an IHEC data hub.
# File paths are hard coded, for now.
import csv
import json
import datetime
import re

class Hub:
    def __init__(self):
        self.data = {
            "datasets": {},
            "hub_description": {},
            "samples": {},
            }

    def jsonify(self):
        # Legible and pretty representation
        return json.dumps(self.data, indent=4, sort_keys=True)

# TODO: Get official values
class Hub_Description:
    def __init__(self):
        now = datetime.datetime.now()
        # now_formatted = now.strftime("%Y-%m-%d")
        self.h_d = {
            "name": "EpiAtlas Data Hub", # New field.
            "assembly": "TBD", 
            "date": "Data Freeze", #TODO
            "description": "IHEC EpiAtlas", 
            "email": "TBD", 
            "publishing_group": "TBD", 
            "releasing_group": "TBD", 
            "taxon_id": 9606
        }

# Converts a Sample metadata csv to an IDP Sample hub.
def makeSamplesJson(sampleCsvFile):
     
    # create a dictionary which will become the hub
    data = {}
     
    # Open a DictReader
    with open(sampleCsvFile, encoding='utf-8') as csvf:
        csvReader = csv.DictReader(csvf)
         
        # Convert each row into a dictionary and add it to data{}
        for row in csvReader:
            # Remove "harmonized_" from key names
            newRow = {}
            for key, value in row.items():
                newRow[key] = value
            # # The data row, without "harmonized_" in all the field names.
            #     if key.startswith("harmonized_"):
            #         newRow[key.replace('harmonized_', '')] = value
            #     else:
            #         newRow[key] = value

            # Column 'EpiRR' is used as the primary key
            key = newRow['EpiRR']
            data[key] = newRow
            # Remove duplicate EpiRR key entry
            del data[key]['EpiRR']
            
    return data

def makeDatasetJson(datasetCsvFile):
    data = {}

    # Open a DictReader
    with open(datasetCsvFile) as csvf:
        csvReader = csv.DictReader(csvf)
        for row in csvReader:
            key = row['epirr_id'] + "_" + row['assay_type'] + "_" + row['experiment_type']
            blankDict = {# "analysis_attributes": {},
                         "browser": {},
                         "experiment_attributes": {}
                         }
            data[key] = blankDict
            # Parse out analysis software
            match = re.search(r'[^v]*', row['software_version'])
            anal_soft = match.group()
            # Parse out analysis software version
            ver_match = re.search(r'[\d\.]+', row['software_version'])
            anal_soft_v = ver_match.group()
            analysis_attributes = {
                "analysis_group": "IHEC IA - TBD",
                # "alignment_software": "TBD", # Was there any alignment, or is it done before the IA?
                "alignment_software_version": "TBD",
                "analysis_software": anal_soft,
                "analysis_software_version": anal_soft_v
            }
            data[key]['analysis_attributes'] = analysis_attributes
            # Link to a sample.
            data[key]["sample_id"] = row['epirr_id']
    return data

def main():
    # Define in and out file paths, relative to this script.
    sampleCsvFile = r'../openrefine/v1.1/IHEC_metadata_harmonization.v1.1.csv'
    datasetCsvFile = r'../openrefine/v1.1/epiatlas_metadata.csv'
    hubOutFilePath = r'IA_hub.json'
    
    # Put the pieces together
    h = Hub()
    h.data["hub_description"] = Hub_Description().h_d
    h.data['samples'] = makeSamplesJson(sampleCsvFile)
    h.data['datasets'] = makeDatasetJson(datasetCsvFile)

    # Write to a file.
    with open(hubOutFilePath, 'w') as jsonf:
        jsonf.write(h.jsonify())

if __name__ == "__main__":
  main()
