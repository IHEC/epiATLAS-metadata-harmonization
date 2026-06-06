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
        
        fields_to_include = get_fields_to_include()
        # Convert each row into a dictionary and add it to data{}
        for row in csvReader:
            # Remove "harmonized_" from key names
            newRow = {}
            for key, value in row.items():
                # newRow[key] = value
                # # Include donor_id (without 'harmonized_') duplicate, for IDP legacy compatibility.
                # if key.startswith("harmonized_donor_id"):
                #     newRow['donor_id'] = value
                if key in fields_to_include:
                    # The data row, without "harmonized_" in all the field names.
                    if key.startswith("harmonized_"):
                        newRow[key.replace('harmonized_', '')] = value
                    # But keep the automated_ in the key name.
                    elif "_harmonized" in key:
                        newRow[key.replace('_harmonized', '')] = value    
                    else:
                        newRow[key] = value
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

# Not all sample metadata from the _extended .csv are necessary.  Only include these ones.
def get_fields_to_include():
    fields_to_include = [
        "EpiRR",
        "project",
        "harmonized_biomaterial_type",
        "harmonized_sample_label",
        "harmonized_sample_ontology_intermediate",
        "harmonized_sample_disease_high",
        "harmonized_sample_disease_intermediate",
        "harmonized_EpiRR_status",
        "epiATLAS_status",
        "harmonized_cell_type",
        "harmonized_cell_line",
        "harmonized_tissue_type",
        "harmonized_sample_ontology_curie",
        "harmonized_cell_markers",
        "automated_harmonized_sample_ontology",
        "automated_harmonized_sample_ontology_term",
        "harmonized_sample_disease",
        "harmonized_sample_disease_ontology_curie",
        "automated_harmonized_sample_disease_ontology_curie_ncit",
        "harmonized_donor_type",
        "harmonized_donor_id",
        "harmonized_donor_age",
        "harmonized_donor_age_unit",
        "automated_harmonized_donor_age_in_years",
        "harmonized_donor_life_stage",
        "harmonized_donor_life_stage_uncorrected",
        "harmonized_donor_sex",
        "harmonized_donor_sex_uncorrected",
        "harmonized_donor_health_status",
        "harmonized_donor_health_status_ontology_curie",
        "automated_harmonized_donor_health_status_ontology_curie_ncit",
        "epirr_id_without_version",
    ]
    return fields_to_include

def main():
    # Define in and out file paths, relative to this script.
    # sampleCsvFile = r'../openrefine/v1.2/IHEC_metadata_harmonization.v1.2.csv'
    sampleCsvFile = r'../openrefine/v2.0/IHEC_sample_metadata_harmonization.v2.0_extended.csv'
    datasetCsvFile = r'../openrefine/v1.1/epiatlas_metadata.csv' # No update found.
    # hubOutFilePath = r'IHEC_metadata_harmonization.v1.2.json'
    hubOutFilePath = r'IHEC_metadata_harmonization.v2.0.json'
    
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
