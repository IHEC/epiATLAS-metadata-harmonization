#!/usr/bin/python3
# Converts IHEC_metadata_hormonization.vX.X.csv file to json format, with EpiRR used as the primary key.
# This was the initial script version that evolved into ./makeIAhub.py
import csv
import json
 
# Function to convert a CSV to JSON
# Takes the file paths as arguments
def make_json(csvFilePath, jsonFilePath):
     
    # create a dictionary
    data = {}
     
    # Open a csv reader called DictReader
    with open(csvFilePath, encoding='utf-8') as csvf:
        csvReader = csv.DictReader(csvf)
         
        # Convert each row into a dictionary and add it to data
        for row in csvReader:
            # Remove "harmonized_" from key names
            newRow = {} # The data row, without "harmonized_" in all the key names.
            for key, value in row.items():
                if key.startswith("harmonized_"):
                    newRow[key.replace('harmonized_', '')] = value
                else:
                    newRow[key] = value
            # Column 'EpiRR' is the primary key
            key = newRow['EpiRR']
            data[key] = newRow
            # Remove duplicate EpiRR key entry
            del data[key]['EpiRR']
            
    # Open a json writer, and use the json.dumps() function to dump data
    with open(jsonFilePath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))
         
def main():
    # Define in and out file paths, relative to this file.
    csvFilePath = r'../openrefine/v1.1/IHEC_metadata_harmonization.v1.1.csv'
    jsonFilePath = r'asJSON.json'
    
    # execute the function
    make_json(csvFilePath, jsonFilePath)

if __name__ == "__main__":
  main()
