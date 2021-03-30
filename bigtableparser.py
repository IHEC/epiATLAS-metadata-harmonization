import pandas as pd
import json

def to_json(e):
	return json.dumps(e, indent=4, sort_keys=True, default=str)


def to_jsonfile(f, data):
	with open(f, "w") as outfile:
		outfile.write(to_json(data))
	return f

def uniq(es):
	if len(es) == 1:
		for e in es:
			return e
	else: raise Exception(es)

def hash_bigtable(f=None):
	if not f: f = "raw/IHEC_metadata_summary.xlsx"
	frame = pd.read_excel(f, sheet_name="BigTable")
	data = frame.to_dict()
	columns = data.keys()
	n = uniq({len(data[column]) for column in columns})
	hashed = [{column: data[column][i] for column in columns} for i in range(n)] 
	jsoned = to_json(hashed) # clean up any unknown types
	parsed = json.loads(jsoned)
	return to_jsonfile(f + ".json", parsed)


if __name__ == "__main__":
	print(hash_bigtable())



