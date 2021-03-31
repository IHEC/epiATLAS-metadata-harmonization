import pandas as pd
import json
import sys
import merging

uniq = merging.uniq


def to_json(e):
	return json.dumps(e, indent=4, sort_keys=True, default=str)


def to_jsonfile(f, data):
	with open(f, "w") as outfile:
		outfile.write(to_json(data))
	return f


def from_jsonfile(f):
	with open(f) as infile:
		return json.load(infile)



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



def merge_attributes(f=None):
	if not f: f = "config/merges.json"
	dbfile = 'raw/IHEC_metadata_summary.xlsx.json' 
	cfg = from_jsonfile(f)
	records = from_jsonfile(dbfile)
	updated = list()
	minimal = list()
	for record in records:
		sm_record = dict()
		for rule in cfg:
			strategy = uniq(rule['strategy']) 
			opts = rule.get("options")
			f = getattr(merging, strategy)
			term = rule['harmonized']
			record[term] = f(record, rule['strategy'][strategy], opts)
			sm_record[term] = record[term]
		updated.append(record)
		minimal.append(sm_record)	

	print(to_jsonfile('merged/IHEC_metadata.merged.json', updated))
	print(to_jsonfile('merged/IHEC_metadata.merged_minimal.json', minimal))
	return 'ok'

if __name__ == "__main__":
	cfg = sys.argv[1:]
	if '-hash' in cfg:
		print(hash_bigtable())
	elif '-merge' in cfg:
		print(merge_attributes())
	else:
		print('...')



