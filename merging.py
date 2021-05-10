import math
import re

curie_re = re.compile(r'([A-Za-z]+)[_:](\d+)$')


def uniq(es):
	if len(es) == 1:
		for e in es: return e
	else: raise Exception(es)

def try_int(e):
	x = e
	try: x = int(e)
	except:
		try: x = float(e)
		except: pass
	return x


def try_curie(e):
    m = re.search(curie_re, e)
    if m:
        return ':'.join([m.group(1).lower(), m.group(2)])
    else:
        return e


def merge_delimited(record, cfg, opts=None):
	def clean(e):
		if isinstance(e, float):
			if math.isnan(e): e = ''
		return str(e).strip()

	if not opts: opts = dict()
	delimit = opts.get("delimit", "::")
	found = [clean(record.get(e, '')) for e in cfg]
	for e in found:
		assert e.find(delimit) == -1


	cleanedup = {e for e in found if e}
	if 'lowercase' in opts:
		cleanedup = {e.lower() for e in cleanedup}
	elif 'make_int' in opts:
		cleanedup = {try_int(e) for e in cleanedup}
	elif 'merge_ontology' in opts:
		cleanedup = {try_curie(e) for e in cleanedup for e in e.split(';')}


	if len(cleanedup) == 1: return uniq(cleanedup)
	else:
		return delimit.join(cleanedup) if cleanedup else ''



