
import json
from multiprocessing import Pool
from pytools.persistent_dict import PersistentDict as PD

patient = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"]



# def uniques(line)
# 	vs = []  
# 	all_entries = it.chain.from_iterable(line)
#     vs += [re.split(',|\*|\|', entry['tag'])[0] for entry in all_entries]
#     unique_vs = list(set(vs))
#     unique_vs.sort()

def file_opener(i):
	clones = ["1", "2", "3", "4", "5", "6"]
	vs = []
	for c in clones:
		with open('../data/input/%s_%s_clone.json' % (i,c)) as f:
	 		data = json.load(f)
	 		all_entries = it.chain.from_iterable(data)	
	 		vs += [re.split(',|\*|\|', entry['tag'])[0] for entry in all_entries]	
			unique_vs = list(set(vs))
    		rv = unique_vs.sort()	
   			return rv
	 		# for line in data:
	 		# 	p = Pool(64)
	   #          results = p.map(uniques, data)


for i in patient:
	result = file_opener(i)
	print(result)






	# for c in clones:
	# 	# print(type(i))
	# 	# print(type(c))
	# 	# print(i,c)
	# 	with open('../data/input/%s_%s_clone.json' % (i,c)) as f:
	# 		data = json.load(f)
	# 		for line in data:
	# 			p = Pool(64)
 #                results = p.map(line_worker, data)