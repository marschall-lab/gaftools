import os
import sys
import glob


if len(sys.argv) < 2:
	print("You need to give a directory where the graph components are")
	sys.exit()

if not sys.argv[1].endswith("/"):
	sys.argv[1] += "/"
for f in glob.glob(f"{sys.argv[1]}*.gfa"):
	sn_name_counts = dict()
	with open(f, "r") as infile:
		for l in infile:
			if l.startswith("S"):
				if "SN" in l:
					l = l.strip().split()
					i = 0
					for idx, c in enumerate(l):
						if c.startswith("SN"):
							i = idx
					sn_tag = l[i].split(":")[-1]
					if sn_tag in sn_name_counts:
						sn_name_counts[sn_tag] += 1
					else:
						sn_name_counts[sn_tag] = 1

	# Taking a majority vote
	majority = (0, "")
	for sn, count in sn_name_counts.items():
		if count > majority[0]:
			majority = (count, sn)
	os.rename(f, os.path.join(sys.argv[1], f"{majority[1]}.gfa"))
