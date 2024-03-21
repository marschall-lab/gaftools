#!/bin/bash

# This script takes a GFA and split it into its individual components
# to be able to split the graph, you need to have GFASubgraph package installed
# This can be found here https://github.com/fawaz-dabbaghieh/gfa_subgraphs/

# then the script runs the renaming script to rename each component


in_gfa=$1
out_dir=$2

if [ -z "$1" ]
then
	echo "You need to give the input GFA as first argument"
	exit 1
fi


if [ -z "$2" ]
then
	echo "You need to give the output directory as second argument"
	exit 1
fi

if [ ! -f $1 ]
then
    echo "The GFA" $1 "was not found"
    exit
fi


echo "separating the components from the graph " $1
echo `date`

/usr/bin/time -v GFASubgraph -g $1 output_comps --output_dir $2 2> splitting_time.log

echo "Done separatig the components"
echo `date`



echo "Now, renaming the components"
echo `date`

echo 'import os
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
	os.rename(f, os.path.join(sys.argv[1], f"{majority[1]}.gfa"))' | /usr/bin/time -v python3 - $2 2> renaming_time.log

echo "Finished renaming the components"
echo `date`
