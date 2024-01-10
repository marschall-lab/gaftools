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

GFASubgraph -g $1 output_comps --output_dir $2

echo "Done separatig the components"
echo `date`



echo "Now, renaming the components"
echo `date`

python3 rename_comps.py $2

echo "Finished renaming the components"
echo `date`