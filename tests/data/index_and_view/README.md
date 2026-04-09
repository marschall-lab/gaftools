Keeping track of the files being used for view and index testing

Note:

* We do not test conversion validity in these test cases.
* We test subsetting of file and if the conversion works.
* The ultimate conversion test happens using real data and with the two modes of minigraph.

Note: the graphs used to align these reads are `../gfa2rgfa/reference-graph.gfa` and `../gfa2rgfa/reference-graph-FOO#2.gfa`

# FASTA

reads.fa

* all the reads for the test cases.

# GAF

graphaligner.gaf(.gz)

* Aligning `reads.fa` to `../gfa2rgfa/reference-graph.gfa` using GraphAligner (and its gzip version)

graphaligner-stable.gaf(.gz)

* Converted graphaligner.gaf to stable format using gaftools view
* NOTE: The ultimate conversion test happens using real data and with the two modes of minigraph.

# Viewing Index files

get_offsets.py

* Given a file, it outputs the offsets of the line starting for all the lines (the same way gaftools index finds it).
* These offsets are used for make_indices.py.

make_indices.py

* simple python script which makes the pickled viewing indices using the information of offsets from get_offsets.py

view-index-stable.gvi

* View index for graphaligner-stable.gaf(.gz)

view-index-unstable.gvi

* View index for graphaligner.gaf(.gz)
