Keeping track of the files being used for gfa2rgfa testing

# FASTA

assembly-*.fa
* These are the fasta files for the dummy assemblies (along with their corresponding indices which will be generated).

# GFA

graph.gfa
* This is the base GFA used for all the conversions to rGFA.
* There are some tweaked forms of this graph.

graph.gfa.gz
* BGZip version of graph.gfa

graph-invertednodes.gfa
* some nodes (eg. s24) are present in their reverse complement form.
* used to test if the gfa2rgfa will automatically fix them.

graph-partial-tagged.gfa
* the reference nodes of graph.gfa have been tagged.
* this is the standard form in which minigraph-cactus outputs its graph.

graph-partial-tagged-alternateSN.gfa
* differs to graph-partial-tagged.gfa because its tagged nodes have SN tags of form `REF#CONTIG1` instead of the standard `REF#0#CONTIG1`
* used to test if the gfa2rgfa fixes them back to standard form.

graph-wrong-Wline-order.gfa
* the W line order has been mixed up and seqfile provided to see if it still produces correct output.

# rGFA

reference-graph.gfa
* rGFA produced based on GFA and order of W-lines present in it.

reference-graph-FOO#1.gfa
* rGFA where the reference is FOO#1.

reference-graph-seqfile.gfa
* rGFA where the seqfile is used for custom assembly order.

reference-graph-partial-seqfile.gfa
* rGFA where the seqfile consists a subset of assemblies


# seqfile

samples.seqfile
* seqfile with all the assemblies

samples-partial.seqfile
* seqfile with a subset of assemblies
