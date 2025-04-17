Changes
=======

v1.1.2 (17/04/2025)
-------------------

* Better BGZF compatibility for all subcommands.
* Formatting updates with logger.
* :issue:`34`: Added subcommand :code:`gfa2rgfa` to convert GFA to rGFA format. Tested with HPRC Minigraph-Cactus Graph (v1.1).
* Updated documentation to include :code:`gfa2rgfa` subcommand.


v1.1.1 (23/01/2025)
-------------------

* Updated documentation to include conda installation instructions.
* Updated documentation to include preprint link and citation text.
* :issue:`33`: Fixed error when viewing regions spanning multiple nodes.
* :issue:`32`: Fixed error when trying to extract whole chromosomes.


v1.1.0 (09/10/2024)
-------------------

* Adding :code:`find_path` subcommand to take node path information as input and outputting the sequences.
* Updated :code:`order_gfa` subcommand to output chromosome-wise graphs only with :code:`--by-chrom` flag.
* Bug fix in :code:`view` subcommand and minor changes to make it faster.
* Updated documentation.


v1.0.0 (24/06/2024)
-------------------

* Initial commit.
