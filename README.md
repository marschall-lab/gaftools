# gaftools

General purpose utility related to GAF files

## Installation

### Normal Installation

Follow these instructions to install gaftools.

```sh
git clone git@github.com:marschall-lab/gaftools.git
cd gaftools
pip install -e .
```

### Developer Installation

Using these steps you can install gaftools in a virtual environment.

#### Using Python Virtual Environment

```sh
git clone git@github.com:marschall-lab/gaftools.git
cd gaftools
python -m venv venv
source venv/bin/activate
pip install -e .
```

#### Using Conda Environment

```sh
conda create -n gaftools-dev python=3.10
conda activate gaftools-dev
pip install -e .
```

In order to use gaftools, the virtual environment has to be activated with `source venv/bin/activate`.

## Usage

gaftools has a command-line interface which can be accessed after installed. There are multiple subcommands that are available for the user. Below, the help window is shown.

```sh
$ gaftools --help
usage: gaftools [-h] [--version] [--debug] {convert,index,order_gfa,phase,realign,sort,stat,view} ...

positional arguments:
  {convert,index,order_gfa,phase,realign,sort,stat,view}
    convert             Convert Coordinate Systems between the stable system and unstable system
    index               Index the GAF File
    order_gfa           Adds BO and NO tags to GFA
    phase               Add phasing information to the bam file
    realign             Realign GAF file using wavefront alignment algorithm (WFA)
    sort                Sorting GAF alignments using BO and NO tags of the corresponding graph
    stat                Statistics of a GAF File
    view                View the GAF File based on parameters

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --debug               Print debug messages
```

The following subcommands are currently available:

1. [convert](#convert): This command can be used to convert the formatting of GAF files. The GAF file can be converted from the vertex-based formatting (from here on, this will be referred to as unstable-coordinated based formatting) to stable-coordinate based formatting.

2. [index](#index): This command indexes the GAF file which is required for viewing. Currently, an inverted search list is being used to access parts of the GAF file which contain particular nodes.

3. [order_gfa](#order_gfa): This commands tags the GFA using the BO and NO notations required for sorting.

4. [phase](#phase): This command adds WhatsHap-generated phase tags to the GAF.

5. [realign](#realign): This command realigns the GAF alignments using wavefront alignment algorithm.

6. [sort](#sort): This command sorts the GAF file based on the BO and NO tagging scheme.

7. [stat](#stat): This command produces certain stats for the GAF file.

8. [view](#view): This command allows users to view the GAF under different filters. The user can give nodes/regions in the genome and command lists out the alignments that involve those nodes/regions. There are many viewing formats that can specified to tailor the output. By default, the viewing output shows the nodes in a human readable form (instead of node IDs, it shows the contig, start and end position associated with the node). But note that this human readable format is not the stable coordinate system.

### <a id="convert"></a> Subcommand: convert

```sh
$ gaftools convert --help
usage: gaftools convert [-h] [-o OUTPUT] [--unstable] [--stable] GAF rGFA

Convert Coordinate Systems between the stable system and unstable system

positional arguments:
  GAF                   GAF File whose coordinates have to be changed
  rGFA                  Input rGFA file to convert the coordinates

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output GAF file. If omitted, use standard output.
  --unstable            Convert to Unstable Coordinates
  --stable              Convert to Stable Coordinates

```

The `convert` command requires two positional arguments: the path to the GAF file and the path to the rGFA file (which was used for the alignment).

Options:

--unstable - Provide this flag when the GAF has to be converted to unstable coordinate system.

--stable - Provide this flag when the GAF has to be converted to stable coordinate system.

### <a id="index"></a> Subcommand: index

```sh
$ gaftools index --help
usage: gaftools index [-h] [-o OUTPUT] GAF rGFA

Index the GAF File

positional arguments:
  GAF                   Input GAF file (can be gzip-compressed)
  rGFA                  Reference rGFA file has to be input.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output Indexed GAF file. If omitted, use <GAF File>.gai.
```

The `index` command produces a file which has information of the positions of nodes in the GAF file. The same index file should work for both the stable and unstable GAF files.

### <a id="order_gfa"></a> Subcommand: order_gfa

```sh
$ gaftools order_gfa --help
usage: gaftools order_gfa [-h] [--chromosome_order CHROMOSOME_ORDER] [--with-sequence] [--outdir OUTDIR] GRAPH

Adds BO and NO tags to GFA

positional arguments:
  GRAPH                 Input GFA file

optional arguments:
  -h, --help            show this help message and exit
  --chromosome_order CHROMOSOME_ORDER
                        Order in which to arrange chromosomes in terms of BO sorting. Expecting comma-separated list. Default: chr1,...,chr22,chrX,chrY,chrM
  --with-sequence       Retain sequences in output (default is to strip sequences)
  --outdir OUTDIR       Output Directory to store all the GFA and CSV files. Default location is a "out" folder from the directory of execution.
```

### <a id="phase"></a> Subcommand: phase

```sh
$ gaftools phase --help
usage: gaftools phase [-h] [-o OUTPUT] GAF phase

Add phasing information to the bam file

positional arguments:
  GAF                   GAF File
  phase                 WhatsHap Haplotag file (.tsv)

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output GAF file. If omitted, use standard output.
```

### <a id="realign"></a> Subcommand: realign

```sh
$ gaftools realign --help
usage: gaftools realign [-h] GAF GFA FASTA

Realign GAF file using wavefront alignment algorithm (WFA)

positional arguments:
  GAF         GAF File
  GFA         Input GFA file
  FASTA       Input FASTA file

optional arguments:
  -h, --help  show this help message and exit
```

### <a id="sort"></a> Subcommand: sort

```sh
$ gaftools sort --help
usage: gaftools sort [-h] [--outgaf OUTGAF] [--outind OUTIND] [--bgzip] GAF GFA

Sorting GAF alignments using BO and NO tags of the corresponding graph

positional arguments:
  GAF              Input GAF File
  GFA              GFA file with the sort keys (BO and NO tagged)

optional arguments:
  -h, --help       show this help message and exit
  --outgaf OUTGAF  Output GAF File path (Default: sys.stdout)
  --outind OUTIND  Output Index File path for the GAF file. (When --outgaf is not given, no index is created. If it is given and --outind is not specified, it will have same file name with
                   .gai extension)
  --bgzip          Flag to bgzip the output. Can only be given with --output.
```

### <a id="stat"></a> Subcommand: stat

```sh
$ gaftools stat --help
usage: gaftools stat [-h] [--cigar] GAF

Statistics of a GAF File

positional arguments:
  GAF         Input GAF file

optional arguments:
  -h, --help  show this help message and exit
  --cigar     Outputs cigar related statistics (takes a long)
```

### <a id="view"></a> Subcommand: view

```sh
$ gaftools view --help
usage: gaftools view [-h] [-i INDEX] [-o OUTPUT] [--node NODE] [--only-alignment] [--full-alignment] [--remove-cigar] [--remove-read-id] [--show-node-id] GAF

View the GAF File based on parameters

positional arguments:
  GAF                   Input GAF file (can be gzip-compressed)

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        Path to GAF Index file. If not provided, it is assumed to be in the same directory as GAF file with the same name and .gaf.gai extension
  -o OUTPUT, --output OUTPUT
                        Output Indexed GAF file. If omitted, use <GAF File>.gai.
  --node NODE           Specify nodes to filter alignments. Instead of node ID, regions can also be specified. Can be used multiple times. When multiple nodes are specified, output
                        contains partial alignment between first and last node. Entire alignment can be shown with --full-alignment flag.
  --only-alignment      Show alignments which contain the list of nodes.
  --full-alignment      Show the entire alignment with multiple nodes given.
  --remove-cigar        Option to remove CIGAR strings (if --only-alignment has not been chosen).
  --remove-read-id      Option to remove read IDs from output.
  --show-node-id        Show list of nodes as node IDs. Option only available with --only-alignment (Default: Show in human readable form.)
```

The `view` command requires an index file to be provided along with the GAF file. It has the following formatting/filtering options:

--node NODE -> The user can provide node IDs/regions and the command outputs the lines which have an intersection of all the node IDs/regions specified. The regions can be provided using the format [contig name]:[start pos]-[end pos] or simply [contig name]:[position]. The --node option can be given multiple times. 

If two or more nodes have been specified, or the region that has been specified spans multiple nodes, then the output will contain a part of the path matching (or alignment) between the first encountered node and last encountered node given by the user. The full path matching (or alignment) can be given in the output using the --full-alignment flag. 

--only-alignment -> The user can give this flag to restrict the out to just the path matching column and the read ID (this can be removed with --remove-read-id flag). By default, the entire line is otherwise shown in the output.

--remove-cigar -> When --only-alignment flag has not been provided, this flag can be provided to omit the cigar string from the output.

--remove-read-id -> With the --only-alignment or without it, the read ID is given as the first column of the output. This flag removes that column.

--show-node-id -> By default, the alignments or path matching is shown in a human-readable format (This is not the stable coordinates exactly. It just makes it easier to read the alignment). Through this flag, instead of the human-readable format, the node IDs are instead shown in the output.

## Editing

Since no C++ code has been implemented yet (the .cpp scripts in src directory are not utilised anywhere), you can simply edit the python scripts and then run the command. When C++ scripts get involved, to make edits and see the changes, the `pip install -e .` command has to be run again.

