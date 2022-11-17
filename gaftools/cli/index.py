"""
Index the GAF File
"""

import logging
import sys
import platform

from gaftools import __version__
from gaftools.cli import log_memory_usage
from gaftools.cli import CommandLineError


logger = logging.getLogger(__name__)

def run(
    phase_input_files,
    variant_file,
    reference=None,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    ignore_read_groups=False,
    indels=True,
    mapping_quality=20,
    max_coverage=15,
    nopriors=False,
    ped=None,
    recombrate=1.26,
    genmap=None,
    gt_qual_threshold=0,
    prioroutput=None,
    constant=0.0,
    overhang=10,
    affine_gap=False,
    gap_start=10,
    gap_extend=7,
    mismatch=15,
    write_command_line_header=True,
    use_ped_samples=False,
    eff_pop_size = 10
):
    exit


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be gzip-compressed)')
    arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
        help='BAM or VCF file(s) with phase information, either through sequencing reads (BAM) or through phased blocks (VCF)')

    arg('-o', '--output', default=sys.stdout,
        help='Output VCF file. Add .gz to the file name to get compressed output. '
        'If omitted, use standard output.')
    arg('--reference', '-r', metavar='FASTA',
        help='Reference file. Provide this to detect alleles through re-alignment. '
        'If no index (.fai) exists, it will be created')

    arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
    arg('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
        help='Reduce coverage to at most MAXCOV (default: %(default)s).')
    arg('--mapping-quality', '--mapq', metavar='QUAL',
        default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
    arg('--indels', dest='indels', default=False, action='store_true',
        help='Also genotype indels (default: genotype only SNPs)')
    arg('--ignore-read-groups', default=False, action='store_true',
        help='Ignore read groups in BAM header and assume all reads come '
        'from the same sample.')
    arg('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
        help='Name of a sample to genotype. If not given, all samples in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
        help='Name of chromosome to genotyped. If not given, all chromosomes in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--gt-qual-threshold', metavar='GTQUALTHRESHOLD', type=float, default=0,
        help='Phred scaled error probability threshold used for genotyping (default: %(default)s). Must be at least 0. '
        'If error probability of genotype is higher, genotype ./. is output.')
    arg('--no-priors', dest='nopriors', default=False, action='store_true',
        help='Skip initial prior genotyping and use uniform priors (default: %(default)s).')
    arg('-p', '--prioroutput', default=None,
        help='output prior genotype likelihoods to the given file.')
    arg('--overhang', metavar='OVERHANG', default=10, type=int,
        help='When --reference is used, extend alignment by this many bases to left and right when realigning (default: %(default)s).')
    arg('--constant', metavar='CONSTANT', default=0, type=float,
        help='This constant is used to regularize the priors (default: %(default)s).')
    arg('--affine-gap', default=False, action='store_true',
        help='When detecting alleles through re-alignment, use affine gap costs (EXPERIMENTAL).')
    arg('--gap-start', metavar='GAPSTART', default=10, type=float,
        help='gap starting penalty in case affine gap costs are used (default: %(default)s).')
    arg('--gap-extend', metavar='GAPEXTEND', default=7, type=float,
        help='gap extend penalty in case affine gap costs are used (default: %(default)s).')
    arg('--mismatch', metavar='MISMATCH', default=15, type=float,
        help='mismatch cost in case affine gap costs are used (default: %(default)s)')
    arg('--eff-pop-size', metavar='EFFPOPSIZE', default = 10, type = int,
        help="Parameter for transition probability computing (default: %(default)s)")

    arg = parser.add_argument_group('Pedigree genotyping').add_argument
    arg('--ped', metavar='PED/FAM',
        help='Use pedigree information in PED file to improve genotyping '
        '(switches to PedMEC algorithm). Columns 2, 3, 4 must refer to child, '
        'mother, and father sample names as used in the VCF and BAM. Other '
        'columns are ignored (EXPERIMENTAL).')
    arg('--recombrate', metavar='RECOMBRATE', type=float, default=1.26,
        help='Recombination rate in cM/Mb (used with --ped). If given, a constant recombination '
        'rate is assumed (default: %(default)gcM/Mb).')
    arg('--genmap', metavar='FILE',
        help='File with genetic map (used with --ped) to be used instead of constant recombination '
        'rate, i.e. overrides option --recombrate.')
    arg('--use-ped-samples', dest='use_ped_samples',
        action='store_true', default=False,
        help='Only work on samples mentioned in the provided PED file.')
    
# fmt: on


def validate(args, parser):
    if args.ignore_read_groups and args.ped:
        parser.error("Option --ignore-read-groups cannot be used together with --ped")
    if args.genmap and not args.ped:
        parser.error("Option --genmap can only be used together with --ped")
    if args.genmap and len(args.chromosomes) != 1:
        parser.error(
            "Option --genmap can only be used when working on exactly one chromosome (use --chromosome)"
        )
    if len(args.phase_input_files) == 0:
        parser.error("Not providing any PHASEINPUT files not allowed for genotyping.")
    if args.gt_qual_threshold < 0:
        parser.error("Genotype quality threshold (gt-qual-threshold) must be at least 0.")
    if args.prioroutput is not None and args.nopriors:
        parser.error("Genotype priors are only computed if --no-priors is NOT set.")
    if args.constant != 0 and args.nopriors:
        parser.error("--constant can only be used if --no-priors is NOT set..")
    if args.affine_gap and not args.reference:
        parser.error("Option --affine-gap can only be used together with --reference.")
    if args.use_ped_samples and not args.ped:
        parser.error("Option --use-ped-samples can only be used when PED file is provided (--ped).")
    if args.use_ped_samples and args.samples:
        parser.error("Option --use-ped-samples cannot be used together with --samples")


def main(args):
    run(**vars(args))
