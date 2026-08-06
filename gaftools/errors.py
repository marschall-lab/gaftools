"""
All the expected errors with appropriate names.
"""


class CommandLineError(Exception):
    """Super class for remoting errors in the CLI"""


class IncorrectHaplotypeError(CommandLineError):
    """Incorrect format of haplotype. Should be H followed by a number or 'none'"""


class IncorrectPhaseSetError(CommandLineError):
    """Incorrect format of phase set. Should be a number or 'none'"""


class DuplicateHaplotagError(CommandLineError):
    """Multiple instances of the same read in the haplotag file."""


class IncompatibleHaplotagError(CommandLineError):
    """Incompatible haplotype and phaseset information.
    Either both are 'none' or have appropriate values."""


class IncorrectGfaFormatError(CommandLineError):
    """GFA file with formatting errors"""


class IncorrectGafFormatError(CommandLineError):
    """GAF file with formatting errors"""


class PathNotFoundError(CommandLineError):
    """GFA does not have the path specified"""


class IncorrectPathFormatError(CommandLineError):
    """Path format provided is incorrect"""


class ReadNotFoundError(CommandLineError):
    """Query read not found in GAF"""


class IncorrectReadFormatError(CommandLineError):
    """Sequence file is neither FASTA or FASTQ"""


class IndexNotFoundError(CommandLineError):
    """Index file not found"""


class ChromosomeNotFoundError(CommandLineError):
    """Chromosome not found in the list of components identified"""


class BranchedGfaComponentError(CommandLineError):
    """A GFA component has multiple branches"""
