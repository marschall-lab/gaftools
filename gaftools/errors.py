"""
All the expected errors with appropriate names.
"""


class CommandLineError(Exception):
    """Super class for remoting errors in the CLI"""

    pass


class IncorrectHaplotypeError(CommandLineError):
    """Incorrect format of haplotype. Should be H followed by a number or 'none'"""

    pass


class IncorrectPhaseSetError(CommandLineError):
    """Incorrect format of phase set. Should be a number or 'none'"""

    pass


class DuplicateHaplotagError(CommandLineError):
    """Multiple instances of the same read in the haplotag file."""

    pass


class IncompatibleHaplotagError(CommandLineError):
    """Incompatible haplotype and phaseset information.
    Either both are 'none' or have appropriate values."""

    pass


class IncorrectGfaFormatError(CommandLineError):
    """GFA file with formatting errors"""

    pass


class IncorrectGafFormatError(CommandLineError):
    """GAF file with formatting errors"""

    pass


class PathNotFoundError(CommandLineError):
    """GFA does not have the path specified"""

    pass


class IncorrectPathFormatError(CommandLineError):
    """Path format provided is incorrect"""

    pass


class ReadNotFoundError(CommandLineError):
    """Query read not found in GAF"""

    pass


class IncorrectReadFormatError(CommandLineError):
    """Sequence file is neither FASTA or FASTQ"""

    pass


class IndexNotFoundError(CommandLineError):
    """Index file not found"""

    pass


class ChromosomeNotFoundError(CommandLineError):
    """Chromosome not found in the list of components identified"""

    pass


class BranchedGfaComponentError(CommandLineError):
    """A GFA component has multiple branches"""

    pass
