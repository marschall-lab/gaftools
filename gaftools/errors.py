"""
All the expected errors with appropriate names.
"""


class CommandLineError(Exception):
    pass


class IncorrectHaplotypeError(CommandLineError):
    pass


class IncorrectPhaseSetError(CommandLineError):
    pass


class DuplicateHaplotagError(CommandLineError):
    pass


class IncompatibleHaplotagError(CommandLineError):
    pass


class IncorrectGfaFormatError(CommandLineError):
    pass
