"""
The code has been taken from WhatsHap.
Link to WhatsHap: https://github.com/whatshap/whatshap
"""

import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter


class HelpfulArgumentParser(ArgumentParser):
    """An ArgumentParser that prints full help on errors."""

    def __init__(self, *args, **kwargs):
        if "formatter_class" not in kwargs:
            kwargs["formatter_class"] = RawDescriptionHelpFormatter
        super().__init__(*args, **kwargs)

    def error(self, message):
        self.print_help(sys.stderr)
        args = {"prog": "gaftools", "message": message}
        self.exit(2, "{prog} error: {message}\n".format(**args))
