"""Command-line entrypoint for the PICurv conductor."""

from .cli import build_main_parser, dispatch_command


def main():
    """!
    @brief Parse command-line arguments and dispatch the requested command.
    """
    dispatch_command(build_main_parser().parse_args())
