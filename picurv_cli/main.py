"""Command-line entrypoint for the PICurv conductor."""

import argparse

from .cli import build_main_parser, dispatch_command


def _selected_subparser(parser, args):
    """!
    @brief The deepest sub-command parser the parsed arguments selected.
    @param[in] parser Top-level parser.
    @param[in] args Namespace from parsing with @p parser.
    @return The parser whose usage describes the selected command.
    """
    while True:
        actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
        name = getattr(args, actions[0].dest, None) if actions else None
        if name not in (actions[0].choices if actions else {}):
            return parser
        parser = actions[0].choices[name]


def main():
    """!
    @brief Parse command-line arguments and dispatch the requested command.
    @details Commands that forward arguments to make take its options after `--`; an option
             given before it is reported by the command it was given to, with that form.
    """
    parser = build_main_parser()
    args, unknown = parser.parse_known_args()
    if unknown:
        message = f"unrecognized arguments: {' '.join(unknown)}"
        if hasattr(args, "make_args"):
            message += f"; pass make options after '--', as in 'picurv {args.command} -- {unknown[0]}'"
        _selected_subparser(parser, args).error(message)
    if getattr(args, "make_args", None) and args.make_args[0] == "--":
        args.make_args = args.make_args[1:]
    dispatch_command(args)
