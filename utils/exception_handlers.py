"""exception_handlers.py: module for handling exceptions."""
import logging
from utils.util import multi_log


def argument_exception(error_args, logger=None):
    """
    Handles mistakes in key arguments by letting the user re-insert the argument.
    :param error_args: the argument's name that had the mistake. Str
    :param logger: must log the mistake and the new argument. logging object
    :return: the new argument
    """
    arg_new = input("Argument {0} was not correct, please submit another argument: ".format(error_args))
    if logger:
        multi_log(logger, "Argument {0} was not entered correctly".format(error_args),
                  "New argument {0}: {1}".format(error_args, arg_new))
    return arg_new
