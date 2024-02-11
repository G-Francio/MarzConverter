# Logging, see here: https://stackoverflow.com/questions/3383865/how-to-log-error-to-file-and-not-fail-on-exception
# and here: https://stackoverflow.com/questions/8050775/using-pythons-logging-module-to-log-all-exceptions-and-errors
# and here: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python/16993115#16993115

# And some other sources I forgot to note...

# Many thanks to this random post though
#   I came up with that after reading the answers to
#       http://stackoverflow.com/questions/5875225/
#   which pointed me to
#       http://bugs.python.org/issue6435

import logging
import sys
import tempfile
from pathlib import Path


class CachelessFormatter(logging.Formatter):
    """
    A custom logging formatter that disables the caching of the exception text.

    Inherits from logging.Formatter.
    """

    def format(self, record):
        """
        Format the log record, disabling the caching of the exception text.

        Parameters:
            record (LogRecord): The log record to format.

        Returns:
            str: The formatted log record.
        """
        # Disable the caching of the exception text.
        backup = record.exc_text
        record.exc_text = None
        s = logging.Formatter.format(self, record)
        record.exc_text = backup
        return s


class NoTracebackFormatter(CachelessFormatter):
    """
    A custom logging formatter that omits the traceback information when formatting exceptions.

    Inherits from CachelessFormatter.
    """

    def formatException(self, ei):
        """
        Format the exception information without including the traceback.

        Parameters:
            ei (tuple): The exception information tuple (type, value, traceback).

        Returns:
            str: The formatted exception information without the traceback.
        """
        return f"{ei[0].__name__} - {ei[1]}"


class FitsParseError(Exception):
    """
    Custom exception class for handling fits parsing errors.
    """


class NoInputError(Exception):
    """
    Custom exception class for handling cases where no input is provided.
    """


# Get logger
logger = logging.Logger("MyLogger")
logger.setLevel(logging.INFO)

# Clear handlers
if logger.hasHandlers():
    logger.handlers.clear()

logger.propagate = False

# Add Handlers
handler_file = logging.FileHandler(Path(tempfile.gettempdir()) / "MarzConverter.log")
handler_stream = logging.StreamHandler()

logger.addHandler(handler_file)
logger.addHandler(handler_stream)

handler_stream.setFormatter(NoTracebackFormatter())


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Custom exception handler function for logging uncaught exceptions.

    Parameters:
        exc_type (type): The type of the exception.
        exc_value (Exception): The exception instance.
        exc_traceback (traceback): The traceback information.

    Returns:
        None
    """
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.critical(
        "Uncaught exception, consider reporting the log (/tmp/MarzConverter.log).",
        exc_info=(exc_type, exc_value, exc_traceback),
    )


# Custom logging
sys.excepthook = handle_exception
